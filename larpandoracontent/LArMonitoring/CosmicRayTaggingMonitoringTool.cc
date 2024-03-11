/**
 *  @file   larpandoracontent/LArMonitoring/CosmicRayTaggingMonitoringTool.cc
 *
 *  @brief  Implementation of the cosmic-ray tagging monitoring tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/PdgTable.h"

#include "larpandoracontent/LArMonitoring/CosmicRayTaggingMonitoringTool.h"

#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"

using namespace pandora;

namespace lar_content
{

CosmicRayTaggingMonitoringTool::CosmicRayTaggingMonitoringTool() :
    m_minHitsToConsiderTagging(15),
    m_minPurity(0.95),
    m_minImpurity(0.95),
    m_minSignificance(0.1)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTaggingMonitoringTool::FindAmbiguousPfos(const PfoList &parentCosmicRayPfos, PfoList &ambiguousPfos, const MasterAlgorithm *const pAlgorithm)
{
    if (this->GetPandora().GetSettings()->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    const MCParticleList *pMCParticleList = nullptr;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*pAlgorithm, pMCParticleList));

    const CaloHitList *pCaloHitList = nullptr;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*pAlgorithm, m_caloHitList2D, pCaloHitList));

    // Identify reconstructable MCParticles, and get mappings to their good hits
    LArMCParticleHelper::MCContributionMap nuMCParticlesToGoodHitsMap;
    LArMCParticleHelper::MCContributionMap beamMCParticlesToGoodHitsMap;
    LArMCParticleHelper::MCContributionMap crMCParticlesToGoodHitsMap;

    LArMCParticleHelper::SelectReconstructableMCParticles(
        pMCParticleList, pCaloHitList, m_parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, nuMCParticlesToGoodHitsMap);
    LArMCParticleHelper::SelectReconstructableMCParticles(
        pMCParticleList, pCaloHitList, m_parameters, LArMCParticleHelper::IsBeamParticle, beamMCParticlesToGoodHitsMap);
    LArMCParticleHelper::SelectReconstructableMCParticles(
        pMCParticleList, pCaloHitList, m_parameters, LArMCParticleHelper::IsCosmicRay, crMCParticlesToGoodHitsMap);

    // Get the hit sharing maps between Pfos and reconstructable MCParticles
    LArMCParticleHelper::MCContributionMapVector mcParticlesToGoodHitsMaps(
        {nuMCParticlesToGoodHitsMap, beamMCParticlesToGoodHitsMap, crMCParticlesToGoodHitsMap});

    LArMCParticleHelper::PfoContributionMap pfoToReconstructable2DHitsMap;
    LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(
        parentCosmicRayPfos, mcParticlesToGoodHitsMaps, pfoToReconstructable2DHitsMap, m_parameters.m_foldBackHierarchy);

    LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCParticleHitSharingMap;
    LArMCParticleHelper::MCParticleToPfoHitSharingMap mcParticleToPfoHitSharingMap;
    LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(
        pfoToReconstructable2DHitsMap, mcParticlesToGoodHitsMaps, pfoToMCParticleHitSharingMap, mcParticleToPfoHitSharingMap);

    // Calculate the purity and significane and classification of each Pfo
    PfoToFloatMap pfoSignificanceMap;
    PfoToFloatMap pfoPurityMap;
    PfoClassificationMap pfoClassificationMap;
    LArMCParticleHelper::MCContributionMapVector targetsToGoodHitsMaps({nuMCParticlesToGoodHitsMap, beamMCParticlesToGoodHitsMap});
    this->CalculatePfoMetrics(pfoToMCParticleHitSharingMap, pfoToReconstructable2DHitsMap, targetsToGoodHitsMaps, pfoSignificanceMap,
        pfoPurityMap, pfoClassificationMap);

    // -------------------------------------------------------------------------------------------------------------------------------------

    // Print the monte-carlo information for this event
    MCParticleVector orderedMCParticleVector;
    LArMonitoringHelper::GetOrderedMCParticleVector(mcParticlesToGoodHitsMaps, orderedMCParticleVector);

    LArFormattingHelper::PrintHeader("MC : Reconstructable neutrino final state particles");
    LArMonitoringHelper::PrintMCParticleTable(nuMCParticlesToGoodHitsMap, orderedMCParticleVector);

    LArFormattingHelper::PrintHeader("MC : Reconstructable primary beam particles");
    LArMonitoringHelper::PrintMCParticleTable(beamMCParticlesToGoodHitsMap, orderedMCParticleVector);

    LArFormattingHelper::PrintHeader("MC : Reconstructable primary cosmic-rays");
    LArMonitoringHelper::PrintMCParticleTable(crMCParticlesToGoodHitsMap, orderedMCParticleVector);

    LArFormattingHelper::PrintHeader("Reco : Primary cosmic-ray candidates");
    std::cout << "Columns with headers [n] are the number of shared hits between the Pfo and the target MCParticle with ID n." << std::endl;

    PfoVector orderedPfoVector;
    LArMonitoringHelper::GetOrderedPfoVector(pfoToReconstructable2DHitsMap, orderedPfoVector);
    this->PrintPfoTable(orderedPfoVector, pfoToReconstructable2DHitsMap, pfoPurityMap, pfoSignificanceMap, pfoClassificationMap, ambiguousPfos);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTaggingMonitoringTool::CalculatePfoMetrics(const LArMCParticleHelper::PfoToMCParticleHitSharingMap &hitSharingMap,
    const LArMCParticleHelper::PfoContributionMap &pfoToCaloHitListMap, const LArMCParticleHelper::MCContributionMapVector &targetsToGoodHitsMaps,
    PfoToFloatMap &pfoSignificanceMap, PfoToFloatMap &pfoPurityMap, PfoClassificationMap &pfoClassificationMap) const
{
    PfoVector sortedPfos;
    for (const auto &mapEntry : hitSharingMap)
        sortedPfos.push_back(mapEntry.first);
    std::sort(sortedPfos.begin(), sortedPfos.end(), LArPfoHelper::SortByNHits);

    for (const ParticleFlowObject *const pPfo : sortedPfos)
    {
        if (pfoToCaloHitListMap.find(pPfo) == pfoToCaloHitListMap.end())
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);

        const unsigned int n2DHits(pfoToCaloHitListMap.at(pPfo).size());
        float significance(0);
        float purity(0);

        if (n2DHits != 0)
        {
            // Sum over all target/shared hits pairs
            for (const LArMCParticleHelper::MCParticleCaloHitListPair &targetHitsShared : hitSharingMap.at(pPfo))
            {
                bool foundTarget(false);
                unsigned int nMCHits(std::numeric_limits<unsigned int>::max());

                // ATTN This map is unordered, but this does not impact search for specific target hit
                for (const LArMCParticleHelper::MCContributionMap &mcContributionMap : targetsToGoodHitsMaps)
                {
                    if (mcContributionMap.find(targetHitsShared.first) != mcContributionMap.end())
                    {
                        foundTarget = true;
                        nMCHits = mcContributionMap.at(targetHitsShared.first).size();
                        break;
                    }
                }

                if (!foundTarget)
                    continue;

                significance += static_cast<float>(targetHitsShared.second.size()) / static_cast<float>(nMCHits);
                purity += static_cast<float>(targetHitsShared.second.size()) / static_cast<float>(n2DHits);
            }
        }

        if (!pfoSignificanceMap.insert(PfoToFloatMap::value_type(pPfo, significance)).second)
            throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);

        if (!pfoPurityMap.insert(PfoToFloatMap::value_type(pPfo, purity)).second)
            throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);

        Classification classification(this->ClassifyPfo(n2DHits, significance, purity, this->IsMainMCParticleMuon(pPfo)));
        if (!pfoClassificationMap.insert(PfoClassificationMap::value_type(pPfo, classification)).second)
            throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CosmicRayTaggingMonitoringTool::IsMainMCParticleMuon(const ParticleFlowObject *const /*pPfo*/) const
{
    bool isMuon(false);
    try
    {
        isMuon = false; // TODO Local treatment is being developed, specific to this tool (std::abs(LArMCParticleHelper::GetMainMCParticle(pPfo)->GetParticleId()) == MU_MINUS);
    }
    catch (const StatusCodeException &)
    {
    }

    return isMuon;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CosmicRayTaggingMonitoringTool::Classification CosmicRayTaggingMonitoringTool::ClassifyPfo(
    const unsigned int &nHits, const float &significance, const float &purity, const bool isMuon) const
{
    if (nHits < m_minHitsToConsiderTagging)
        return SPARSE;

    bool isPure(purity > m_minPurity);
    bool isImpure((1 - purity) > m_minImpurity);
    bool isSignificant(significance > m_minSignificance);

    if (!isPure && !isImpure)
        return MIXED;

    if (isPure && isSignificant)
        return TARGET;

    if (isPure && !isSignificant)
        return FRAGMENTED;

    if (!isPure && isSignificant)
        return ABSORBED;

    if (!isPure && !isSignificant && isMuon)
        return CR_MUON;

    // !isPure && !isSignificant && !isMuon
    return CR_OTHER;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTaggingMonitoringTool::PrintPfoTable(const PfoVector &orderedPfoVector,
    const LArMCParticleHelper::PfoContributionMap &pfoToReconstructable2DHitsMap, const PfoToFloatMap &pfoPurityMap,
    const PfoToFloatMap &pfoSignificanceMap, const PfoClassificationMap &pfoClassificationMap, const PfoList &ambiguousPfos) const
{
    if (orderedPfoVector.empty())
    {
        std::cout << "No Pfos supplied." << std::endl;
        return;
    }

    LArFormattingHelper::Table table(
        {"ID", "PID", "", "nHits", "U", "V", "W", "", "nGoodHits", "U", "V", "W", "", "Purity", "Significance", "Classification", "", "Tagged?"});

    for (unsigned int id = 0; id < orderedPfoVector.size(); ++id)
    {
        const ParticleFlowObject *const pPfo(orderedPfoVector.at(id));

        LArMCParticleHelper::PfoContributionMap::const_iterator it = pfoToReconstructable2DHitsMap.find(pPfo);
        if (pfoToReconstructable2DHitsMap.end() == it)
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);

        if (pfoPurityMap.end() == pfoPurityMap.find(pPfo))
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);

        if (pfoSignificanceMap.end() == pfoSignificanceMap.find(pPfo))
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);

        if (pfoClassificationMap.end() == pfoClassificationMap.find(pPfo))
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);

        table.AddElement(id);
        table.AddElement(pPfo->GetParticleId());

        CaloHitList all2DCaloHits;
        LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_U, all2DCaloHits);
        LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_V, all2DCaloHits);
        LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, all2DCaloHits);

        table.AddElement(all2DCaloHits.size());
        table.AddElement(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, all2DCaloHits));
        table.AddElement(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, all2DCaloHits));
        table.AddElement(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, all2DCaloHits));

        table.AddElement(it->second.size());
        table.AddElement(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, it->second));
        table.AddElement(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, it->second));
        table.AddElement(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, it->second));

        table.AddElement(pfoPurityMap.at(pPfo));
        table.AddElement(pfoSignificanceMap.at(pPfo));

        const Classification classification(pfoClassificationMap.at(pPfo));
        table.AddElement(this->GetClassificationName(classification), LArFormattingHelper::INVERTED, this->GetClassificationColor(classification));

        const bool isTagged(std::find(ambiguousPfos.begin(), ambiguousPfos.end(), pPfo) == ambiguousPfos.end());
        const bool isGoodTag(isTagged && (classification == CR_MUON || classification == CR_OTHER));
        const bool isBadTag(isTagged && (classification == TARGET));
        const LArFormattingHelper::Color tagColor(
            isGoodTag ? LArFormattingHelper::LIGHT_GREEN : (isBadTag ? LArFormattingHelper::LIGHT_RED : LArFormattingHelper::LIGHT_YELLOW));
        table.AddElement(isTagged ? "yes" : "no", LArFormattingHelper::INVERTED, tagColor);
    }

    table.Print();
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArFormattingHelper::Color CosmicRayTaggingMonitoringTool::GetClassificationColor(const Classification &classification) const
{
    switch (classification)
    {
        case TARGET:
            return LArFormattingHelper::LIGHT_GREEN;
        case CR_MUON:
            return LArFormattingHelper::LIGHT_RED;
        case CR_OTHER:
            return LArFormattingHelper::LIGHT_YELLOW;
        case FRAGMENTED:
            return LArFormattingHelper::LIGHT_BLUE;
        case ABSORBED:
            return LArFormattingHelper::LIGHT_MAGENTA;
        case MIXED:
            return LArFormattingHelper::LIGHT_CYAN;
        default:
            return LArFormattingHelper::LIGHT_GRAY;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::string CosmicRayTaggingMonitoringTool::GetClassificationName(const Classification &classification) const
{
    switch (classification)
    {
        case TARGET:
            return "TARGET";
        case CR_MUON:
            return "CR_MUON";
        case CR_OTHER:
            return "CR_OTHER";
        case FRAGMENTED:
            return "FRAGMENTED";
        case ABSORBED:
            return "ABSORBED";
        case MIXED:
            return "MIXED";
        case SPARSE:
            return "SPARSE";
        default:
            return "UNCLASSIFIED";
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayTaggingMonitoringTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitList2D", m_caloHitList2D));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinPrimaryGoodHits", m_parameters.m_minPrimaryGoodHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinHitsForGoodView", m_parameters.m_minHitsForGoodView));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinPrimaryGoodViews", m_parameters.m_minPrimaryGoodViews));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SelectInputHits", m_parameters.m_selectInputHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxPhotonPropagation", m_parameters.m_maxPhotonPropagation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinHitSharingFraction", m_parameters.m_minHitSharingFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinHitsToConsiderTagging", m_minHitsToConsiderTagging));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinPurity", m_minPurity));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinImpurity", m_minImpurity));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinSignificance", m_minSignificance));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
