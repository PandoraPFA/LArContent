/**
 *  @file   larpandoracontent/LArCheating/ValidationAlgorithm.cc
 *
 *  @brief  Implementation of the validation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArHierarchyHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArMetrics/ValidationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ValidationAlgorithm::ValidationAlgorithm() :
    m_caloHitListName("CaloHitList2D"),
    m_mcParticleListName("Input"),
    m_pfoListNames({"TrackParticles3D", "ShowerParticles3D", "NeutrinoParticles3D"}),
    m_fileName("Validation.root"),
    m_minPurity(0.5f),
    m_minCompleteness(0.1f),
    m_minRecoHits(15),
    m_minRecoHitsPerView(5),
    m_minRecoGoodViews(2),
    m_removeRecoNeutrons(true),
    m_selectRecoHits(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

ValidationAlgorithm::~ValidationAlgorithm()
{
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "EventTree", m_fileName, "UPDATE"));
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "HierarchyTree", m_fileName, "UPDATE"));
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "PFPTree", m_fileName, "UPDATE"));
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "TrackTree", m_fileName, "UPDATE"));
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "ShowerTree", m_fileName, "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ValidationAlgorithm::Run()
{
    // Get Lists
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));
    PfoList pfoList;
    for (const std::string &pfoListName : m_pfoListNames)
    {
        const PfoList *pPfoList(nullptr);

        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, pfoListName, pPfoList))
            continue;

        pfoList.insert(pfoList.begin(), pPfoList->begin(), pPfoList->end());
    }
    if (pfoList.empty())
        return STATUS_CODE_SUCCESS;

    // Folding options - want to fold showers
    LArHierarchyHelper::FoldingParameters foldParameters;
    foldParameters.m_foldToLeadingShowers = true;

    // RecoCriteria - want to catch reconstructed 'non-targets'
    const LArHierarchyHelper::MCHierarchy::ReconstructabilityCriteria recoCriteria(1, 1, 1, m_removeRecoNeutrons);

    // Matching
    LArHierarchyHelper::MCHierarchy mcHierarchy(recoCriteria);
    LArHierarchyHelper::FillMCHierarchy(*pMCParticleList, *pCaloHitList, foldParameters, mcHierarchy);
    LArHierarchyHelper::RecoHierarchy recoHierarchy;
    LArHierarchyHelper::FillRecoHierarchy(pfoList, foldParameters, recoHierarchy);
    const LArHierarchyHelper::QualityCuts quality(m_minPurity, m_minCompleteness, m_selectRecoHits);
    LArHierarchyHelper::MatchInfo matchInfo(mcHierarchy, recoHierarchy, quality);
    LArHierarchyHelper::MatchHierarchies(matchInfo);

    // Get hierarchy parents
    MCParticleList nuParticles;
    mcHierarchy.GetRootMCParticles(nuParticles);
    MCParticleVector nuParticlesVec;
    nuParticlesVec.insert(nuParticlesVec.begin(), nuParticles.begin(), nuParticles.end());
    std::sort(nuParticlesVec.begin(), nuParticlesVec.end(), LArMCParticleHelper::SortByMomentum);

    // Input entry for each hierarchy
    for (unsigned int i = 0; i < nuParticlesVec.size(); ++i)
    {
        MCParticleVector targetMC; PfoVector bestRecoMatch; IntVector isTarget;
        // ATTN: matches are NOT sorted
        LArHierarchyHelper::MCMatchesVector mcMatchesVec(matchInfo.GetMatches(nuParticlesVec.at(i)));

        for (const LArHierarchyHelper::MCMatches &mcMatches : mcMatchesVec)
        {
            // Is MCParticle a target?
            const LArHierarchyHelper::MCHierarchy::Node *pMCNode(mcMatches.GetMC());
            const bool thisIsTarget(this->IsReconstructable(pMCNode));

            // ATTN: consider in metrics if target OR if reconstructed
            const int nMatches(mcMatches.GetRecoMatches().size());
            if (!thisIsTarget && !mcMatches.IsQuality(quality))
                continue;

            isTarget.push_back(thisIsTarget);
            targetMC.push_back(pMCNode->GetMCParticles().front());

            // Determine best match pfo (if it exists)
            if (nMatches == 0)
            {
                bestRecoMatch.push_back(nullptr);
            }
            else
            {
                // Find the best match
                const ParticleFlowObject *pBestMatch(nullptr);
                float bestCompleteness(0.0);
                for (const auto pRecoNode : mcMatches.GetRecoMatches())
                {
                    const float thisCompleteness(mcMatches.GetCompleteness(pRecoNode, false));
                    if (thisCompleteness > bestCompleteness)
                    {
                        bestCompleteness = thisCompleteness;
                        pBestMatch = pRecoNode->GetRecoParticles().front();
                    }
                }
                bestRecoMatch.push_back(pBestMatch);
            }
        }

        // Ignore if no targets
        if (targetMC.empty())
            continue;

        // Run tools
        for (BaseValidationTool *const pValidationTool : m_validationToolVector)
            pValidationTool->Run(this, nuParticlesVec.at(i), mcMatchesVec, targetMC, bestRecoMatch);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ValidationAlgorithm::IsReconstructable(const LArHierarchyHelper::MCHierarchy::Node *pMCNode)
{
    unsigned int nHitsU(0), nHitsV(0), nHitsW(0);
    const CaloHitList caloHitList(pMCNode->GetCaloHits());

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        if (pCaloHit->GetHitType() == TPC_VIEW_U)
        {
            ++nHitsU;
        }
        else if (pCaloHit->GetHitType() == TPC_VIEW_V)
        {
            ++nHitsV;
        }
        else if (pCaloHit->GetHitType() == TPC_VIEW_W)
        {
            ++nHitsW;
        }
    }

    if ((nHitsU + nHitsV + nHitsW) < m_minRecoHits)
        return false;

    unsigned int nAboveThresholdViews(0);
    for (unsigned int nHits : {nHitsU, nHitsV, nHitsW})
    {
        if (nHits >= m_minRecoHitsPerView)
            ++nAboveThresholdViews;
    }

    return (nAboveThresholdViews >= m_minRecoGoodViews);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "PfoListNames", m_pfoListNames));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FileName", m_fileName));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinPurity", m_minPurity));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinCompleteness", m_minCompleteness));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinRecoHits", m_minRecoHits));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinRecoHitsPerView", m_minRecoHitsPerView));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinRecoGoodViews", m_minRecoGoodViews));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "RemoveRecoNeutrons", m_removeRecoNeutrons));

    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(
        STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "ValidationTools", algorithmToolVector));

    for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
    {
        BaseValidationTool *const pValidationTool(dynamic_cast<BaseValidationTool *>(pAlgorithmTool));

        if (!pValidationTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_validationToolVector.push_back(pValidationTool);
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
