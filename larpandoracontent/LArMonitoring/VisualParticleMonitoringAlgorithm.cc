/**
 *  @file   larpandoracontent/LArMonitoring/VisualParticleMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/VisualParticleMonitoringAlgorithm.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"

using namespace pandora;

namespace lar_content
{

VisualParticleMonitoringAlgorithm::VisualParticleMonitoringAlgorithm() :
    m_groupByPdg(false),
    m_isTestBeam{false}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

VisualParticleMonitoringAlgorithm::~VisualParticleMonitoringAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VisualParticleMonitoringAlgorithm::Run()
{
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));

    LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
    this->MakeSelection(pMCParticleList, pCaloHitList, targetMCParticleToHitsMap);

    if (m_groupByPdg)
        this->VisualizeMCByPdgCode(targetMCParticleToHitsMap);
    else
        this->VisualizeIndependentMC(targetMCParticleToHitsMap);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VisualParticleMonitoringAlgorithm::VisualizeIndependentMC(const LArMCParticleHelper::MCContributionMap &mcMap)
{
    std::map<int, const std::string> keys = {{13, "mu"}, {11, "e"}, {22, "gamma"}, {321, "kaon"}, {211, "pi"}, {2212, "p"}};
    std::map<int, Color> colors = {{0, RED}, {1, BLACK}, {2, BLUE}, {3, CYAN}, {4, MAGENTA}, {5, GREEN}, {6, ORANGE}, {7, GRAY}};
    /*
    std::map<std::string, Color> colors = {{"mu", MAGENTA}, {"e", RED}, {"gamma", ORANGE}, {"kaon", BLACK}, {"pi", GREEN}, {"p", BLUE},
        {"other", GRAY}};
    */
    MCParticleList linearisedMC;
    if (mcMap.empty())
        return;

    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    LArMCParticleHelper::GetBreadthFirstHierarchyRepresentation(mcMap.begin()->first, linearisedMC);

    int colorIdx{0}, mcIdx{0};
    for (const MCParticle *pMC : linearisedMC)
    {
        const auto iter{mcMap.find(pMC)};
        if (iter == mcMap.end())
            continue;
        std::string key("other");
        try
        {
            const int pdg{std::abs(pMC->GetParticleId())};
            if (keys.find(pdg) != keys.end())
                key = keys[pdg];
        }
        catch (const StatusCodeException&)
        {
            key = "unknown";
        }

        CaloHitList uHits, vHits, wHits;
        for (const CaloHit *pCaloHit : iter->second)
        {
            const HitType view{pCaloHit->GetHitType()};
            if (view == HitType::TPC_VIEW_U)
                uHits.emplace_back(pCaloHit);
            else if (view == HitType::TPC_VIEW_V)
                vHits.emplace_back(pCaloHit);
            else
                wHits.emplace_back(pCaloHit);
        }
        std::string suffix{std::to_string(mcIdx) + "_" + key};
        if (!uHits.empty())
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &uHits, "u_" + suffix, colors[colorIdx]));
        if (!vHits.empty())
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &vHits, "v_" + suffix, colors[colorIdx]));
        if (!wHits.empty())
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &wHits, "w_" + suffix, colors[colorIdx]));
        colorIdx = colorIdx >= colors.size() ? 0 : colorIdx + 1;
        ++mcIdx;
    }

    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VisualParticleMonitoringAlgorithm::VisualizeMCByPdgCode(const LArMCParticleHelper::MCContributionMap &mcMap)
{
    std::map<int, const std::string> keys = {{13, "mu"}, {11, "e"}, {22, "gamma"}, {321, "kaon"}, {211, "pi"}, {2212, "p"}};
    std::map<std::string, Color> colors = {{"mu", MAGENTA}, {"e", RED}, {"gamma", ORANGE}, {"kaon", BLACK}, {"pi", GREEN}, {"p", BLUE},
        {"other", GRAY}};

    std::map<std::string, CaloHitList> uHits, vHits, wHits;
    for (const auto [ key, value ] : keys)
    {
        uHits[value] = CaloHitList();
        vHits[value] = CaloHitList();
        wHits[value] = CaloHitList();
    }
    uHits["other"] = CaloHitList();
    vHits["other"] = CaloHitList();
    wHits["other"] = CaloHitList();

    for (const auto [ pMC, pCaloHits ] : mcMap)
    {
        for (const CaloHit *pCaloHit : pCaloHits)
        {
            const HitType view{pCaloHit->GetHitType()};

            try
            {
                const int pdg{std::abs(pMC->GetParticleId())};
                std::string key("other");
                if (keys.find(pdg) != keys.end())
                    key = keys[pdg];

                if (view == HitType::TPC_VIEW_U)
                    uHits[key].emplace_back(pCaloHit);
                else if (view == HitType::TPC_VIEW_V)
                    vHits[key].emplace_back(pCaloHit);
                else
                    wHits[key].emplace_back(pCaloHit);
            }
            catch (const StatusCodeException&)
            {
                continue;
            }
        }
    }

    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));

    for (const auto [ key, value ] : keys)
        if (!uHits[value].empty())
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &uHits[value], "u_" + value, colors[value]));
    if (!uHits["other"].empty())
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &uHits["other"], "u_other", colors["other"]));

    for (const auto [ key, value ] : keys)
        if (!vHits[value].empty())
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &vHits[value], "v_" + value, colors[value]));
    if (!vHits["other"].empty())
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &uHits["other"], "v_other", colors["other"]));

    for (const auto [ key, value ] : keys)
        if (!wHits[value].empty())
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &wHits[value], "w_" + value, colors[value]));
    if (!wHits["other"].empty())
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &uHits["other"], "w_other", colors["other"]));

    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VisualParticleMonitoringAlgorithm::MakeSelection(const MCParticleList *pMCList, const CaloHitList *pCaloHitList,
    LArMCParticleHelper::MCContributionMap &mcMap)
{
    // Default reconstructability criteria are very liberal to allow for unfolded hierarchy
    LArMCParticleHelper::PrimaryParameters parameters;
    parameters.m_minPrimaryGoodHits = 2;
    parameters.m_minHitsForGoodView = 1;
    parameters.m_maxPhotonPropagation = std::numeric_limits<float>::max();
    parameters.m_minHitSharingFraction = 0;
    parameters.m_foldBackHierarchy = false;

    if (!m_isTestBeam)
    {
        LArMCParticleHelper::SelectReconstructableMCParticles(pMCList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, mcMap);
        LArMCParticleHelper::SelectReconstructableMCParticles(pMCList, pCaloHitList, parameters, LArMCParticleHelper::IsCosmicRay, mcMap);
    }
    else
    {
        LArMCParticleHelper::SelectReconstructableMCParticles(pMCList, pCaloHitList, parameters, LArMCParticleHelper::IsBeamParticle, mcMap);
        LArMCParticleHelper::SelectReconstructableMCParticles(pMCList, pCaloHitList, parameters, LArMCParticleHelper::IsCosmicRay, mcMap);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VisualParticleMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    if (m_caloHitListName.empty())
        m_caloHitListName = "CaloHitList2D";

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "GroupMCByPDG", m_groupByPdg));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

