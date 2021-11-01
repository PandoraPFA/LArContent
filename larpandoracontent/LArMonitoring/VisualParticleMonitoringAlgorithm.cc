/**
 *  @file   larpandoracontent/LArMonitoring/VisualParticleMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/VisualParticleMonitoringAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"

using namespace pandora;

namespace lar_content
{

VisualParticleMonitoringAlgorithm::VisualParticleMonitoringAlgorithm() :
    m_visualizeMC(false),
    m_visualizePfo(false),
    m_groupMCByPdg(false),
    m_showPfoByPid(false),
    m_showPfoMatchedMC(false),
    m_isTestBeam{false},
    m_transparencyThresholdE{-1.f},
    m_energyScaleThresholdE{1.f},
    m_scalingFactor{1.f}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

VisualParticleMonitoringAlgorithm::~VisualParticleMonitoringAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VisualParticleMonitoringAlgorithm::Run()
{
#ifdef MONITORING
    LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
    if (m_visualizeMC || m_showPfoMatchedMC)
    {
        const CaloHitList *pCaloHitList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
        const MCParticleList *pMCParticleList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));
        this->MakeSelection(pMCParticleList, pCaloHitList, targetMCParticleToHitsMap);
    }

    if (m_visualizeMC)
    {
        if (m_groupMCByPdg)
            this->VisualizeMCByPdgCode(targetMCParticleToHitsMap);
        else
            this->VisualizeIndependentMC(targetMCParticleToHitsMap);
    }
    if (m_visualizePfo)
    {
        const PfoList *pPfoList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_pfoListName, pPfoList));
        if (m_showPfoByPid)
        {
            this->VisualizePfoByParticleId(*pPfoList);
        }
        else
        {
            if (m_showPfoMatchedMC)
                this->VisualizeIndependentPfo(*pPfoList, targetMCParticleToHitsMap);
            else
                this->VisualizeIndependentPfo(*pPfoList);
        }
    }
#endif // MONITORING
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

#ifdef MONITORING

void VisualParticleMonitoringAlgorithm::VisualizeIndependentMC(const LArMCParticleHelper::MCContributionMap &mcMap) const
{
    const std::map<int, const std::string> keys = {{13, "mu"}, {11, "e"}, {22, "gamma"}, {321, "kaon"}, {211, "pi"}, {2212, "p"}};
    const std::map<int, Color> colors = {{0, RED}, {1, BLACK}, {2, BLUE}, {3, CYAN}, {4, MAGENTA}, {5, GREEN}, {6, ORANGE}, {7, GRAY}};
    MCParticleList linearisedMC;
    if (mcMap.empty())
        return;

    PANDORA_MONITORING_API(
        SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, m_transparencyThresholdE, m_energyScaleThresholdE, m_scalingFactor));
    LArMCParticleHelper::GetBreadthFirstHierarchyRepresentation(mcMap.begin()->first, linearisedMC);

    size_t colorIdx{0};
    int mcIdx{0};
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
                key = keys.at(pdg);
        }
        catch (const StatusCodeException &)
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
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &uHits, "u_" + suffix, colors.at(colorIdx)));
        if (!vHits.empty())
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &vHits, "v_" + suffix, colors.at(colorIdx)));
        if (!wHits.empty())
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &wHits, "w_" + suffix, colors.at(colorIdx)));
        colorIdx = (colorIdx + 1) >= colors.size() ? 0 : colorIdx + 1;
        ++mcIdx;
    }

    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VisualParticleMonitoringAlgorithm::VisualizeMCByPdgCode(const LArMCParticleHelper::MCContributionMap &mcMap) const
{
    const std::map<int, const std::string> keys = {{13, "mu"}, {11, "e"}, {22, "gamma"}, {321, "kaon"}, {211, "pi"}, {2212, "p"}};
    const std::map<std::string, Color> colors = {
        {"mu", MAGENTA}, {"e", RED}, {"gamma", ORANGE}, {"kaon", BLACK}, {"pi", GREEN}, {"p", BLUE}, {"other", GRAY}};

    std::map<std::string, CaloHitList> uHits, vHits, wHits;
    for (const auto & [key, value] : keys)
    {
        (void)key; // GCC 7 support, 8+ doesn't need this
        uHits[value] = CaloHitList();
        vHits[value] = CaloHitList();
        wHits[value] = CaloHitList();
    }
    uHits["other"] = CaloHitList();
    vHits["other"] = CaloHitList();
    wHits["other"] = CaloHitList();

    for (const auto & [pMC, pCaloHits] : mcMap)
    {
        for (const CaloHit *pCaloHit : pCaloHits)
        {
            const HitType view{pCaloHit->GetHitType()};

            try
            {
                const int pdg{std::abs(pMC->GetParticleId())};
                std::string key("other");
                if (keys.find(pdg) != keys.end())
                    key = keys.at(pdg);

                if (view == HitType::TPC_VIEW_U)
                    uHits[key].emplace_back(pCaloHit);
                else if (view == HitType::TPC_VIEW_V)
                    vHits[key].emplace_back(pCaloHit);
                else
                    wHits[key].emplace_back(pCaloHit);
            }
            catch (const StatusCodeException &)
            {
                continue;
            }
        }
    }

    PANDORA_MONITORING_API(
        SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, m_transparencyThresholdE, m_energyScaleThresholdE, m_scalingFactor));

    for (const auto & [key, value] : keys)
    {
        (void)key; // GCC 7 support, 8+ doesn't need this
        if (!uHits[value].empty())
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &uHits[value], "u_" + value, colors.at(value)));
    }
    if (!uHits["other"].empty())
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &uHits["other"], "u_other", colors.at("other")));

    for (const auto & [key, value] : keys)
    {
        (void)key; // GCC 7 support, 8+ doesn't need this
        if (!vHits[value].empty())
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &vHits[value], "v_" + value, colors.at(value)));
    }
    if (!vHits["other"].empty())
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &uHits["other"], "v_other", colors.at("other")));

    for (const auto & [key, value] : keys)
    {
        (void)key; // GCC 7 support, 8+ doesn't need this
        if (!wHits[value].empty())
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &wHits[value], "w_" + value, colors.at(value)));
    }
    if (!wHits["other"].empty())
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &uHits["other"], "w_other", colors.at("other")));

    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VisualParticleMonitoringAlgorithm::VisualizeIndependentPfo(const PfoList &pfoList) const
{
    // ATTN - If we aren't showing matched MC, just pass in an empty MC to hits map
    LArMCParticleHelper::MCContributionMap mcMap;
    this->VisualizeIndependentPfo(pfoList, mcMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VisualParticleMonitoringAlgorithm::VisualizeIndependentPfo(const PfoList &pfoList, const LArMCParticleHelper::MCContributionMap &mcMap) const
{
    const std::map<int, Color> colors = {{0, RED}, {1, BLACK}, {2, BLUE}, {3, CYAN}, {4, MAGENTA}, {5, GREEN}, {6, ORANGE}, {7, GRAY}};
    PfoList linearisedPfo;
    if (pfoList.empty())
        return;

    PANDORA_MONITORING_API(
        SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, m_transparencyThresholdE, m_energyScaleThresholdE, m_scalingFactor));
    LArPfoHelper::GetBreadthFirstHierarchyRepresentation(pfoList.front(), linearisedPfo);

    size_t colorIdx{0};
    int pfoIdx{0};
    for (const ParticleFlowObject *pPfo : linearisedPfo)
    {
        CaloHitList uHits, vHits, wHits;
        const bool isTrack{LArPfoHelper::IsTrack(pPfo)};
        CaloHitList caloHits;
        for (const auto view : {HitType::TPC_VIEW_U, HitType::TPC_VIEW_V, HitType::TPC_VIEW_W})
        {
            LArPfoHelper::GetCaloHits(pPfo, view, caloHits);
            LArPfoHelper::GetIsolatedCaloHits(pPfo, view, caloHits);
        }
        for (const CaloHit *pCaloHit : caloHits)
        {
            const HitType view{pCaloHit->GetHitType()};
            if (view == HitType::TPC_VIEW_U)
                uHits.emplace_back(pCaloHit);
            else if (view == HitType::TPC_VIEW_V)
                vHits.emplace_back(pCaloHit);
            else if (view == HitType::TPC_VIEW_W)
                wHits.emplace_back(pCaloHit);
        }
        std::string suffix{std::to_string(pfoIdx)};
        suffix += isTrack ? "_T" : "_S";
        if (!uHits.empty())
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &uHits, "u_" + suffix, colors.at(colorIdx)));
        if (!vHits.empty())
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &vHits, "v_" + suffix, colors.at(colorIdx)));
        if (!wHits.empty())
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &wHits, "w_" + suffix, colors.at(colorIdx)));
        colorIdx = (colorIdx + 1) >= colors.size() ? 0 : colorIdx + 1;
        if (m_showPfoMatchedMC)
        {
            try
            {
                const MCParticle *pMC{LArMCParticleHelper::GetMainMCParticle(pPfo)};
                if (pMC)
                {
                    const auto iter{mcMap.find(pMC)};
                    if (iter != mcMap.end())
                    {
                        CaloHitList uHitsMC, vHitsMC, wHitsMC;
                        for (const CaloHit *pCaloHit : iter->second)
                        {
                            const HitType view{pCaloHit->GetHitType()};
                            if (view == HitType::TPC_VIEW_U)
                                uHitsMC.emplace_back(pCaloHit);
                            else if (view == HitType::TPC_VIEW_V)
                                vHitsMC.emplace_back(pCaloHit);
                            else if (view == HitType::TPC_VIEW_W)
                                wHitsMC.emplace_back(pCaloHit);
                        }
                        std::string mcSuffix(suffix + "_MC");
                        if (!uHitsMC.empty())
                            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &uHitsMC, "u_" + mcSuffix, colors.at(colorIdx)));
                        if (!vHitsMC.empty())
                            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &vHitsMC, "v_" + mcSuffix, colors.at(colorIdx)));
                        if (!wHitsMC.empty())
                            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &wHitsMC, "w_" + mcSuffix, colors.at(colorIdx)));
                    }
                }
            }
            catch (const StatusCodeException &)
            { // No matched MC, move on
            }
        }
        ++pfoIdx;
    }

    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VisualParticleMonitoringAlgorithm::VisualizePfoByParticleId(const PfoList &pfoList) const
{
    PfoList linearisedPfo;
    if (pfoList.empty())
        return;

    PANDORA_MONITORING_API(
        SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, m_transparencyThresholdE, m_energyScaleThresholdE, m_scalingFactor));
    LArPfoHelper::GetBreadthFirstHierarchyRepresentation(pfoList.front(), linearisedPfo);

    int pfoIdx{0};
    for (const ParticleFlowObject *pPfo : linearisedPfo)
    {
        CaloHitList uTrackHits, vTrackHits, wTrackHits, uShowerHits, vShowerHits, wShowerHits;
        const bool isTrack{LArPfoHelper::IsTrack(pPfo)};
        CaloHitList caloHits;
        for (const auto view : {HitType::TPC_VIEW_U, HitType::TPC_VIEW_V, HitType::TPC_VIEW_W})
        {
            LArPfoHelper::GetCaloHits(pPfo, view, caloHits);
            LArPfoHelper::GetIsolatedCaloHits(pPfo, view, caloHits);
        }
        for (const CaloHit *pCaloHit : caloHits)
        {
            const HitType view{pCaloHit->GetHitType()};
            if (view == HitType::TPC_VIEW_U)
            {
                if (isTrack)
                    uTrackHits.emplace_back(pCaloHit);
                else
                    uShowerHits.emplace_back(pCaloHit);
            }
            else if (view == HitType::TPC_VIEW_V)
            {
                if (isTrack)
                    vTrackHits.emplace_back(pCaloHit);
                else
                    vShowerHits.emplace_back(pCaloHit);
            }
            else
            {
                if (isTrack)
                    wTrackHits.emplace_back(pCaloHit);
                else
                    wShowerHits.emplace_back(pCaloHit);
            }
        }
        if (isTrack)
        {
            std::string suffix{std::to_string(pfoIdx) + "_T"};
            if (!uTrackHits.empty())
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &uTrackHits, "u_" + suffix, BLUE));
            if (!vTrackHits.empty())
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &vTrackHits, "v_" + suffix, BLUE));
            if (!wTrackHits.empty())
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &wTrackHits, "w_" + suffix, BLUE));
        }
        else
        {
            std::string suffix{std::to_string(pfoIdx) + "_S"};
            if (!uShowerHits.empty())
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &uShowerHits, "u_" + suffix, RED));
            if (!vShowerHits.empty())
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &vShowerHits, "v_" + suffix, RED));
            if (!wShowerHits.empty())
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &wShowerHits, "w_" + suffix, RED));
        }
        ++pfoIdx;
    }

    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VisualParticleMonitoringAlgorithm::MakeSelection(
    const MCParticleList *pMCList, const CaloHitList *pCaloHitList, LArMCParticleHelper::MCContributionMap &mcMap) const
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
#endif // MONITORING

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VisualParticleMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    if (m_caloHitListName.empty())
        m_caloHitListName = "CaloHitList2D";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    if (m_pfoListName.empty())
        m_pfoListName = "RecreatedPfos";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VisualizeMC", m_visualizeMC));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VisualizePFO", m_visualizePfo));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "GroupMCByPDG", m_groupMCByPdg));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ShowPFOByPID", m_showPfoByPid));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ShowPFOMatchedMC", m_showPfoMatchedMC));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TransparencyThresholdE", m_transparencyThresholdE));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "EnergyScaleThresholdE", m_energyScaleThresholdE));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ScalingFactor", m_scalingFactor));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
