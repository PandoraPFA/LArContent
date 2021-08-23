/**
 *  @file   larpandoracontent/LArMonitoring/HierarchyMonitoringAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/HierarchyMonitoringAlgorithm.h"

using namespace pandora;

namespace lar_content
{

HierarchyMonitoringAlgorithm::HierarchyMonitoringAlgorithm() :
    m_visualizeMC(false),
    m_visualizeReco(false),
    m_visualizeDistinct(false),
    m_match(false),
    m_transparencyThresholdE{-1.f},
    m_energyScaleThresholdE{1.f},
    m_scalingFactor{1.f}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

HierarchyMonitoringAlgorithm::~HierarchyMonitoringAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HierarchyMonitoringAlgorithm::Run()
{
    PANDORA_MONITORING_API(
        SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, m_transparencyThresholdE, m_energyScaleThresholdE, m_scalingFactor));

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));
    const PfoList *pPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_pfoListName, pPfoList));

    LArHierarchyHelper::MCHierarchy mcHierarchy;
    LArHierarchyHelper::RecoHierarchy recoHierarchy;
    LArHierarchyHelper::FoldingParameters foldParameters(1);

    if (m_visualizeMC || m_match)
        LArHierarchyHelper::FillMCHierarchy(*pMCParticleList, *pCaloHitList, foldParameters, mcHierarchy);
    if (m_visualizeReco || m_match)
        LArHierarchyHelper::FillRecoHierarchy(*pPfoList, foldParameters, recoHierarchy);
    if (m_match)
    {
        LArHierarchyHelper::MatchInfo matchInfo;
        LArHierarchyHelper::MatchHierarchies(mcHierarchy, recoHierarchy, matchInfo);
        matchInfo.Print(mcHierarchy);
        this->VisualizeMatches(matchInfo);
    }
    else
    {
        if (m_visualizeMC)
        {
            if (m_visualizeDistinct)
                this->VisualizeMCDistinct(mcHierarchy);
            else
                this->VisualizeMC(mcHierarchy);
        }
        if (m_visualizeReco)
            this->VisualizeReco(recoHierarchy);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HierarchyMonitoringAlgorithm::VisualizeMC(const LArHierarchyHelper::MCHierarchy &hierarchy) const
{
    const std::map<int, const std::string> keys = {{13, "mu"}, {11, "e"}, {22, "gamma"}, {321, "kaon"}, {211, "pi"}, {2212, "p"}};
    const std::map<std::string, int> colors = {{"mu", 5}, {"e", 2}, {"gamma", 9}, {"kaon", 1}, {"pi", 3}, {"p", 4}, {"other", 14}};

    LArHierarchyHelper::MCHierarchy::NodeVector nodes;
    hierarchy.GetFlattenedNodes(nodes);

    int nodeIdx{0};
    for (const LArHierarchyHelper::MCHierarchy::Node *pNode : nodes)
    {
        std::string key("other");
        const int pdg{std::abs(pNode->GetParticleId())};
        if (keys.find(pdg) != keys.end())
            key = keys.at(pdg);

        CaloHitList uHits, vHits, wHits;
        this->FillHitLists(pNode->GetCaloHits(), uHits, vHits, wHits);
        std::string suffix{std::to_string(nodeIdx) + "_" + key};
        this->Visualize(uHits, "u_" + suffix, colors.at(key));
        this->Visualize(vHits, "v_" + suffix, colors.at(key));
        this->Visualize(wHits, "w_" + suffix, colors.at(key));
        ++nodeIdx;
    }

    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HierarchyMonitoringAlgorithm::VisualizeMCDistinct(const LArHierarchyHelper::MCHierarchy &hierarchy) const
{
    const std::map<int, const std::string> keys = {{13, "mu"}, {11, "e"}, {22, "gamma"}, {321, "kaon"}, {211, "pi"}, {2212, "p"}};
    const int nColours{9};
    const int colors[nColours] = {1, 2, 3, 4, 5, 6, 7, 8, 9};

    LArHierarchyHelper::MCHierarchy::NodeVector nodes;
    hierarchy.GetFlattenedNodes(nodes);

    int nodeIdx{0}, colorIdx{0};
    for (const LArHierarchyHelper::MCHierarchy::Node *pNode : nodes)
    {
        std::string key("other");
        const int pdg{std::abs(pNode->GetParticleId())};
        if (keys.find(pdg) != keys.end())
            key = keys.at(pdg);

        CaloHitList uHits, vHits, wHits;
        this->FillHitLists(pNode->GetCaloHits(), uHits, vHits, wHits);
        std::string suffix{std::to_string(nodeIdx) + "_" + key};
        this->Visualize(uHits, "u_" + suffix, colors[colorIdx]);
        this->Visualize(vHits, "v_" + suffix, colors[colorIdx]);
        this->Visualize(wHits, "w_" + suffix, colors[colorIdx]);
        colorIdx = (colorIdx + 1) >= nColours ? 0 : colorIdx + 1;
        ++nodeIdx;
    }

    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HierarchyMonitoringAlgorithm::VisualizeMCProcess(const LArHierarchyHelper::MCHierarchy &hierarchy) const
{
    const std::map<int, const std::string> keys = {{13, "mu"}, {11, "e"}, {22, "gamma"}, {321, "kaon"}, {211, "pi"}, {2212, "p"}};
    const int nColours{9};
    const int colors[nColours] = {1, 2, 3, 4, 5, 6, 7, 8, 9};

    LArHierarchyHelper::MCHierarchy::NodeVector nodes;
    hierarchy.GetFlattenedNodes(nodes);

    int nodeIdx{0}, colorIdx{0};
    for (const LArHierarchyHelper::MCHierarchy::Node *pNode : nodes)
    {
        std::string key("other");
        const int pdg{std::abs(pNode->GetParticleId())};
        if (keys.find(pdg) != keys.end())
            key = keys.at(pdg);

        CaloHitList uHits, vHits, wHits;
        this->FillHitLists(pNode->GetCaloHits(), uHits, vHits, wHits);
        std::string suffix{std::to_string(nodeIdx) + "_" + key};
        this->Visualize(uHits, "u_" + suffix, colors[colorIdx]);
        this->Visualize(vHits, "v_" + suffix, colors[colorIdx]);
        this->Visualize(wHits, "w_" + suffix, colors[colorIdx]);
        colorIdx = (colorIdx + 1) >= nColours ? 0 : colorIdx + 1;
        ++nodeIdx;
    }

    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HierarchyMonitoringAlgorithm::VisualizeReco(const LArHierarchyHelper::RecoHierarchy &hierarchy) const
{
    const int nColors{7};
    const int colors[nColors] = {5, 2, 9, 1, 3, 4, 14};

    LArHierarchyHelper::RecoHierarchy::NodeVector nodes;
    hierarchy.GetFlattenedNodes(nodes);

    int colorIdx{0};
    int pfoIdx{0};
    for (const LArHierarchyHelper::RecoHierarchy::Node *pNode : nodes)
    {
        CaloHitList uHits, vHits, wHits;
        this->FillHitLists(pNode->GetCaloHits(), uHits, vHits, wHits);
        const int pdg{pNode->GetParticleId()};
        const std::string key{pdg == MU_MINUS ? "T" : pdg == E_MINUS ? "S" : "?"};
        const std::string suffix{std::to_string(pfoIdx) + "_" + key};
        this->Visualize(uHits, "u_" + suffix, colors[colorIdx]);
        this->Visualize(vHits, "v_" + suffix, colors[colorIdx]);
        this->Visualize(wHits, "w_" + suffix, colors[colorIdx]);
        colorIdx = (colorIdx + 1) >= nColors ? 0 : colorIdx + 1;
        ++pfoIdx;
    }

    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HierarchyMonitoringAlgorithm::VisualizeMatches(const LArHierarchyHelper::MatchInfo &matchInfo) const
{
    std::map<const LArHierarchyHelper::MCHierarchy::Node *, int> mcIdxMap;
    int mcIdx{0};
    for (const LArHierarchyHelper::MCMatches &matches : matchInfo.GetGoodMatches())
        if (mcIdxMap.find(matches.GetMC()) == mcIdxMap.end())
            mcIdxMap[matches.GetMC()] = mcIdx++;
    for (const LArHierarchyHelper::MCMatches &matches : matchInfo.GetAboveThresholdMatches())
        if (mcIdxMap.find(matches.GetMC()) == mcIdxMap.end())
            mcIdxMap[matches.GetMC()] = mcIdx++;
    for (const LArHierarchyHelper::MCMatches &matches : matchInfo.GetSubThresholdMatches())
        if (mcIdxMap.find(matches.GetMC()) == mcIdxMap.end())
            mcIdxMap[matches.GetMC()] = mcIdx++;
    for (const LArHierarchyHelper::MCHierarchy::Node *pNode : matchInfo.GetUnmatchedMC())
        if (mcIdxMap.find(pNode) == mcIdxMap.end())
            mcIdxMap[pNode] = mcIdx++;

    for (const LArHierarchyHelper::MCMatches &matches : matchInfo.GetGoodMatches())
        this->VisualizeMatchedMC(matches, mcIdxMap.at(matches.GetMC()));
    for (const LArHierarchyHelper::MCMatches &matches : matchInfo.GetAboveThresholdMatches())
        this->VisualizeMatchedMC(matches, mcIdxMap.at(matches.GetMC()));
    for (const LArHierarchyHelper::MCMatches &matches : matchInfo.GetSubThresholdMatches())
        this->VisualizeMatchedMC(matches, mcIdxMap.at(matches.GetMC()));

    for (const LArHierarchyHelper::MCHierarchy::Node *pNode : matchInfo.GetUnmatchedMC())
        this->VisualizeUnmatchedMC(pNode, mcIdxMap.at(pNode));
    for (const LArHierarchyHelper::RecoHierarchy::Node *pNode : matchInfo.GetUnmatchedReco())
        this->VisualizeUnmatchedReco(pNode);

    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HierarchyMonitoringAlgorithm::VisualizeMatchedMC(const LArHierarchyHelper::MCMatches &matches, const int mcIdx) const
{
    const std::map<int, const std::string> keys = {{13, "mu"}, {11, "e"}, {22, "gamma"}, {321, "kaon"}, {211, "pi"}, {2212, "p"}};
    const int green{3};
    const std::map<int, int> colors = {{0, 5}, {1, 2}, {2, 9}, {3, 1}, {4, 4}};

    const LArHierarchyHelper::MCHierarchy::Node *pMCNode{matches.GetMC()};
    const int pdg{pMCNode->GetParticleId()};

    std::string key("other");
    if (keys.find(pdg) != keys.end())
        key = keys.at(pdg);

    CaloHitList mcUHits, mcVHits, mcWHits;
    this->FillHitLists(pMCNode->GetCaloHits(), mcUHits, mcVHits, mcWHits);
    std::string leading{pMCNode->IsLeadingLepton() ? "_leading" : ""};
    std::string suffix{std::to_string(mcIdx) + "_" + key + leading};
    this->Visualize(mcUHits, "u_" + suffix, green);
    this->Visualize(mcVHits, "v_" + suffix, green);
    this->Visualize(mcWHits, "w_" + suffix, green);

    int recoIdx{0};
    size_t colorIdx{0};
    for (const LArHierarchyHelper::RecoHierarchy::Node *pRecoNode : matches.GetRecoMatches())
    {
        CaloHitList recoUHits, recoVHits, recoWHits;
        this->FillHitLists(pRecoNode->GetCaloHits(), recoUHits, recoVHits, recoWHits);
        const int characterisation{pRecoNode->GetParticleId()};
        const std::string recoKey{characterisation == MU_MINUS ? "T" : characterisation == E_MINUS ? "S" : "?"};
        const std::string recoSuffix{std::to_string(mcIdx) + "->" + std::to_string(recoIdx) + "_" + recoKey};
        if (!recoUHits.empty())
        {
            const float purity{matches.GetPurity(pRecoNode, TPC_VIEW_U, true)};
            const float completeness{matches.GetCompleteness(pRecoNode, TPC_VIEW_U, true)};
            std::ostringstream metrics;
            metrics.precision(2);
            metrics << std::fixed << " p: " << purity << " c: " << completeness;
            PANDORA_MONITORING_API(
                VisualizeCaloHits(this->GetPandora(), &recoUHits, "u_" + recoSuffix + metrics.str(), static_cast<Color>(colors.at(colorIdx))));
        }
        if (!recoVHits.empty())
        {
            const float purity{matches.GetPurity(pRecoNode, TPC_VIEW_V, true)};
            const float completeness{matches.GetCompleteness(pRecoNode, TPC_VIEW_V, true)};
            std::ostringstream metrics;
            metrics.precision(2);
            metrics << std::fixed << " p: " << purity << " c: " << completeness;
            PANDORA_MONITORING_API(
                VisualizeCaloHits(this->GetPandora(), &recoVHits, "v_" + recoSuffix + metrics.str(), static_cast<Color>(colors.at(colorIdx))));
        }
        if (!recoWHits.empty())
        {
            const float purity{matches.GetPurity(pRecoNode, TPC_VIEW_W, true)};
            const float completeness{matches.GetCompleteness(pRecoNode, TPC_VIEW_W, true)};
            std::ostringstream metrics;
            metrics.precision(2);
            metrics << std::fixed << " p: " << purity << " c: " << completeness;
            PANDORA_MONITORING_API(
                VisualizeCaloHits(this->GetPandora(), &recoWHits, "w_" + recoSuffix + metrics.str(), static_cast<Color>(colors.at(colorIdx))));
        }
        colorIdx = (colorIdx + 1) >= colors.size() ? 0 : colorIdx + 1;
        ++recoIdx;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HierarchyMonitoringAlgorithm::VisualizeUnmatchedMC(const LArHierarchyHelper::MCHierarchy::Node *pNode, const int mcIdx) const
{
    const std::map<int, const std::string> keys = {{13, "mu"}, {11, "e"}, {22, "gamma"}, {321, "kaon"}, {211, "pi"}, {2212, "p"}};
    const int gray{14};

    const int pdg{pNode->GetParticleId()};

    std::string key("other");
    if (keys.find(pdg) != keys.end())
        key = keys.at(pdg);

    CaloHitList uHits, vHits, wHits;
    this->FillHitLists(pNode->GetCaloHits(), uHits, vHits, wHits);
    std::string suffix{std::to_string(mcIdx) + "_" + key + "_unmatched"};
    this->Visualize(uHits, "u_" + suffix, gray);
    this->Visualize(vHits, "v_" + suffix, gray);
    this->Visualize(wHits, "w_" + suffix, gray);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HierarchyMonitoringAlgorithm::VisualizeUnmatchedReco(const LArHierarchyHelper::RecoHierarchy::Node *pNode) const
{
    const int gray{14};

    CaloHitList uHits, vHits, wHits;
    this->FillHitLists(pNode->GetCaloHits(), uHits, vHits, wHits);
    const int pdg{pNode->GetParticleId()};
    const std::string key{pdg == MU_MINUS ? "T" : pdg == E_MINUS ? "S" : "?"};
    const std::string suffix{"unmatched_reco"};
    this->Visualize(uHits, "u_" + suffix, gray);
    this->Visualize(vHits, "v_" + suffix, gray);
    this->Visualize(wHits, "w_" + suffix, gray);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HierarchyMonitoringAlgorithm::Visualize(const CaloHitList &hits, const std::string &label, const int color) const
{
    if (!hits.empty() && hits.size() > 1)
    {
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &hits, label, static_cast<Color>(color)));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HierarchyMonitoringAlgorithm::FillHitLists(const CaloHitList &hits, CaloHitList &uHits, CaloHitList &vHits, CaloHitList &wHits) const
{
    for (const CaloHit *pCaloHit : hits)
    {
        const HitType view{pCaloHit->GetHitType()};
        if (view == HitType::TPC_VIEW_U)
            uHits.emplace_back(pCaloHit);
        else if (view == HitType::TPC_VIEW_V)
            vHits.emplace_back(pCaloHit);
        else if (view == HitType::TPC_VIEW_W)
            wHits.emplace_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HierarchyMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    if (m_caloHitListName.empty())
        m_caloHitListName = "CaloHitList2D";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    if (m_pfoListName.empty())
        m_pfoListName = "RecreatedPfos";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VisualizeMC", m_visualizeMC));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VisualizeReco", m_visualizeReco));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VisualizeDistinct", m_visualizeDistinct));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PerformMatching", m_match));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TransparencyThresholdE", m_transparencyThresholdE));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "EnergyScaleThresholdE", m_energyScaleThresholdE));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ScalingFactor", m_scalingFactor));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
