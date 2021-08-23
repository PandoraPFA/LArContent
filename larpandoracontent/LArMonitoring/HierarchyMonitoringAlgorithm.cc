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
    m_visualizeProcess{false},
    m_match(false),
    m_collectionOnly{false},
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

    LArHierarchyHelper::MCHierarchy::ReconstructabilityCriteria criteria(15, 5, 2, false);
    LArHierarchyHelper::MCHierarchy mcHierarchy(criteria);
    LArHierarchyHelper::RecoHierarchy recoHierarchy;
    LArHierarchyHelper::FoldingParameters foldParameters; //(1);

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
            else if (m_visualizeProcess)
                this->VisualizeMCProcess(mcHierarchy);
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
    const std::map<MCProcess, std::string> procToCategoryMap = {{MC_PROC_INCIDENT_NU, "invisible"}, {MC_PROC_UNKNOWN, "invisible"},
        {MC_PROC_PRIMARY, "primary"}, {MC_PROC_COMPT, "compt"}, {MC_PROC_PHOT, "phot"}, {MC_PROC_ANNIHIL, "annihil"}, {MC_PROC_E_IONI, "ioni"},
        {MC_PROC_E_BREM, "brem"}, {MC_PROC_CONV, "conv"}, {MC_PROC_MU_IONI, "ioni"}, {MC_PROC_MU_MINUS_CAPTURE_AT_REST, "capture"},
        {MC_PROC_NEUTRON_INELASTIC, "inelastic"}, {MC_PROC_N_CAPTURE, "capture"}, {MC_PROC_HAD_ELASTIC, "elastic"}, {MC_PROC_DECAY, "decay"},
        {MC_PROC_COULOMB_SCAT, "coulomb"}, {MC_PROC_MU_BREM, "brem"}, {MC_PROC_MU_PAIR_PROD, "pair_prod"}, {MC_PROC_PHOTON_INELASTIC, "inelastic"},
        {MC_PROC_HAD_IONI, "ioni"}, {MC_PROC_PROTON_INELASTIC, "inelastic"}, {MC_PROC_PI_PLUS_INELASTIC, "inelastic"},
        {MC_PROC_CHIPS_NUCLEAR_CAPTURE_AT_REST, "capture"}, {MC_PROC_PI_MINUS_INELASTIC, "inelastic"}, {MC_PROC_TRANSPORTATION, "transport"},
        {MC_PROC_RAYLEIGH, "rayleigh"}, {MC_PROC_HAD_BREM, "brem"}, {MC_PROC_HAD_PAIR_PROD, "pair_prod"}, {MC_PROC_ION_IONI, "ioni"},
        {MC_PROC_NEUTRON_KILLER, "kill"}, {MC_PROC_ION_INELASTIC, "inelastic"}, {MC_PROC_HE3_INELASTIC, "inelastic"},
        {MC_PROC_ALPHA_INELASTIC, "inelastic"}, {MC_PROC_ANTI_HE3_INELASTIC, "inelastic"}, {MC_PROC_ANTI_ALPHA_INELASTIC, "inelastic"},
        {MC_PROC_HAD_FRITIOF_CAPTURE_AT_REST, "inelastic"}, {MC_PROC_ANTI_DEUTERON_INELASTIC, "inelastic"},
        {MC_PROC_ANTI_NEUTRON_INELASTIC, "inelastic"}, {MC_PROC_ANTI_PROTON_INELASTIC, "inelastic"},
        {MC_PROC_ANTI_TRITON_INELASTIC, "inelastic"}, {MC_PROC_DEUTERON_INELASTIC, "inelastic"}, {MC_PROC_ELECTRON_NUCLEAR, "nuclear"},
        {MC_PROC_PHOTON_NUCLEAR, "nuclear"}, {MC_PROC_KAON_PLUS_INELASTIC, "inelastic"}, {MC_PROC_KAON_MINUS_INELASTIC, "inelastic"},
        {MC_PROC_HAD_BERTINI_CAPTURE_AT_REST, "capture"}, {MC_PROC_LAMBDA_INELASTIC, "inelastic"}, {MC_PROC_MU_NUCLEAR, "nuclear"},
        {MC_PROC_TRITON_INELASTIC, "inelastic"}, {MC_PROC_PRIMARY_BACKGROUND, "background"}};

    const std::map<std::string, int> categoryToColorMap = {{"invisible", 0}, {"primary", 1}, {"compt", 2}, {"phot", 3}, {"annihil", 4},
        {"ioni", 5}, {"brem", 6}, {"conv", 3}, {"capture", 6}, {"inelastic", 9}, {"elastic", 8}, {"decay", 7}, {"coulomb", 9},
        {"pair_prod", 4}, {"transport", 1}, {"rayleigh", 9}, {"kill", 2}, {"nuclear", 5}, {"background", 7}};

    LArHierarchyHelper::MCHierarchy::NodeVector nodes;
    hierarchy.GetFlattenedNodes(nodes);

    int nodeIdx{0};
    for (const LArHierarchyHelper::MCHierarchy::Node *pNode : nodes)
    {
        const LArMCParticle *pMC{dynamic_cast<const LArMCParticle *>(pNode->GetLeadingMCParticle())};
        if (!pMC)
            continue;
        const MCProcess process{pMC->GetProcess()};
        const std::string category{procToCategoryMap.at(process)};
        const int pdg{std::abs(pNode->GetParticleId())};
        const int tier{LArMCParticleHelper::GetHierarchyTier(pMC)};

        CaloHitList uHits, vHits, wHits;
        this->FillHitLists(pNode->GetCaloHits(), uHits, vHits, wHits);
        std::string suffix{std::to_string(nodeIdx) + " (" + std::to_string(tier) + ") " + std::to_string(pdg) + " " + category + " " +
                           std::to_string(process)};
        if (process == MC_PROC_DECAY)
        {
            const MCParticleList &parentList{pMC->GetParentList()};
            if (!parentList.empty())
            {
                const MCParticle *pParent{parentList.front()};
                const int parentPdg{std::abs(pParent->GetParticleId())};
                suffix += " from " + std::to_string(parentPdg);
            }
        }
        this->Visualize(uHits, "u_" + suffix, categoryToColorMap.at(category));
        this->Visualize(vHits, "v_" + suffix, categoryToColorMap.at(category));
        this->Visualize(wHits, "w_" + suffix, categoryToColorMap.at(category));
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
        if (m_collectionOnly && hits.front()->GetHitType() != TPC_VIEW_W)
            return;
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
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VisualizeProcess", m_visualizeProcess));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PerformMatching", m_match));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CollectionOnly", m_collectionOnly));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TransparencyThresholdE", m_transparencyThresholdE));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "EnergyScaleThresholdE", m_energyScaleThresholdE));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ScalingFactor", m_scalingFactor));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
