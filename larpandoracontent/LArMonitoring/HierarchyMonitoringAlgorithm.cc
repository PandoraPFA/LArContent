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
    m_foldToPrimaries{false},
    m_foldDynamic{true},
    m_minPurity{0.8f},
    m_minCompleteness{0.65f},
    m_minMatchCompleteness{0.1f},
    m_transparencyThresholdE{-1.f},
    m_energyScaleThresholdE{1.f},
    m_scalingFactor{1.f}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

HierarchyMonitoringAlgorithm::~HierarchyMonitoringAlgorithm()
{
    if (!m_rootFileName.empty())
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "processes", m_rootFileName, "UPDATE"));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HierarchyMonitoringAlgorithm::Run()
{
#ifdef MONITORING
    PANDORA_MONITORING_API(
        SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, m_transparencyThresholdE, m_energyScaleThresholdE, m_scalingFactor));

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));

    LArHierarchyHelper::MCHierarchy::ReconstructabilityCriteria criteria(15, 5, 2, false);
    LArHierarchyHelper::MCHierarchy mcHierarchy(criteria);
    LArHierarchyHelper::RecoHierarchy recoHierarchy;
    LArHierarchyHelper::FoldingParameters foldParameters;
    if (m_foldToPrimaries)
        foldParameters.m_foldToTier = true;
    else if (m_foldDynamic)
        foldParameters.m_foldDynamic = true;

    if (m_visualizeMC || m_match)
    {
        LArHierarchyHelper::FillMCHierarchy(*pMCParticleList, *pCaloHitList, foldParameters, mcHierarchy);
        std::cout << mcHierarchy.ToString() << std::endl;
    }
    if (m_visualizeReco || m_match)
    {
        const PfoList *pPfoList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_pfoListName, pPfoList));
        LArHierarchyHelper::FillRecoHierarchy(*pPfoList, foldParameters, recoHierarchy);
        std::cout << recoHierarchy.ToString() << std::endl;
    }
    if (m_match)
    {
        LArHierarchyHelper::QualityCuts quality(m_minPurity, m_minCompleteness);
        LArHierarchyHelper::MatchInfo matchInfo(mcHierarchy, recoHierarchy, quality);
        LArHierarchyHelper::MatchHierarchies(matchInfo);
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
#endif

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
#ifdef MONITORING
void HierarchyMonitoringAlgorithm::VisualizeMC(const LArHierarchyHelper::MCHierarchy &hierarchy) const
{
    const std::map<int, const std::string> keys = {{13, "mu"}, {11, "e"}, {22, "gamma"}, {321, "kaon"}, {211, "pi"}, {2212, "p"}};
    const std::map<std::string, int> colors = {{"mu", 5}, {"e", 2}, {"gamma", 9}, {"kaon", 1}, {"pi", 3}, {"p", 4}, {"other", 14}};

    MCParticleList rootMCParticles;
    hierarchy.GetRootMCParticles(rootMCParticles);
    for (const MCParticle *const pRoot : rootMCParticles)
    {
        LArHierarchyHelper::MCHierarchy::NodeVector nodes;
        hierarchy.GetFlattenedNodes(pRoot, nodes);

        bool depositedHits{false};
        int nodeIdx{0};
        for (const LArHierarchyHelper::MCHierarchy::Node *pNode : nodes)
        {
            std::string key("other");
            const int pdg{std::abs(pNode->GetParticleId())};
            if (keys.find(pdg) != keys.end())
                key = keys.at(pdg);

            CaloHitList uHits, vHits, wHits;
            this->FillHitLists(pNode->GetCaloHits(), uHits, vHits, wHits);
            if (!(uHits.empty() && vHits.empty() && wHits.empty()))
            {
                depositedHits = true;
                std::string suffix{std::to_string(nodeIdx) + "_" + key};
                this->Visualize(uHits, "u_" + suffix, colors.at(key));
                this->Visualize(vHits, "v_" + suffix, colors.at(key));
                this->Visualize(wHits, "w_" + suffix, colors.at(key));
            }
            ++nodeIdx;
        }
        if (depositedHits)
        {
            PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HierarchyMonitoringAlgorithm::VisualizeMCDistinct(const LArHierarchyHelper::MCHierarchy &hierarchy) const
{
    const std::map<int, const std::string> keys = {{13, "mu"}, {11, "e"}, {22, "gamma"}, {321, "kaon"}, {211, "pi"}, {2212, "p"}};
    const int nColours{9};
    const int colors[nColours] = {1, 2, 3, 4, 5, 6, 7, 8, 9};

    MCParticleList rootMCParticles;
    hierarchy.GetRootMCParticles(rootMCParticles);
    for (const MCParticle *const pRoot : rootMCParticles)
    {
        LArHierarchyHelper::MCHierarchy::NodeVector nodes;
        hierarchy.GetFlattenedNodes(pRoot, nodes);

        bool depositedHits{false};
        int nodeIdx{0}, colorIdx{0};
        for (const LArHierarchyHelper::MCHierarchy::Node *pNode : nodes)
        {
            std::string key("other");
            const int pdg{std::abs(pNode->GetParticleId())};
            if (keys.find(pdg) != keys.end())
                key = keys.at(pdg);

            CaloHitList uHits, vHits, wHits;
            this->FillHitLists(pNode->GetCaloHits(), uHits, vHits, wHits);
            if (!(uHits.empty() && vHits.empty() && wHits.empty()))
            {
                depositedHits = true;
                std::string suffix{std::to_string(nodeIdx) + "_" + key};
                this->Visualize(uHits, "u_" + suffix, colors[colorIdx]);
                this->Visualize(vHits, "v_" + suffix, colors[colorIdx]);
                this->Visualize(wHits, "w_" + suffix, colors[colorIdx]);
                colorIdx = (colorIdx + 1) >= nColours ? 0 : colorIdx + 1;
            }
            ++nodeIdx;
        }

        if (depositedHits)
        {
            PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        }
    }
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

    MCParticleList rootMCParticles;
    hierarchy.GetRootMCParticles(rootMCParticles);
    for (const MCParticle *const pRoot : rootMCParticles)
    {
        LArHierarchyHelper::MCHierarchy::NodeVector nodes;
        hierarchy.GetFlattenedNodes(pRoot, nodes);

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

            this->Visualize(uHits, "U " + suffix, categoryToColorMap.at(category));
            this->Visualize(vHits, "V " + suffix, categoryToColorMap.at(category));
            this->Visualize(wHits, "W " + suffix, categoryToColorMap.at(category));
            ++nodeIdx;

            const int proc{static_cast<int>(process)};
            const float mom{pMC->GetMomentum().GetMagnitude()};
            if (!m_rootFileName.empty())
            {
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "processes", "process", proc));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "processes", "momentum", mom));
                PANDORA_MONITORING_API(FillTree(this->GetPandora(), "processes"));
            }
        }

        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HierarchyMonitoringAlgorithm::VisualizeReco(const LArHierarchyHelper::RecoHierarchy &hierarchy) const
{
    const int nColors{7};
    const int colors[nColors] = {5, 2, 9, 1, 3, 4, 14};

    PfoList rootPfos;
    hierarchy.GetRootPfos(rootPfos);
    for (const ParticleFlowObject *const pRoot : rootPfos)
    {
        LArHierarchyHelper::RecoHierarchy::NodeVector nodes;
        hierarchy.GetFlattenedNodes(pRoot, nodes);

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
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HierarchyMonitoringAlgorithm::VisualizeMatches(const LArHierarchyHelper::MatchInfo &matchInfo) const
{
    const int nColors{8};
    const int colors[nColors] = {2, 3, 4, 5, 6, 7, 8, 9};
    MCParticleList rootMCParticles;
    matchInfo.GetRootMCParticles(rootMCParticles);
    PfoList rootPfos;
    const LArHierarchyHelper::MCHierarchy &mcHierarchy{matchInfo.GetMCHierarchy()};
    const LArHierarchyHelper::RecoHierarchy &recoHierarchy{matchInfo.GetRecoHierarchy()};
    recoHierarchy.GetRootPfos(rootPfos);
    std::map<const LArHierarchyHelper::RecoHierarchy::Node *, int> recoNodeToColorMap;
    std::map<const LArHierarchyHelper::RecoHierarchy::Node *, int> recoNodeToIdMap;
    std::map<const ParticleFlowObject *, int> recoRootToSliceMap;
    std::map<const LArHierarchyHelper::RecoHierarchy::Node *, const ParticleFlowObject *> recoNodeToRootMap;
    std::map<const LArHierarchyHelper::MCHierarchy::Node *, int> mcNodeToIdMap;
    int colorIdx{0}, recoIdx{0}, sliceIdx{0}, mcIdx{0};
    for (const ParticleFlowObject *const pRoot : rootPfos)
    {
        recoRootToSliceMap[pRoot] = sliceIdx;
        LArHierarchyHelper::RecoHierarchy::NodeVector nodes;
        recoHierarchy.GetFlattenedNodes(pRoot, nodes);
        for (const LArHierarchyHelper::RecoHierarchy::Node *pNode : nodes)
        {
            recoNodeToColorMap[pNode] = colors[colorIdx];
            recoNodeToIdMap[pNode] = recoIdx;
            recoNodeToRootMap[pNode] = pRoot;
            ++recoIdx;
        }
        ++sliceIdx;
        ++colorIdx;
        if (colorIdx >= nColors)
            colorIdx = 0;
    }
    for (const MCParticle *const pRoot : rootMCParticles)
    {
        bool isReconstructable{false};
        LArHierarchyHelper::MCHierarchy::NodeVector nodes;
        mcHierarchy.GetFlattenedNodes(pRoot, nodes);

        // Display MC
        for (const LArHierarchyHelper::MCHierarchy::Node *pNode : nodes)
        {
            const MCParticle *pMC{pNode->GetLeadingMCParticle()};
            if (!pMC)
                continue;
            const int pdg{std::abs(pNode->GetParticleId())};

            CaloHitList uHits, vHits, wHits;
            this->FillHitLists(pNode->GetCaloHits(), uHits, vHits, wHits);
            std::string suffix{"MC - ID: " + std::to_string(mcIdx) + " PDG: " + std::to_string(pdg)};
            this->Visualize(uHits, "U " + suffix, 1);
            this->Visualize(vHits, "V " + suffix, 1);
            this->Visualize(wHits, "W " + suffix, 1);
            if (!(uHits.empty() && vHits.empty() && wHits.empty()))
                isReconstructable = true;
            mcNodeToIdMap[pNode] = mcIdx;
            ++mcIdx;
        }

        // Display reco
        const LArHierarchyHelper::MCMatchesVector &matches{matchInfo.GetMatches(pRoot)};
        std::map<const LArHierarchyHelper::RecoHierarchy::Node *, bool> matchedReco;
        std::map<const ParticleFlowObject *, bool> matchedRoots;
        for (const auto &match : matches)
        {
            const LArHierarchyHelper::MCHierarchy::Node *pMC{match.GetMC()};
            for (const LArHierarchyHelper::RecoHierarchy::Node *pReco : match.GetRecoMatches())
            {
                CaloHitList uHits, vHits, wHits;
                const float purity{match.GetPurity(pReco, true)};
                const float completeness{match.GetCompleteness(pReco, true)};
                const std::string purityStr{this->ToStringSF(purity)};
                const std::string completenessStr{this->ToStringSF(completeness)};

                this->FillHitLists(pReco->GetCaloHits(), uHits, vHits, wHits);
                const std::string suffix{"Reco - Slice " + std::to_string(recoRootToSliceMap[recoNodeToRootMap[pReco]]) + " ID " +
                    std::to_string(recoNodeToIdMap[pReco]) + " -> " + std::to_string(mcNodeToIdMap[pMC])};
                const LArHierarchyHelper::QualityCuts &quality{matchInfo.GetQualityCuts()};
                if (purity >= quality.m_minPurity && completeness >= quality.m_minCompleteness)
                {
                    this->Visualize(uHits, "U " + suffix + "(" + purityStr + "," + completenessStr + ")", recoNodeToColorMap[pReco]);
                    this->Visualize(vHits, "V " + suffix + "(" + purityStr + "," + completenessStr + ")", recoNodeToColorMap[pReco]);
                    this->Visualize(wHits, "W " + suffix + "(" + purityStr + "," + completenessStr + ")", recoNodeToColorMap[pReco]);

                    matchedReco[pReco] = true;
                    matchedRoots[recoNodeToRootMap[pReco]] = true;
                }
                else
                {
                    // Only note the reconstructed slice if the completeness passes a minimal threshold
                    if (completeness >= m_minMatchCompleteness)
                    {
                        matchedReco[pReco] = true;
                        matchedRoots[recoNodeToRootMap[pReco]] = true;

                        this->Visualize(uHits, "U " + suffix + "(" + purityStr + "," + completenessStr + ")", 14);
                        this->Visualize(vHits, "V " + suffix + "(" + purityStr + "," + completenessStr + ")", 14);
                        this->Visualize(wHits, "W " + suffix + "(" + purityStr + "," + completenessStr + ")", 14);
                    }
                }
            }
        }

        // Display unmatched reco (for this true slice)
        for (const auto [pRecoRoot, val] : matchedRoots)
        {
            (void)val;
            LArHierarchyHelper::RecoHierarchy::NodeVector recoNodes;
            recoHierarchy.GetFlattenedNodes(pRecoRoot, recoNodes);
            for (const LArHierarchyHelper::RecoHierarchy::Node *pNode : recoNodes)
            {
                if (matchedReco.find(pNode) == matchedReco.end())
                {
                    CaloHitList uHits, vHits, wHits;
                    this->FillHitLists(pNode->GetCaloHits(), uHits, vHits, wHits);
                    const std::string suffix{"Reco - Slice " + std::to_string(recoRootToSliceMap[recoNodeToRootMap[pNode]]) + " ID " +
                        std::to_string(recoNodeToIdMap[pNode]) + " -> #"};
                    this->Visualize(uHits, "U " + suffix, 14);
                    this->Visualize(vHits, "V " + suffix, 14);
                    this->Visualize(wHits, "W " + suffix, 14);
                }
            }
        }

        if (isReconstructable)
        {
            PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HierarchyMonitoringAlgorithm::Visualize(const CaloHitList &hits, const std::string &label, const int color) const
{
    if (!hits.empty())
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
#endif

//------------------------------------------------------------------------------------------------------------------------------------------

std::string HierarchyMonitoringAlgorithm::ToStringSF(const float val, const int sf) const
{
    std::ostringstream out;
    out.precision(sf);
    out << std::fixed << val;
    return out.str();
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
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "RootFileName", m_rootFileName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VisualizeMC", m_visualizeMC));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VisualizeReco", m_visualizeReco));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VisualizeDistinct", m_visualizeDistinct));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VisualizeProcess", m_visualizeProcess));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PerformMatching", m_match));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CollectionOnly", m_collectionOnly));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FoldToPrimaries", m_foldToPrimaries));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FoldDynamic", m_foldDynamic));
    if (m_foldToPrimaries)
        m_foldDynamic = false;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinPurity", m_minPurity));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinCompleteness", m_minCompleteness));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinMatchCompleteness", m_minMatchCompleteness));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TransparencyThresholdE", m_transparencyThresholdE));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "EnergyScaleThresholdE", m_energyScaleThresholdE));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ScalingFactor", m_scalingFactor));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
