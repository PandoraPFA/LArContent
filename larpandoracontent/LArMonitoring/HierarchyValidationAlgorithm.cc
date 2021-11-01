/**
 *  @file   larpandoracontent/LArMonitoring/HierarchyValidationAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/HierarchyValidationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

HierarchyValidationAlgorithm::HierarchyValidationAlgorithm() :
    m_writeTree{false},
    m_foldToPrimaries{false},
    m_foldToLeadingShowers{false},
    m_validateEvent{false},
    m_validateMC{false}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

HierarchyValidationAlgorithm::~HierarchyValidationAlgorithm()
{
    if (m_writeTree)
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treename.c_str(), m_filename.c_str(), "UPDATE"));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HierarchyValidationAlgorithm::Run()
{
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));
    const PfoList *pPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_pfoListName, pPfoList));

    LArHierarchyHelper::FoldingParameters foldParameters;
    if (m_foldToPrimaries)
    {
        foldParameters.m_foldToTier = true;
        foldParameters.m_tier = 1;
    }
    else if (m_foldToLeadingShowers)
    {
        foldParameters.m_foldToLeadingShowers = true;
    }
    LArHierarchyHelper::MCHierarchy mcHierarchy;
    LArHierarchyHelper::FillMCHierarchy(*pMCParticleList, *pCaloHitList, foldParameters, mcHierarchy);
    LArHierarchyHelper::RecoHierarchy recoHierarchy;
    LArHierarchyHelper::FillRecoHierarchy(*pPfoList, foldParameters, recoHierarchy);
    LArHierarchyHelper::MatchInfo matchInfo;
    LArHierarchyHelper::MatchHierarchies(mcHierarchy, recoHierarchy, matchInfo);
    matchInfo.Print(mcHierarchy);

    if (m_validateEvent)
        this->EventValidation(matchInfo);
    else if (m_validateMC)
        this->MCValidation(matchInfo);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HierarchyValidationAlgorithm::EventValidation(const LArHierarchyHelper::MatchInfo &matchInfo) const
{
    if (m_writeTree)
    {
        const int nGoodMatches{static_cast<int>(matchInfo.GetGoodMatches().size())};
        const int nAboveThresholdMatches{static_cast<int>(matchInfo.GetAboveThresholdMatches().size())};
        const int nSubThresholdMatches{static_cast<int>(matchInfo.GetSubThresholdMatches().size())};
        const int nUnmatched{static_cast<int>(matchInfo.GetUnmatchedMC().size())};
        const int nNodes{static_cast<int>(matchInfo.GetNMCNodes())};
        int hasLeadingLepton{0}, isLeadingLeptonCorrect{0};

        std::set<const LArHierarchyHelper::MCHierarchy::Node *> trackNodeSet, showerNodeSet;
        int nGoodTrackMatches{0}, nGoodShowerMatches{0};
        for (const LArHierarchyHelper::MCMatches &mcMatch : matchInfo.GetGoodMatches())
        {
            const LArHierarchyHelper::MCHierarchy::Node *pNode{mcMatch.GetMC()};
            if (pNode->IsLeadingLepton())
            {
                hasLeadingLepton = 1;
                isLeadingLeptonCorrect = 1;
            }
            const int pdg{std::abs(pNode->GetParticleId())};
            if (pdg == PHOTON || pdg == E_MINUS)
            {
                showerNodeSet.insert(pNode);
                ++nGoodShowerMatches;
            }
            else
            {
                trackNodeSet.insert(pNode);
                ++nGoodTrackMatches;
            }
        }

        int nAboveThresholdTrackMatches{0}, nAboveThresholdShowerMatches{0};
        for (const LArHierarchyHelper::MCMatches &mcMatch : matchInfo.GetAboveThresholdMatches())
        {
            const LArHierarchyHelper::MCHierarchy::Node *pNode{mcMatch.GetMC()};
            if (pNode->IsLeadingLepton())
                hasLeadingLepton = 1;
            const int pdg{std::abs(pNode->GetParticleId())};
            if (pdg == PHOTON || pdg == E_MINUS)
            {
                showerNodeSet.insert(pNode);
                ++nAboveThresholdShowerMatches;
            }
            else
            {
                trackNodeSet.insert(pNode);
                ++nAboveThresholdTrackMatches;
            }
        }

        for (const LArHierarchyHelper::MCMatches &mcMatch : matchInfo.GetSubThresholdMatches())
        {
            const LArHierarchyHelper::MCHierarchy::Node *pNode{mcMatch.GetMC()};
            if (pNode->IsLeadingLepton())
                hasLeadingLepton = 1;
            const int pdg{std::abs(pNode->GetParticleId())};
            if (pdg == PHOTON || pdg == E_MINUS)
                showerNodeSet.insert(pNode);
            else
                trackNodeSet.insert(pNode);
        }

        for (const LArHierarchyHelper::MCMatches mcMatch : matchInfo.GetUnmatchedMC())
        {
            const LArHierarchyHelper::MCHierarchy::Node *pNode{mcMatch.GetMC()};
            if (pNode->IsLeadingLepton())
                hasLeadingLepton = 1;
            const int pdg{std::abs(pNode->GetParticleId())};
            if (pdg == PHOTON || pdg == E_MINUS)
                showerNodeSet.insert(pNode);
            else
                trackNodeSet.insert(pNode);
        }

        const int nTrackNodes{static_cast<int>(trackNodeSet.size())}, nShowerNodes{static_cast<int>(showerNodeSet.size())};

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nGoodMatches", nGoodMatches));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nAboveThresholdMatches", nAboveThresholdMatches));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nSubThresholdMatches", nSubThresholdMatches));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nUnmatched", nUnmatched));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nNodes", nNodes));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nGoodTrackMatches", nGoodTrackMatches));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nGoodShowerMatches", nGoodShowerMatches));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nAboveThresholdTrackMatches", nAboveThresholdTrackMatches));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nAboveThresholdShowerMatches", nAboveThresholdShowerMatches));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nTrackNodes", nTrackNodes));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nShowerNodes", nShowerNodes));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "hasLeadingLepton", hasLeadingLepton));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isLeadingLeptonCorrect", isLeadingLeptonCorrect));
        PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HierarchyValidationAlgorithm::MCValidation(const LArHierarchyHelper::MatchInfo &matchInfo) const
{
    if (m_writeTree)
    {
        for (const LArHierarchyHelper::MCMatches &match : matchInfo.GetGoodMatches())
            this->FillMatched(match, true, true);
        for (const LArHierarchyHelper::MCMatches &match : matchInfo.GetAboveThresholdMatches())
            this->FillMatched(match, false, true);
        for (const LArHierarchyHelper::MCMatches &match : matchInfo.GetSubThresholdMatches())
            this->FillMatched(match, false, false);
        for (const LArHierarchyHelper::MCHierarchy::Node *pNode : matchInfo.GetUnmatchedMC())
            this->FillUnmatchedMC(pNode);
        for (const LArHierarchyHelper::RecoHierarchy::Node *pNode : matchInfo.GetUnmatchedReco())
            this->FillUnmatchedReco(pNode);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HierarchyValidationAlgorithm::FillMatched(const LArHierarchyHelper::MCMatches &matches, const bool isGood, const bool isAboveThreshold) const
{
    const LArHierarchyHelper::MCHierarchy::Node *pMCNode{matches.GetMC()};
    const int isTestBeam{pMCNode->IsTestBeamParticle() ? 1 : 0};
    const int isCosmicRay{!isTestBeam && pMCNode->IsCosmicRay() ? 1 : 0};
    const int isNeutrinoInt{!(isTestBeam || isCosmicRay) ? 1 : 0};
    const int pdg{pMCNode->GetParticleId()};
    const int mcHits{static_cast<int>(pMCNode->GetCaloHits().size())};
    const int isLeadingLepton{pMCNode->IsLeadingLepton() ? 1 : 0};

    const LArHierarchyHelper::RecoHierarchy::NodeVector &nodeVector{matches.GetRecoMatches()};
    const int isGoodMatch{isGood};
    const int isAboveThresholdMatch{isAboveThreshold};
    const int nMatches{static_cast<int>(nodeVector.size())};
    IntVector recoIdVector, nRecoHitsVector, nSharedHitsVector;
    FloatVector purityVector, completenessVector;
    FloatVector purityAdcVector, completenessAdcVector;
    FloatVector purityVectorU, purityVectorV, purityVectorW, completenessVectorU, completenessVectorV, completenessVectorW;
    FloatVector purityAdcVectorU, purityAdcVectorV, purityAdcVectorW, completenessAdcVectorU, completenessAdcVectorV, completenessAdcVectorW;
    for (const LArHierarchyHelper::RecoHierarchy::Node *pRecoNode : nodeVector)
    {
        recoIdVector.emplace_back(pRecoNode->GetParticleId());
        nRecoHitsVector.emplace_back(static_cast<int>(pRecoNode->GetCaloHits().size()));
        nSharedHitsVector.emplace_back(static_cast<int>(matches.GetSharedHits(pRecoNode)));
        purityVector.emplace_back(matches.GetPurity(pRecoNode));
        completenessVector.emplace_back(matches.GetCompleteness(pRecoNode));
        purityAdcVector.emplace_back(matches.GetPurity(pRecoNode, true));
        completenessAdcVector.emplace_back(matches.GetCompleteness(pRecoNode, true));
        purityVectorU.emplace_back(matches.GetPurity(pRecoNode, TPC_VIEW_U));
        purityVectorV.emplace_back(matches.GetPurity(pRecoNode, TPC_VIEW_V));
        purityVectorW.emplace_back(matches.GetPurity(pRecoNode, TPC_VIEW_W));
        completenessVectorU.emplace_back(matches.GetCompleteness(pRecoNode, TPC_VIEW_U));
        completenessVectorV.emplace_back(matches.GetCompleteness(pRecoNode, TPC_VIEW_V));
        completenessVectorW.emplace_back(matches.GetCompleteness(pRecoNode, TPC_VIEW_W));
        purityAdcVectorU.emplace_back(matches.GetPurity(pRecoNode, TPC_VIEW_U, true));
        purityAdcVectorV.emplace_back(matches.GetPurity(pRecoNode, TPC_VIEW_V, true));
        purityAdcVectorW.emplace_back(matches.GetPurity(pRecoNode, TPC_VIEW_W, true));
        completenessAdcVectorU.emplace_back(matches.GetCompleteness(pRecoNode, TPC_VIEW_U, true));
        completenessAdcVectorV.emplace_back(matches.GetCompleteness(pRecoNode, TPC_VIEW_V, true));
        completenessAdcVectorW.emplace_back(matches.GetCompleteness(pRecoNode, TPC_VIEW_W, true));
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcPDG", pdg));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcNHits", mcHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isNuInteration", isNeutrinoInt));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isCosmicRay", isCosmicRay));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isTestBeam", isTestBeam));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isLeadingLepton", isLeadingLepton));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isGoodMatch", isGoodMatch));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isAboveThresholdMatch", isAboveThresholdMatch));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nMatches", nMatches));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "recoIdVector", &recoIdVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nRecoHitsVector", &nRecoHitsVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nSharedHitsVector", &nSharedHitsVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityVector", &purityVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessVector", &completenessVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityAdcVector", &purityAdcVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessAdcVector", &completenessAdcVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityVectorU", &purityVectorU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityVectorV", &purityVectorV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityVectorW", &purityVectorW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessVectorU", &completenessVectorU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessVectorV", &completenessVectorV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessVectorW", &completenessVectorW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityAdcVectorU", &purityAdcVectorU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityAdcVectorV", &purityAdcVectorV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityAdcVectorW", &purityAdcVectorW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessAdcVectorU", &completenessAdcVectorU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessAdcVectorV", &completenessAdcVectorV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessAdcVectorW", &completenessAdcVectorW));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HierarchyValidationAlgorithm::FillUnmatchedMC(const LArHierarchyHelper::MCHierarchy::Node *pNode) const
{
    const int isTestBeam{pNode->IsTestBeamParticle() ? 1 : 0};
    const int isCosmicRay{!isTestBeam && pNode->IsCosmicRay() ? 1 : 0};
    const int isNeutrinoInt{!(isTestBeam || isCosmicRay) ? 1 : 0};
    const int pdg{pNode->GetParticleId()};
    const int mcHits{static_cast<int>(pNode->GetCaloHits().size())};
    const int isLeadingLepton{pNode->IsLeadingLepton() ? 1 : 0};

    const int nMatches{0};
    const int isGoodMatch{0};
    const int isAboveThresholdMatch{0};
    IntVector recoIdVector, nRecoHitsVector, nSharedHitsVector;
    FloatVector purityVector, completenessVector;

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcPDG", pdg));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcNHits", mcHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isNuInteration", isNeutrinoInt));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isCosmicRay", isCosmicRay));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isTestBeam", isTestBeam));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isLeadingLepton", isLeadingLepton));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isGoodMatch", isGoodMatch));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isAboveThresholdMatch", isAboveThresholdMatch));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nMatches", nMatches));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "recoIdVector", &recoIdVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nRecoHitsVector", &nRecoHitsVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nSharedHitsVector", &nSharedHitsVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityVector", &purityVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessVector", &completenessVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityAdcVector", &purityVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessAdcVector", &completenessVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityVectorU", &purityVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityVectorV", &purityVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityVectorW", &purityVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessVectorU", &completenessVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessVectorV", &completenessVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessVectorW", &completenessVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityAdcVectorU", &purityVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityAdcVectorV", &purityVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityAdcVectorW", &purityVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessAdcVectorU", &completenessVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessAdcVectorV", &completenessVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessAdcVectorW", &completenessVector));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HierarchyValidationAlgorithm::FillUnmatchedReco(const LArHierarchyHelper::RecoHierarchy::Node *pNode) const
{
    const int isTestBeam{0};
    const int isCosmicRay{0};
    const int isNeutrinoInt{0};
    const int pdg{0};
    const int mcHits{0};
    const int isLeadingLepton{0};

    const int nMatches{0};
    const int isGoodMatch{0};
    const int isAboveThresholdMatch{0};
    IntVector recoIdVector{pNode->GetParticleId()}, nRecoHitsVector{static_cast<int>(pNode->GetCaloHits().size())}, nSharedHitsVector{0};
    FloatVector purityVector{0.f}, completenessVector{0.f};

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcPDG", pdg));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcNHits", mcHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isNuInteration", isNeutrinoInt));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isCosmicRay", isCosmicRay));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isTestBeam", isTestBeam));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isLeadingLepton", isLeadingLepton));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isGoodMatch", isGoodMatch));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isAboveThresholdMatch", isAboveThresholdMatch));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nMatches", nMatches));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "recoIdVector", &recoIdVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nRecoHitsVector", &nRecoHitsVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nSharedHitsVector", &nSharedHitsVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityVector", &purityVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessVector", &completenessVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityAdcVector", &purityVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessAdcVector", &completenessVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityVectorU", &purityVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityVectorV", &purityVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityVectorW", &purityVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessVectorU", &completenessVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessVectorV", &completenessVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessVectorW", &completenessVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityAdcVectorU", &purityVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityAdcVectorV", &purityVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityAdcVectorW", &purityVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessAdcVectorU", &completenessVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessAdcVectorV", &completenessVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessAdcVectorW", &completenessVector));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HierarchyValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    if (m_caloHitListName.empty())
        m_caloHitListName = "CaloHitList2D";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    if (m_pfoListName.empty())
        m_pfoListName = "RecreatedPfos";

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ValidateEvent", m_validateEvent));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ValidateMC", m_validateMC));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteTree", m_writeTree));
    if (m_writeTree)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "FileName", m_filename));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TreeName", m_treename));
        if (!(m_validateEvent || m_validateMC))
        {
            std::cout << "Error: WriteTree requested but no tree names found" << std::endl;
            return STATUS_CODE_NOT_FOUND;
        }
        else if (m_validateEvent && m_validateMC)
        {
            std::cout << "Error: Both event-level and MC-level validation requested simulataneously" << std::endl;
            return STATUS_CODE_INVALID_PARAMETER;
        }
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FoldToPrimaries", m_foldToPrimaries));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FoldToLeadingShowers", m_foldToLeadingShowers));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
