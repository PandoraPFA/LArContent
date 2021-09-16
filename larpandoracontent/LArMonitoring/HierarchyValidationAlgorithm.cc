/**
 *  @file   larpandoracontent/LArMonitoring/HierarchyValidationAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/HierarchyValidationAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"

using namespace pandora;

namespace lar_content
{

HierarchyValidationAlgorithm::HierarchyValidationAlgorithm() :
    m_event{0},
    m_writeTree{false},
    m_foldToPrimaries{false},
    m_foldDynamic{false},
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
        foldParameters.m_foldToTier = true;
    else if (m_foldDynamic)
        foldParameters.m_foldDynamic = true;
    else if (m_foldToLeadingShowers)
        foldParameters.m_foldToLeadingShowers = true;
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

    ++m_event;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HierarchyValidationAlgorithm::EventValidation(const LArHierarchyHelper::MatchInfo &matchInfo) const
{
    if (m_writeTree)
    {
        const LArHierarchyHelper::MCMatchesVector &goodMatches{matchInfo.GetGoodMatches()};
        const LArHierarchyHelper::MCMatchesVector &aboveThresholdMatches{matchInfo.GetAboveThresholdMatches()};
        const LArHierarchyHelper::MCMatchesVector &subThresholdMatches{matchInfo.GetSubThresholdMatches()};
        const LArHierarchyHelper::MCHierarchy::NodeVector &unmatched{matchInfo.GetUnmatchedMC()};
        MCParticleSet primaryMCSet;
        for (const LArHierarchyHelper::MCMatches &mcMatch : goodMatches)
        {
            const MCParticle *const pMC{mcMatch.GetMC()->GetLeadingMCParticle()};
            primaryMCSet.insert(LArMCParticleHelper::GetPrimaryMCParticle(pMC));
        }
        for (const LArHierarchyHelper::MCMatches &mcMatch : aboveThresholdMatches)
        {
            const MCParticle *const pMC{mcMatch.GetMC()->GetLeadingMCParticle()};
            primaryMCSet.insert(LArMCParticleHelper::GetPrimaryMCParticle(pMC));
        }
        for (const LArHierarchyHelper::MCMatches &mcMatch : subThresholdMatches)
        {
            const MCParticle *const pMC{mcMatch.GetMC()->GetLeadingMCParticle()};
            primaryMCSet.insert(LArMCParticleHelper::GetPrimaryMCParticle(pMC));
        }
        for (const LArHierarchyHelper::MCHierarchy::Node *pNode : unmatched)
        {
            const MCParticle *const pMC{pNode->GetLeadingMCParticle()};
            primaryMCSet.insert(LArMCParticleHelper::GetPrimaryMCParticle(pMC));
        }

        MCParticleList primaryMCList;
        for (const MCParticle *const pMC : primaryMCSet)
            primaryMCList.emplace_back(pMC);
        const int interactionType{static_cast<int>(LArInteractionTypeHelper::GetInteractionType(primaryMCList))};

        const int nGoodMatches{static_cast<int>(goodMatches.size())};
        const int nAboveThresholdMatches{static_cast<int>(aboveThresholdMatches.size())};
        const int nSubThresholdMatches{static_cast<int>(subThresholdMatches.size())};
        const int nUnmatched{static_cast<int>(unmatched.size())};
        const int nNodes{static_cast<int>(matchInfo.GetNMCNodes())};
        int hasLeadingMuon{0}, hasLeadingElectron{0}, isLeadingLeptonCorrect{0};

        std::set<const LArHierarchyHelper::MCHierarchy::Node *> trackNodeSet, showerNodeSet;
        int nGoodTrackMatches{0}, nGoodShowerMatches{0};
        for (const LArHierarchyHelper::MCMatches &mcMatch : goodMatches)
        {
            const LArHierarchyHelper::MCHierarchy::Node *pNode{mcMatch.GetMC()};
            const int pdg{std::abs(pNode->GetParticleId())};
            if (pNode->IsLeadingLepton())
            {
                if (pdg == MU_MINUS)
                    hasLeadingMuon = 1;
                else if (pdg == E_MINUS)
                    hasLeadingElectron = 1;
                isLeadingLeptonCorrect = 1;
            }

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
        for (const LArHierarchyHelper::MCMatches &mcMatch : aboveThresholdMatches)
        {
            const LArHierarchyHelper::MCHierarchy::Node *pNode{mcMatch.GetMC()};
            const int pdg{std::abs(pNode->GetParticleId())};
            if (pNode->IsLeadingLepton())
            {
                if (pdg == MU_MINUS)
                    hasLeadingMuon = 1;
                else if (pdg == E_MINUS)
                    hasLeadingElectron = 1;
            }

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

        for (const LArHierarchyHelper::MCMatches &mcMatch : subThresholdMatches)
        {
            const LArHierarchyHelper::MCHierarchy::Node *pNode{mcMatch.GetMC()};
            const int pdg{std::abs(pNode->GetParticleId())};
            if (pNode->IsLeadingLepton())
            {
                if (pdg == MU_MINUS)
                    hasLeadingMuon = 1;
                else if (pdg == E_MINUS)
                    hasLeadingElectron = 1;
            }

            if (pdg == PHOTON || pdg == E_MINUS)
                showerNodeSet.insert(pNode);
            else
                trackNodeSet.insert(pNode);
        }

        for (const LArHierarchyHelper::MCMatches &mcMatch : unmatched)
        {
            const LArHierarchyHelper::MCHierarchy::Node *pNode{mcMatch.GetMC()};
            const int pdg{std::abs(pNode->GetParticleId())};
            if (pNode->IsLeadingLepton())
            {
                if (pdg == MU_MINUS)
                    hasLeadingMuon = 1;
                else if (pdg == E_MINUS)
                    hasLeadingElectron = 1;
            }

            if (pdg == PHOTON || pdg == E_MINUS)
                showerNodeSet.insert(pNode);
            else
                trackNodeSet.insert(pNode);
        }

        const int nTrackNodes{static_cast<int>(trackNodeSet.size())}, nShowerNodes{static_cast<int>(showerNodeSet.size())};
        const CartesianVector &trueVertex{matchInfo.GetMCNeutrino()->GetVertex()};
        const CartesianVector &recoVertex{LArPfoHelper::GetVertex(matchInfo.GetRecoNeutrino())->GetPosition()};
        const float vtxDx{recoVertex.GetX() - trueVertex.GetX()};
        const float vtxDy{recoVertex.GetY() - trueVertex.GetY()};
        const float vtxDz{recoVertex.GetZ() - trueVertex.GetZ()};
        const float vtxDr{std::sqrt(vtxDx * vtxDx + vtxDy * vtxDy + vtxDz * vtxDz)};

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "event", m_event));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "interactionType", interactionType));
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
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "hasLeadingMuon", hasLeadingMuon));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "hasLeadingElectron", hasLeadingElectron));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isLeadingLeptonCorrect", isLeadingLeptonCorrect));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDx", vtxDx));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDy", vtxDy));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDz", vtxDz));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDr", vtxDr));
        PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HierarchyValidationAlgorithm::MCValidation(const LArHierarchyHelper::MatchInfo &matchInfo) const
{
    if (m_writeTree)
    {
        for (const LArHierarchyHelper::MCMatches &match : matchInfo.GetGoodMatches())
            this->FillMatched(match, true, true, matchInfo);
        for (const LArHierarchyHelper::MCMatches &match : matchInfo.GetAboveThresholdMatches())
            this->FillMatched(match, false, true, matchInfo);
        for (const LArHierarchyHelper::MCMatches &match : matchInfo.GetSubThresholdMatches())
            this->FillMatched(match, false, false, matchInfo);
        for (const LArHierarchyHelper::MCHierarchy::Node *pNode : matchInfo.GetUnmatchedMC())
            this->FillUnmatchedMC(pNode, matchInfo);
        for (const LArHierarchyHelper::RecoHierarchy::Node *pNode : matchInfo.GetUnmatchedReco())
            this->FillUnmatchedReco(pNode);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HierarchyValidationAlgorithm::FillMatched(const LArHierarchyHelper::MCMatches &matches, const bool isGood, const bool isAboveThreshold,
    const LArHierarchyHelper::MatchInfo &matchInfo) const
{
    const LArHierarchyHelper::MCHierarchy::Node *pMCNode{matches.GetMC()};
    const int isTestBeam{pMCNode->IsTestBeamParticle() ? 1 : 0};
    const int isCosmicRay{!isTestBeam && pMCNode->IsCosmicRay() ? 1 : 0};
    const int isNeutrinoInt{!(isTestBeam || isCosmicRay) ? 1 : 0};
    const int mcId{pMCNode->GetId()};
    const int pdg{pMCNode->GetParticleId()};
    const int tier{pMCNode->GetHierarchyTier()};
    const int mcHits{static_cast<int>(pMCNode->GetCaloHits().size())};
    const int isLeadingLepton{pMCNode->IsLeadingLepton() ? 1 : 0};

    const MCParticleList &parentList{pMCNode->GetLeadingMCParticle()->GetParentList()};
    const int isElectron{std::abs(pMCNode->GetLeadingMCParticle()->GetParticleId()) == E_MINUS ? 1 : 0};
    const int hasMuonParent{parentList.size() == 1 && std::abs(parentList.front()->GetParticleId()) == MU_MINUS ? 1 : 0};
    const int isMichel{isElectron && hasMuonParent && LArMCParticleHelper::IsDecay(pMCNode->GetLeadingMCParticle()) ? 1 : 0};

    LArHierarchyHelper::MCHierarchy::NodeVector goodMatchesInEvent;
    for (const auto goodMatches : matchInfo.GetGoodMatches())
    {
        if (goodMatches.GetMC() != pMCNode)
            goodMatchesInEvent.emplace_back(goodMatches.GetMC());
    }

    FloatVector expectedChildrenVector;
    IntVector goodRecoChildrenVector;
    for (const LArHierarchyHelper::MCHierarchy::Node *pChildNode : pMCNode->GetChildren())
    {
        expectedChildrenVector.emplace_back(pChildNode->GetId());
        if (std::find(goodMatchesInEvent.begin(), goodMatchesInEvent.end(), pChildNode) != goodMatchesInEvent.end())
            goodRecoChildrenVector.emplace_back(1);
        else
            goodRecoChildrenVector.emplace_back(0);
    }

    const LArHierarchyHelper::RecoHierarchy::NodeVector &nodeVector{matches.GetRecoMatches()};
    const int isGoodMatch{isGood};
    const int isAboveThresholdMatch{isAboveThreshold};
    const int nMatches{static_cast<int>(nodeVector.size())};
    IntVector recoIdVector, nRecoHitsVector, nSharedHitsVector;
    FloatVector purityVector, completenessVector;
    FloatVector purityAdcVector, completenessAdcVector;
    FloatVector purityVectorU, purityVectorV, purityVectorW, completenessVectorU, completenessVectorV, completenessVectorW;
    FloatVector purityAdcVectorU, purityAdcVectorV, purityAdcVectorW, completenessAdcVectorU, completenessAdcVectorV, completenessAdcVectorW;
    const CartesianVector &trueVertex{pMCNode->GetLeadingMCParticle()->GetVertex()};
    float vtxDx{0.f}, vtxDy{0.f}, vtxDz{0.f}, vtxDr{0.f};
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
        if (isGood)
        {
            // Only makes sense to calculate vertex delta if we have a one-to-one match
            const CartesianVector &recoVertex{LArPfoHelper::GetVertex(matchInfo.GetRecoNeutrino())->GetPosition()};
            vtxDx = recoVertex.GetX() - trueVertex.GetX();
            vtxDy = recoVertex.GetY() - trueVertex.GetY();
            vtxDz = recoVertex.GetZ() - trueVertex.GetZ();
            vtxDr = std::sqrt(vtxDx * vtxDx + vtxDy * vtxDy + vtxDz * vtxDz);
        }
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "event", m_event));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcId", mcId));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcPDG", pdg));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcTier", tier));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcNHits", mcHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isNuInteration", isNeutrinoInt));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isCosmicRay", isCosmicRay));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isTestBeam", isTestBeam));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isLeadingLepton", isLeadingLepton));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isMichel", isMichel));
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
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "expectedChildrenVector", &expectedChildrenVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "goodRecoChildrenVector", &goodRecoChildrenVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDx", vtxDx));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDy", vtxDy));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDz", vtxDz));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDr", vtxDr));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HierarchyValidationAlgorithm::FillUnmatchedMC(const LArHierarchyHelper::MCHierarchy::Node *pNode, const LArHierarchyHelper::MatchInfo &matchInfo) const
{
    const int isTestBeam{pNode->IsTestBeamParticle() ? 1 : 0};
    const int isCosmicRay{!isTestBeam && pNode->IsCosmicRay() ? 1 : 0};
    const int isNeutrinoInt{!(isTestBeam || isCosmicRay) ? 1 : 0};
    const int mcId{pNode->GetId()};
    const int pdg{pNode->GetParticleId()};
    const int tier{pNode->GetHierarchyTier()};
    const int mcHits{static_cast<int>(pNode->GetCaloHits().size())};
    const int isLeadingLepton{pNode->IsLeadingLepton() ? 1 : 0};

    const MCParticleList &parentList{pNode->GetLeadingMCParticle()->GetParentList()};
    const int isElectron{std::abs(pNode->GetLeadingMCParticle()->GetParticleId()) == E_MINUS ? 1 : 0};
    const int hasMuonParent{parentList.size() == 1 && std::abs(parentList.front()->GetParticleId()) == MU_MINUS ? 1 : 0};
    const int isMichel{isElectron && hasMuonParent && LArMCParticleHelper::IsDecay(pNode->GetLeadingMCParticle()) ? 1 : 0};

    LArHierarchyHelper::MCHierarchy::NodeVector goodMatchesInEvent;
    for (const auto goodMatches : matchInfo.GetGoodMatches())
    {
        if (goodMatches.GetMC() != pNode)
            goodMatchesInEvent.emplace_back(goodMatches.GetMC());
    }

    FloatVector expectedChildrenVector;
    IntVector goodRecoChildrenVector;
    for (const LArHierarchyHelper::MCHierarchy::Node *pChildNode : pNode->GetChildren())
    {
        expectedChildrenVector.emplace_back(pChildNode->GetId());
        if (std::find(goodMatchesInEvent.begin(), goodMatchesInEvent.end(), pChildNode) != goodMatchesInEvent.end())
            goodRecoChildrenVector.emplace_back(1);
        else
            goodRecoChildrenVector.emplace_back(0);
    }

    const int nMatches{0};
    const int isGoodMatch{0};
    const int isAboveThresholdMatch{0};
    IntVector recoIdVector, nRecoHitsVector, nSharedHitsVector;
    FloatVector purityVector, completenessVector;
    const float vtxDx{0.f}, vtxDy{0.f}, vtxDz{0.f}, vtxDr{0.f};

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "event", m_event));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcId", mcId));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcPDG", pdg));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcTier", tier));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcNHits", mcHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isNuInteration", isNeutrinoInt));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isCosmicRay", isCosmicRay));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isTestBeam", isTestBeam));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isLeadingLepton", isLeadingLepton));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isMichel", isMichel));
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
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "expectedChildrenVector", &expectedChildrenVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "goodRecoChildrenVector", &goodRecoChildrenVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDx", vtxDx));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDy", vtxDy));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDz", vtxDz));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDr", vtxDr));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HierarchyValidationAlgorithm::FillUnmatchedReco(const LArHierarchyHelper::RecoHierarchy::Node *pNode) const
{
    const int isTestBeam{0};
    const int isCosmicRay{0};
    const int isNeutrinoInt{0};
    const int mcId{0};
    const int pdg{0};
    const int tier{0};
    const int mcHits{0};
    const int isLeadingLepton{0};
    const int isMichel{0};
    FloatVector expectedChildrenVector;
    IntVector goodRecoChildrenVector;

    const int nMatches{0};
    const int isGoodMatch{0};
    const int isAboveThresholdMatch{0};
    IntVector recoIdVector{pNode->GetParticleId()}, nRecoHitsVector{static_cast<int>(pNode->GetCaloHits().size())}, nSharedHitsVector{0};
    FloatVector purityVector{0.f}, completenessVector{0.f};
    const float vtxDx{0.f}, vtxDy{0.f}, vtxDz{0.f}, vtxDr{0.f};

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "event", m_event));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcId", mcId));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcPDG", pdg));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcTier", tier));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcNHits", mcHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isNuInteration", isNeutrinoInt));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isCosmicRay", isCosmicRay));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isTestBeam", isTestBeam));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isLeadingLepton", isLeadingLepton));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isMichel", isMichel));
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
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "expectedChildrenVector", &expectedChildrenVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "goodRecoChildrenVector", &goodRecoChildrenVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDx", vtxDx));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDy", vtxDy));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDz", vtxDz));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDr", vtxDr));
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
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FoldDynamic", m_foldDynamic));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FoldToLeadingShowers", m_foldToLeadingShowers));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
