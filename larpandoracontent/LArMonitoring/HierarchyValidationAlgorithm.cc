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
#include "larpandoracontent/LArHelpers/LArVertexHelper.h"

using namespace pandora;

namespace lar_content
{

HierarchyValidationAlgorithm::HierarchyValidationAlgorithm() :
    m_event{-1},
    m_detector{"dune_fd_hd"},
    m_writeTree{false},
    m_foldToPrimaries{false},
    m_foldDynamic{false},
    m_foldToLeadingShowers{false},
    m_validateEvent{false},
    m_validateMC{false},
    m_minPurity{0.5},
    m_minCompleteness{0.5},
    m_minHits{10},
    m_minHitsForGoodView{5}
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
    ++m_event;
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
    StatusCode status{PandoraContentApi::GetList(*this, m_pfoListName, pPfoList)};
    LArHierarchyHelper::QualityCuts quality(m_minPurity, m_minCompleteness);
    LArHierarchyHelper::MCHierarchy::ReconstructabilityCriteria criteria;
    criteria.m_minHits = m_minHits;
    criteria.m_minHitsForGoodView = m_minHitsForGoodView;
    LArHierarchyHelper::MatchInfo matchInfo(mcHierarchy, recoHierarchy, quality);
    LArHierarchyHelper::MatchHierarchies(matchInfo);
    matchInfo.Print(mcHierarchy);

#ifdef MONITORING
    if (m_validateEvent)
        this->EventValidation(matchInfo);
    else if (m_validateMC)
        this->MCValidation(matchInfo);
#endif
    if (status != STATUS_CODE_SUCCESS)
        delete pPfoList;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

#ifdef MONITORING
void HierarchyValidationAlgorithm::EventValidation(const LArHierarchyHelper::MatchInfo &matchInfo) const
{
    if (m_writeTree)
    {
        int interaction{0};
        MCParticleList rootMCParticles;
        matchInfo.GetRootMCParticles(rootMCParticles);

        const LArHierarchyHelper::RecoHierarchy &recoHierarchy{matchInfo.GetRecoHierarchy()};
        PfoList rootPfos;
        recoHierarchy.GetRootPfos(rootPfos);
        std::map<const LArHierarchyHelper::RecoHierarchy::Node *, const ParticleFlowObject *> recoNodeToRootMap;
        for (const ParticleFlowObject *const pRoot : rootPfos)
        {
            LArHierarchyHelper::RecoHierarchy::NodeVector nodes;
            recoHierarchy.GetFlattenedNodes(pRoot, nodes);
            for (const LArHierarchyHelper::RecoHierarchy::Node *pNode : nodes)
                recoNodeToRootMap[pNode] = pRoot;
        }

        std::set<const pandora::ParticleFlowObject *> matchedRecoSliceRoots;
        for (const MCParticle *const pRoot : rootMCParticles)
        {
            const LArHierarchyHelper::MCMatchesVector &matches{matchInfo.GetMatches(pRoot)};

            MCParticleList primaries;
            for (const LArHierarchyHelper::MCMatches &match : matches)
            {
                const LArHierarchyHelper::MCHierarchy::Node *pMCNode{match.GetMC()};
                if (pMCNode->GetHierarchyTier() == 1)
                {
                    const MCParticle *const pLeadingMC{pMCNode->GetLeadingMCParticle()};
                    primaries.emplace_back(pLeadingMC);
                }
            }
            primaries.sort(LArMCParticleHelper::SortByMomentum);
            const InteractionDescriptor descriptor{LArInteractionTypeHelper::GetInteractionDescriptor(primaries)};

            const int isCC{descriptor.IsCC()};
            const int isQE{descriptor.IsQE()};
            const int isResonant{descriptor.IsResonant()};
            const int isDIS{descriptor.IsDIS()};
            const int isCoherent{descriptor.IsCoherent()};
            const int isNuMu{descriptor.IsMuonNeutrino()};
            const int isNuE{descriptor.IsElectronNeutrino()};
            const int nPiZero{static_cast<int>(descriptor.GetNumPiZero())};
            const int nPiPlus{static_cast<int>(descriptor.GetNumPiPlus())};
            const int nPiMinus{static_cast<int>(descriptor.GetNumPiMinus())};
            const int nPhotons{static_cast<int>(descriptor.GetNumPhotons())};
            const int nProtons{static_cast<int>(descriptor.GetNumProtons())};

            std::set<const LArHierarchyHelper::MCHierarchy::Node *> trackNodeSet, showerNodeSet;
            int nGoodTrackMatches{0}, nGoodShowerMatches{0};
            int nGoodMatches{0}, nPoorMatches{0}, nUnmatched{0};
            int nGoodTier1Matches{0}, nTier1Nodes{0};
            int nGoodTier1TrackMatches{0}, nTier1TrackNodes{0};
            int nGoodTier1ShowerMatches{0}, nTier1ShowerNodes{0};
            int hasLeadingMuon{0}, hasLeadingElectron{0}, isLeadingLeptonCorrect{0};
            for (const LArHierarchyHelper::MCMatches &mcMatch : matches)
            {
                const LArHierarchyHelper::MCHierarchy::Node *pNode{mcMatch.GetMC()};

                for (const LArHierarchyHelper::RecoHierarchy::Node *pReco : mcMatch.GetRecoMatches())
                    matchedRecoSliceRoots.insert(recoNodeToRootMap[pReco]);

                const int nReco{static_cast<int>(mcMatch.GetRecoMatches().size())};
                const bool isQuality{mcMatch.IsQuality(matchInfo.GetQualityCuts())};
                if (nReco == 1 && isQuality)
                    ++nGoodMatches;
                else if (nReco >= 1)
                    ++nPoorMatches;
                else
                    ++nUnmatched;
                if (pNode->GetHierarchyTier() == 1)
                {
                    ++nTier1Nodes;
                    if (nReco == 1 && isQuality)
                        ++nGoodTier1Matches;
                }

                const int pdg{std::abs(pNode->GetParticleId())};
                if (pNode->IsLeadingLepton())
                {
                    if (pdg == MU_MINUS)
                        hasLeadingMuon = 1;
                    else if (pdg == E_MINUS)
                        hasLeadingElectron = 1;
                    isLeadingLeptonCorrect = nReco == 1 ? 1 : 0;
                }

                if (pdg == PHOTON || pdg == E_MINUS)
                {
                    showerNodeSet.insert(pNode);
                    if (nReco == 1 && isQuality)
                    {
                        ++nGoodShowerMatches;
                        if (pNode->GetHierarchyTier() == 1)
                            ++nGoodTier1ShowerMatches;
                    }
                    if (pNode->GetHierarchyTier() == 1)
                        ++nTier1ShowerNodes;
                }
                else
                {
                    trackNodeSet.insert(pNode);
                    if (nReco == 1 && isQuality)
                    {
                        ++nGoodTrackMatches;
                        if (pNode->GetHierarchyTier() == 1)
                            ++nGoodTier1TrackMatches;
                    }
                    if (pNode->GetHierarchyTier() == 1)
                        ++nTier1TrackNodes;
                }
            }

            const int nNodes{static_cast<int>(matchInfo.GetNMCNodes(pRoot))};
            const int nTrackNodes{static_cast<int>(trackNodeSet.size())}, nShowerNodes{static_cast<int>(showerNodeSet.size())};
            const int nRecoSlices{static_cast<int>(matchedRecoSliceRoots.size())};
            const CartesianVector &trueVertex{pRoot->GetVertex()};
            float vtxDx{std::numeric_limits<float>::max()};
            float vtxDy{std::numeric_limits<float>::max()};
            float vtxDz{std::numeric_limits<float>::max()};
            float vtxDr{std::numeric_limits<float>::max()};
            const int isFiducial{LArVertexHelper::IsInFiducialVolume(this->GetPandora(), trueVertex, m_detector)};

            for (const ParticleFlowObject * pRootPfo : matchedRecoSliceRoots)
            {
                const CartesianVector &recoVertex{LArPfoHelper::GetVertex(pRootPfo)->GetPosition()};
                const float dx{recoVertex.GetX() - trueVertex.GetX()};
                const float dy{recoVertex.GetY() - trueVertex.GetY()};
                const float dz{recoVertex.GetZ() - trueVertex.GetZ()};
                const float dr{std::sqrt(dx * dx + dy * dy + dz * dz)};
                if (dr < vtxDr)
                {
                    vtxDx = dx;
                    vtxDy = dy;
                    vtxDz = dz;
                    vtxDr = dr;
                }
            }

            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "event", m_event));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "interaction", interaction));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nRecoSlices", nRecoSlices));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isCC", isCC));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isQE", isQE));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isResonant", isResonant));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isDIS", isDIS));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isCoherent", isCoherent));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isNuMu", isNuMu));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isNuE", isNuE));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nPiZero", nPiZero));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nPiPlus", nPiPlus));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nPiMinus", nPiMinus));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nPhotons", nPhotons));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nProtons", nProtons));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isFiducial", isFiducial));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDx", vtxDx));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDy", vtxDy));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDz", vtxDz));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDr", vtxDr));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nGoodMatches", nGoodMatches));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nPoorMatches", nPoorMatches));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nUnmatched", nUnmatched));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nNodes", nNodes));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nGoodTier1Matches", nGoodTier1Matches));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nTier1Nodes", nTier1Nodes));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nGoodTrackMatches", nGoodTrackMatches));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nGoodShowerMatches", nGoodShowerMatches));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nTrackNodes", nTrackNodes));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nShowerNodes", nShowerNodes));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nGoodTier1TrackMatches", nGoodTier1TrackMatches));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nTier1TrackNodes", nTier1TrackNodes));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nGoodTier1ShowerMatches", nGoodTier1ShowerMatches));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nTier1ShowerNodes", nTier1ShowerNodes));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "hasLeadingMuon", hasLeadingMuon));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "hasLeadingElectron", hasLeadingElectron));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isLeadingLeptonCorrect", isLeadingLeptonCorrect));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
            ++interaction;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HierarchyValidationAlgorithm::MCValidation(const LArHierarchyHelper::MatchInfo &matchInfo) const
{
    if (m_writeTree)
    {
        int interaction{0};
        MCParticleList rootMCParticles;
        matchInfo.GetRootMCParticles(rootMCParticles);
        PfoList rootPfos;
        const LArHierarchyHelper::RecoHierarchy &recoHierarchy{matchInfo.GetRecoHierarchy()};
        recoHierarchy.GetRootPfos(rootPfos);
        std::map<const LArHierarchyHelper::RecoHierarchy::Node *, const ParticleFlowObject *> recoNodeToRootMap;
        for (const ParticleFlowObject *const pRoot : rootPfos)
        {
            LArHierarchyHelper::RecoHierarchy::NodeVector nodes;
            recoHierarchy.GetFlattenedNodes(pRoot, nodes);
            for (const LArHierarchyHelper::RecoHierarchy::Node *pNode : nodes)
                recoNodeToRootMap[pNode] = pRoot;
        }

        for (const MCParticle *const pRoot : rootMCParticles)
        {
            MCParticleList primaries;
	    CaloHitList allHits;
            for (const LArHierarchyHelper::MCMatches &matches : matchInfo.GetMatches(pRoot))
            {
                const LArHierarchyHelper::MCHierarchy::Node *pMCNode{matches.GetMC()};
                if (pMCNode->GetHierarchyTier() == 1)
                {
                    const MCParticle *const pLeadingMC{pMCNode->GetLeadingMCParticle()};
                    primaries.emplace_back(pLeadingMC);
                }
		//allHits.emplace_back(pMCNode->GetCaloHits());
		//std::cout << "Number of hits" << allHits.size() << std::endl;
            }
            primaries.sort(LArMCParticleHelper::SortByMomentum);
            const InteractionDescriptor descriptor{LArInteractionTypeHelper::GetInteractionDescriptor(primaries)};

            for (const LArHierarchyHelper::MCMatches &matches : matchInfo.GetMatches(pRoot))
            {
                const LArHierarchyHelper::MCHierarchy::Node *pMCNode{matches.GetMC()};
                const int isTestBeam{pMCNode->IsTestBeamParticle() ? 1 : 0};
                const int isCosmicRay{!isTestBeam && pMCNode->IsCosmicRay() ? 1 : 0};
                const int isNeutrinoInt{!(isTestBeam || isCosmicRay) ? 1 : 0};
		//const int mcEnergy{pLeadingMC->GetEnergy()};
                const int mcId{pMCNode->GetId()};
                const int pdg{pMCNode->GetParticleId()};
                const int tier{pMCNode->GetHierarchyTier()};
                const int mcHits{static_cast<int>(pMCNode->GetCaloHits().size())};
                const int isLeadingLepton{pMCNode->IsLeadingLepton() ? 1 : 0};

                const MCParticle *const pLeadingMC{pMCNode->GetLeadingMCParticle()};
                const MCParticleList &parentList{pLeadingMC->GetParentList()};
                const int isElectron{std::abs(pLeadingMC->GetParticleId()) == E_MINUS ? 1 : 0};
                const int hasMuonParent{parentList.size() == 1 && std::abs(parentList.front()->GetParticleId()) == MU_MINUS ? 1 : 0};
                const int isMichel{isElectron && hasMuonParent && LArMCParticleHelper::IsDecay(pLeadingMC) ? 1 : 0};
                const float mcMomentum{pLeadingMC->GetMomentum().GetMagnitude()};

                const LArHierarchyHelper::RecoHierarchy::NodeVector &nodeVector{matches.GetRecoMatches()};
                const int nMatches{static_cast<int>(nodeVector.size())};
                IntVector recoSliceIdVector, recoIdVector, nRecoHitsVector, nSharedHitsVector;
                FloatVector purityVector, completenessVector;
                FloatVector purityAdcVector, completenessAdcVector;
                FloatVector purityVectorU, purityVectorV, purityVectorW, completenessVectorU, completenessVectorV, completenessVectorW;
                FloatVector purityAdcVectorU, purityAdcVectorV, purityAdcVectorW, completenessAdcVectorU, completenessAdcVectorV, completenessAdcVectorW;
                const CartesianVector &trueVertex{pLeadingMC->GetVertex()};
                float vtxDx{0.f}, vtxDy{0.f}, vtxDz{0.f}, vtxDr{0.f};

                const int isCC{descriptor.IsCC()};
                const int isQE{descriptor.IsQE()};
                const int isResonant{descriptor.IsResonant()};
                const int isDIS{descriptor.IsDIS()};
                const int isCoherent{descriptor.IsCoherent()};
                const int isNuMu{descriptor.IsMuonNeutrino()};
                const int isNuE{descriptor.IsElectronNeutrino()};
                const int nPiZero{static_cast<int>(descriptor.GetNumPiZero())};
                const int nPiPlus{static_cast<int>(descriptor.GetNumPiPlus())};
                const int nPiMinus{static_cast<int>(descriptor.GetNumPiMinus())};
                const int nPhotons{static_cast<int>(descriptor.GetNumPhotons())};
                const int nProtons{static_cast<int>(descriptor.GetNumProtons())};

                for (const LArHierarchyHelper::RecoHierarchy::Node *pRecoNode : nodeVector)
                {
                    const int sliceId{static_cast<int>(std::distance(rootPfos.begin(), std::find(rootPfos.begin(), rootPfos.end(),
                        recoNodeToRootMap[pRecoNode])))};
                    recoSliceIdVector.emplace_back(sliceId);
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
                    if (nMatches > 0)
                    {
                        const CartesianVector recoVertex{LArPfoHelper::GetVertex(pRecoNode->GetLeadingPfo())->GetPosition()};
                        vtxDx = recoVertex.GetX() - trueVertex.GetX();
                        vtxDy = recoVertex.GetY() - trueVertex.GetY();
                        vtxDz = recoVertex.GetZ() - trueVertex.GetZ();
                        vtxDr = std::sqrt(vtxDx * vtxDx + vtxDy * vtxDy + vtxDz * vtxDz);
                    }
                }
  
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "event", m_event));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "interaction", interaction));
                //PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcEnergy", mcEnergy));
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcId", mcId));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcPDG", pdg));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcTier", tier));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcNHits", mcHits));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcMomentum", mcMomentum));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isNuInteration", isNeutrinoInt));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isCosmicRay", isCosmicRay));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isTestBeam", isTestBeam));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isLeadingLepton", isLeadingLepton));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isMichel", isMichel));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nMatches", nMatches));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "recoSliceIdVector", &recoSliceIdVector));
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
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDx", vtxDx));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDy", vtxDy));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDz", vtxDz));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDr", vtxDr));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isCC", isCC));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isQE", isQE));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isResonant", isResonant));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isDIS", isDIS));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isCoherent", isCoherent));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isNuMu", isNuMu));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isNuE", isNuE));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nPiZero", nPiZero));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nPiPlus", nPiPlus));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nPiMinus", nPiMinus));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nPhotons", nPhotons));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nProtons", nProtons));
                PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
            }
            ++interaction;
        }
    }
}
#endif

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HierarchyValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    if (m_caloHitListName.empty())
        m_caloHitListName = "CaloHitList2D";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    if (m_pfoListName.empty())
        m_pfoListName = "RecreatedPfos";

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Detector", m_detector));

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

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinPurity", m_minPurity));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinCompleteness", m_minCompleteness));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinHits", m_minHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinHitsForGoodView", m_minHitsForGoodView));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
