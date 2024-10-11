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
    m_writeEventTree{false},
    m_writeMCTree{false},
    m_foldToPrimaries{false},
    m_foldDynamic{false},
    m_foldToLeadingShowers{false},
    m_validateEvent{false},
    m_validateMC{false},
    m_minPurity{0.8f},
    m_minCompleteness{0.65f},
    m_minRecoHits{30},
    m_minRecoHitsPerView{10},
    m_minRecoGoodViews{2},
    m_removeRecoNeutrons{true},
    m_selectRecoHits{false}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

HierarchyValidationAlgorithm::~HierarchyValidationAlgorithm()
{
    if (m_writeEventTree)
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_eventTreeName.c_str(), m_eventFileName.c_str(), "UPDATE"));
    }
    if (m_writeMCTree)
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_MCTreeName.c_str(), m_MCFileName.c_str(), "UPDATE"));
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
    const LArHierarchyHelper::MCHierarchy::ReconstructabilityCriteria recoCriteria(
        m_minRecoHits, m_minRecoHitsPerView, m_minRecoGoodViews, m_removeRecoNeutrons);

    LArHierarchyHelper::MCHierarchy mcHierarchy(recoCriteria);
    LArHierarchyHelper::FillMCHierarchy(*pMCParticleList, *pCaloHitList, foldParameters, mcHierarchy);
    LArHierarchyHelper::RecoHierarchy recoHierarchy;
    LArHierarchyHelper::FillRecoHierarchy(*pPfoList, foldParameters, recoHierarchy);
    const LArHierarchyHelper::QualityCuts quality(m_minPurity, m_minCompleteness, m_selectRecoHits);
    LArHierarchyHelper::MatchInfo matchInfo(mcHierarchy, recoHierarchy, quality);
    LArHierarchyHelper::MatchHierarchies(matchInfo);
    matchInfo.Print(mcHierarchy);

#ifdef MONITORING
    if (m_validateEvent)
        this->EventValidation(matchInfo);
    if (m_validateMC)
        this->MCValidation(matchInfo);
#endif

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

#ifdef MONITORING
void HierarchyValidationAlgorithm::EventValidation(const LArHierarchyHelper::MatchInfo &matchInfo) const
{
    if (m_writeEventTree)
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
            // Check primaries is not empty before proceeding
            if (primaries.size() == 0)
                continue;
            primaries.sort(LArMCParticleHelper::SortByMomentum);
	    
	    int isCC{0}; isQE{0}, isResonant{0}, isDIS{0}, isCoherent{0}, isNuMu{0}, isNuE{0}, nPiZero{0}, nPiPlus{0}, nPiMinus{0}, 
		nPhotons{0}, nProtons{0};

	    try
            {
                const InteractionDescriptor descriptor{LArInteractionTypeHelper::GetInteractionDescriptor(primaries)};
                isCC = descriptor.IsCC();
                isQE = descriptor.IsQE();
                isResonant = descriptor.IsResonant();
                isDIS = descriptor.IsDIS();
                isCoherent = descriptor.IsCoherent();
                isNuMu = descriptor.IsMuonNeutrino();
                isNuE = descriptor.IsElectronNeutrino();
                nPiZero = static_cast<int>(descriptor.GetNumPiZero());
                nPiPlus = static_cast<int>(descriptor.GetNumPiPlus());
                nPiMinus = static_cast<int>(descriptor.GetNumPiMinus());
                nPhotons = static_cast<int>(descriptor.GetNumPhotons());
                nProtons = static_cast<int>(descriptor.GetNumProtons());
            }
            catch (...)
            {
            }

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
            const float trueVtxX{trueVertex.GetX()};
            const float trueVtxY{trueVertex.GetY()};
            const float trueVtxZ{trueVertex.GetZ()};
            float recoVtxX{std::numeric_limits<float>::max()};
            float recoVtxY{std::numeric_limits<float>::max()};
            float recoVtxZ{std::numeric_limits<float>::max()};
            float vtxDx{std::numeric_limits<float>::max()};
            float vtxDy{std::numeric_limits<float>::max()};
            float vtxDz{std::numeric_limits<float>::max()};
            float vtxDr{std::numeric_limits<float>::max()};
            const int isFiducial{LArVertexHelper::IsInFiducialVolume(this->GetPandora(), trueVertex, m_detector)};

            for (const ParticleFlowObject *pRootPfo : matchedRecoSliceRoots)
            {
                try
                {
                    const CartesianVector &recoVertex{LArPfoHelper::GetVertex(pRootPfo)->GetPosition()};
                    const float recoVertexX{recoVertex.GetX()};
                    const float recoVertexY{recoVertex.GetY()};
                    const float recoVertexZ{recoVertex.GetZ()};
                    const float dx{recoVertexX - trueVtxX};
                    const float dy{recoVertexY - trueVtxY};
                    const float dz{recoVertexZ - trueVtxZ};
                    const float dr{std::sqrt(dx * dx + dy * dy + dz * dz)};
                    if (dr < vtxDr)
                    {
                        recoVtxX = recoVertexX;
                        recoVtxY = recoVertexY;
                        recoVtxZ = recoVertexZ;
                        vtxDx = dx;
                        vtxDy = dy;
                        vtxDz = dz;
                        vtxDr = dr;
                    }
                }
                catch (StatusCodeException &)
                {
                    if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                        std::cout << "HierarchyValidationAlgorithm::EventValidation could not retrieve recoVertex" << std::endl;
                }
            }

            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "event", m_event));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "interaction", interaction));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "nRecoSlices", nRecoSlices));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "isCC", isCC));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "isQE", isQE));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "isResonant", isResonant));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "isDIS", isDIS));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "isCoherent", isCoherent));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "isNuMu", isNuMu));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "isNuE", isNuE));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "nPiZero", nPiZero));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "nPiPlus", nPiPlus));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "nPiMinus", nPiMinus));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "nPhotons", nPhotons));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "nProtons", nProtons));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "isFiducial", isFiducial));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "recoVtxX", recoVtxX));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "recoVtxY", recoVtxY));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "recoVtxZ", recoVtxZ));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "trueVtxX", recoVtxX));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "trueVtxY", recoVtxY));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "trueVtxZ", recoVtxZ));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "vtxDx", vtxDx));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "vtxDy", vtxDy));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "vtxDz", vtxDz));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "vtxDr", vtxDr));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "nGoodMatches", nGoodMatches));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "nPoorMatches", nPoorMatches));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "nUnmatched", nUnmatched));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "nNodes", nNodes));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "nGoodTier1Matches", nGoodTier1Matches));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "nTier1Nodes", nTier1Nodes));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "nGoodTrackMatches", nGoodTrackMatches));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "nGoodShowerMatches", nGoodShowerMatches));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "nTrackNodes", nTrackNodes));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "nShowerNodes", nShowerNodes));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "nGoodTier1TrackMatches", nGoodTier1TrackMatches));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "nTier1TrackNodes", nTier1TrackNodes));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "nGoodTier1ShowerMatches", nGoodTier1ShowerMatches));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "nTier1ShowerNodes", nTier1ShowerNodes));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "hasLeadingMuon", hasLeadingMuon));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "hasLeadingElectron", hasLeadingElectron));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_eventTreeName.c_str(), "isLeadingLeptonCorrect", isLeadingLeptonCorrect));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_eventTreeName.c_str()));
            ++interaction;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HierarchyValidationAlgorithm::MCValidation(const LArHierarchyHelper::MatchInfo &matchInfo) const
{
    if (m_writeMCTree)
    {
        int interaction{0};
        MCParticleList rootMCParticles;
	const CartesianVector null(0.f,0.f,0.f);
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
            for (const LArHierarchyHelper::MCMatches &matches : matchInfo.GetMatches(pRoot))
            {
                const LArHierarchyHelper::MCHierarchy::Node *pMCNode{matches.GetMC()};
                if (pMCNode->GetHierarchyTier() == 1)
                {
                    const MCParticle *const pLeadingMC{pMCNode->GetLeadingMCParticle()};
                    primaries.emplace_back(pLeadingMC);
                }
            }
            // Check primaries is not empty before proceeding
            if (primaries.size() == 0)
                continue;
            primaries.sort(LArMCParticleHelper::SortByMomentum);

            for (const LArHierarchyHelper::MCMatches &matches : matchInfo.GetMatches(pRoot))
            {
                const LArHierarchyHelper::MCHierarchy::Node *pMCNode{matches.GetMC()};
                const int isSignal{LArMCParticleHelper::IsNeutrino(LArMCParticleHelper::GetParentMCParticle(pMCNode->GetLeadingMCParticle())) ? 1 : 0};
		const int isTestBeam{pMCNode->IsTestBeamParticle() ? 1 : 0};
                const int isCosmicRay{!isTestBeam && pMCNode->IsCosmicRay() ? 1 : 0};
                const int isNeutrinoInt{!(isTestBeam || isCosmicRay) ? 1 : 0};
                const int mcId{pMCNode->GetId()};
                const int pdg{pMCNode->GetParticleId()};
                const int tier{pMCNode->GetHierarchyTier()};
                const int mcHits{static_cast<int>(pMCNode->GetCaloHits().size())};
                const int isLeadingLepton{pMCNode->IsLeadingLepton() ? 1 : 0};
		const CaloHitList mcCaloHits{pMCNode->GetCaloHits()};
		float mcTotalMipEnergy{0};
		for (const CaloHit *pCaloHit : mcCaloHits)
		{
		    mcTotalMipEnergy += pCaloHit->GetMipEquivalentEnergy();
		}

                const MCParticle *const pLeadingMC{pMCNode->GetLeadingMCParticle()};
                const MCParticleList &parentList{pLeadingMC->GetParentList()};
                const int isElectron{std::abs(pLeadingMC->GetParticleId()) == E_MINUS ? 1 : 0};
                const int hasMuonParent{parentList.size() == 1 && std::abs(parentList.front()->GetParticleId()) == MU_MINUS ? 1 : 0};
                const int isMichel{isElectron && hasMuonParent && LArMCParticleHelper::IsDecay(pLeadingMC) ? 1 : 0};
                const CartesianVector mcMomentum{pLeadingMC->GetMomentum()};
                const float mcMtm{mcMomentum.GetMagnitude()};
                const float mcMtmX{mcMomentum.GetX()};
                const float mcMtmY{mcMomentum.GetY()};
                const float mcMtmZ{mcMomentum.GetZ()};
                const float mcEnergy{pLeadingMC->GetEnergy()};

                const LArHierarchyHelper::RecoHierarchy::NodeVector &nodeVector{matches.GetRecoMatches()};
                const int nMatches{static_cast<int>(nodeVector.size())};
                IntVector recoSliceIdVector, recoIdVector, nRecoHitsVector, nSharedHitsVector;
		FloatVector recoTotalAdcVector;
                FloatVector purityVector, completenessVector;
                FloatVector purityAdcVector, completenessAdcVector;
                FloatVector purityVectorU, purityVectorV, purityVectorW, completenessVectorU, completenessVectorV, completenessVectorW;
                FloatVector purityAdcVectorU, purityAdcVectorV, purityAdcVectorW, completenessAdcVectorU, completenessAdcVectorV, completenessAdcVectorW;
		const CartesianVector &trueVertex{isSignal == 1 ? pLeadingMC->GetVertex() : null};
                const float trueVtxX{trueVertex.GetX()};
                const float trueVtxY{trueVertex.GetY()};
                const float trueVtxZ{trueVertex.GetZ()};
                float recoVtxX{std::numeric_limits<float>::max()};
                float recoVtxY{std::numeric_limits<float>::max()};
                float recoVtxZ{std::numeric_limits<float>::max()};
                float vtxDx{std::numeric_limits<float>::max()};
                float vtxDy{std::numeric_limits<float>::max()};
                float vtxDz{std::numeric_limits<float>::max()};
                float vtxDr{std::numeric_limits<float>::max()};

                int isCC{0}; isQE{0}, isResonant{0}, isDIS{0}, isCoherent{0}, isNuMu{0}, isNuE{0}, nPiZero{0}, nPiPlus{0}, nPiMinus{0},
                    nPhotons{0}, nProtons{0};
                try
                {
                    const InteractionDescriptor descriptor{LArInteractionTypeHelper::GetInteractionDescriptor(primaries)};
                    isCC = descriptor.IsCC();
                    isQE = descriptor.IsQE();
                    isResonant = descriptor.IsResonant();
                    isDIS = descriptor.IsDIS();
                    isCoherent = descriptor.IsCoherent();
                    isNuMu = descriptor.IsMuonNeutrino();
                    isNuE = descriptor.IsElectronNeutrino();
                    nPiZero = static_cast<int>(descriptor.GetNumPiZero());
                    nPiPlus = static_cast<int>(descriptor.GetNumPiPlus());
                    nPiMinus = static_cast<int>(descriptor.GetNumPiMinus());
                    nPhotons = static_cast<int>(descriptor.GetNumPhotons());
                    nProtons = static_cast<int>(descriptor.GetNumProtons());
                }
                catch (...)
                {
                }

                for (const LArHierarchyHelper::RecoHierarchy::Node *pRecoNode : nodeVector)
                {
                    const int sliceId{static_cast<int>(
                        std::distance(rootPfos.begin(), std::find(rootPfos.begin(), rootPfos.end(), recoNodeToRootMap[pRecoNode])))};
                    recoSliceIdVector.emplace_back(sliceId);
                    recoIdVector.emplace_back(pRecoNode->GetParticleId());
                    nRecoHitsVector.emplace_back(static_cast<int>(matches.GetSelectedRecoHits(pRecoNode).size()));
		    const CaloHitList recoCaloHits{pRecoNode->GetCaloHits()};
                    float recoTotalMipEnergy{0};
                    for (const CaloHit *pCaloHit : recoCaloHits)
                    {
                        recoTotalMipEnergy += pCaloHit->GetMipEquivalentEnergy();
                    }
                    recoTotalAdcVector.emplace_back(recoTotalMipEnergy);
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
                        try
                        {
			    const CartesianVector recoVertex{isSignal == 1 ? LArPfoHelper::GetVertex(pRecoNode->GetLeadingPfo())->GetPosition() : null};
                            recoVtxX = recoVertex.GetX();
                            recoVtxY = recoVertex.GetY();
                            recoVtxZ = recoVertex.GetZ();
                            vtxDx = recoVtxX - trueVtxX;
                            vtxDy = recoVtxY - trueVtxY;
                            vtxDz = recoVtxZ - trueVtxZ;
                            vtxDr = std::sqrt(vtxDx * vtxDx + vtxDy * vtxDy + vtxDz * vtxDz);
                        }
                        catch (StatusCodeException &)
                        {
                            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                                std::cout << "HierarchyValidationAlgorithm::MCValidation could not retrieve recoVertex" << std::endl;
                        }
                    }
                }

                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "event", m_event));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "interaction", interaction));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "mcId", mcId));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "mcPDG", pdg));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "mcTier", tier));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "mcNHits", mcHits));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "mcMtm", mcMtm));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "mcMtmX", mcMtmX));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "mcMtmY", mcMtmY));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "mcMtmZ", mcMtmZ));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "mcEnergy", mcEnergy));
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "recoTotalAdc", &recoTotalAdcVector));
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "mcTotalAdc", mcTotalMipEnergy));
		PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "isSignal", isSignal));       
	       	PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "isNuInteraction", isNeutrinoInt));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "isCosmicRay", isCosmicRay));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "isTestBeam", isTestBeam));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "isLeadingLepton", isLeadingLepton));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "isMichel", isMichel));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "nMatches", nMatches));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "recoSliceIdVector", &recoSliceIdVector));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "recoIdVector", &recoIdVector));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "nRecoHitsVector", &nRecoHitsVector));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "nSharedHitsVector", &nSharedHitsVector));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "purityVector", &purityVector));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "completenessVector", &completenessVector));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "purityAdcVector", &purityAdcVector));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "completenessAdcVector", &completenessAdcVector));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "purityVectorU", &purityVectorU));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "purityVectorV", &purityVectorV));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "purityVectorW", &purityVectorW));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "completenessVectorU", &completenessVectorU));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "completenessVectorV", &completenessVectorV));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "completenessVectorW", &completenessVectorW));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "purityAdcVectorU", &purityAdcVectorU));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "purityAdcVectorV", &purityAdcVectorV));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "purityAdcVectorW", &purityAdcVectorW));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "completenessAdcVectorU", &completenessAdcVectorU));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "completenessAdcVectorV", &completenessAdcVectorV));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "completenessAdcVectorW", &completenessAdcVectorW));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "recoVtxX", recoVtxX));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "recoVtxY", recoVtxY));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "recoVtxZ", recoVtxZ));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "trueVtxX", trueVtxX));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "trueVtxY", trueVtxY));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "trueVtxZ", trueVtxZ));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "vtxDx", vtxDx));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "vtxDy", vtxDy));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "vtxDz", vtxDz));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "vtxDr", vtxDr));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "isCC", isCC));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "isQE", isQE));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "isResonant", isResonant));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "isDIS", isDIS));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "isCoherent", isCoherent));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "isNuMu", isNuMu));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "isNuE", isNuE));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "nPiZero", nPiZero));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "nPiPlus", nPiPlus));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "nPiMinus", nPiMinus));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "nPhotons", nPhotons));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_MCTreeName.c_str(), "nProtons", nProtons));
                PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_MCTreeName.c_str()));
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

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteEventTree", m_writeEventTree));
    if (m_writeEventTree)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "EventFileName", m_eventFileName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "EventTreeName", m_eventTreeName));
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteMCTree", m_writeMCTree));
    if (m_writeMCTree)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCFileName", m_MCFileName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCTreeName", m_MCTreeName));
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FoldToPrimaries", m_foldToPrimaries));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FoldDynamic", m_foldDynamic));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FoldToLeadingShowers", m_foldToLeadingShowers));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinPurity", m_minPurity));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinCompleteness", m_minCompleteness));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinRecoHits", m_minRecoHits));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinRecoHitsPerView", m_minRecoHitsPerView));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinRecoGoodViews", m_minRecoGoodViews));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "RemoveRecoNeutrons", m_removeRecoNeutrons));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SelectRecoHits", m_selectRecoHits));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
