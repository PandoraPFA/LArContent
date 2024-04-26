/**
 *  @file   LArReco/src/UGProjectAlgorithm.cxx
 *
 *  @brief  Implementation of the UGProjectAlgorithm class.
 *
 *  $Log: $
 */
#include "Pandora/AlgorithmHeaders.h"
#include "Helpers/MCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArObjectHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "UGProjectAlgorithm.h"

using namespace pandora;

namespace lar_content
{

UGProjectAlgorithm::UGProjectAlgorithm():
    m_writeToTree(),
    m_eventId(-1),
    m_treeName(),
    m_fileName(),
    m_mcParticleListName(),
    m_caloHitListName(),
    m_inputPfoListName()
{
}

UGProjectAlgorithm::~UGProjectAlgorithm()
{
    if (m_writeToTree)
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode UGProjectAlgorithm::Run()
{
    // ATTN - m_eventId is initialised to -1, so the first event is event zero
    ++m_eventId;

    // Input lists
    const PfoList *pPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputPfoListName, pPfoList));
    
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));
    
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    // True vertex extraction
    MCParticleVector primaries;
    LArMCParticleHelper::GetPrimaryMCParticleList(pMCParticleList, primaries);
    const MCParticle *pTrueNeutrino{nullptr};
    if (primaries.empty())
    {
        return STATUS_CODE_SUCCESS;
    }
    const MCParticle *primary{primaries.front()};
    const MCParticleList &parents{primary->GetParentList()};
    if (parents.size() == 1 && LArMCParticleHelper::IsNeutrino(parents.front()))
        pTrueNeutrino = parents.front();
    else
        return STATUS_CODE_SUCCESS;

    float trueNeutrinoVertexX{0.f}, trueNeutrinoVertexY{0.f}, trueNeutrinoVertexZ{0.f}, trueNeutrinoVertexU{0.f}, trueNeutrinoVertexV{0.f},
        trueNeutrinoVertexW{0.f}, trueNeutrinoEnergy{0.f};
    if (pTrueNeutrino)
    {
        const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()}; 
        const CartesianVector &trueVertex{pTrueNeutrino->GetVertex()};
        trueNeutrinoVertexX = trueVertex.GetX();
        trueNeutrinoVertexY = trueVertex.GetY();
        trueNeutrinoVertexZ = trueVertex.GetZ();
        trueNeutrinoVertexU = static_cast<float>(transform->YZtoU(trueVertex.GetY(), trueVertex.GetZ()));
        trueNeutrinoVertexV = static_cast<float>(transform->YZtoV(trueVertex.GetY(), trueVertex.GetZ()));
        trueNeutrinoVertexW = static_cast<float>(transform->YZtoW(trueVertex.GetY(), trueVertex.GetZ()));
        trueNeutrinoEnergy = pTrueNeutrino->GetEnergy();
    }

    // Get MC hits
    LArMCParticleHelper::MCContributionMap mcToHitsMap, mcToHitsMapU, mcToHitsMapV, mcToHitsMapW;
    for (const CaloHit *pCaloHit : *pCaloHitList)
    {
        try
        {
            const MCParticle *pMC{MCParticleHelper::GetMainMCParticle(pCaloHit)};
            mcToHitsMap[pMC].emplace_back(pCaloHit);
            switch (pCaloHit->GetHitType())
            {
                case TPC_VIEW_U:
                    mcToHitsMapU[pMC].emplace_back(pCaloHit);
                    break;
                case TPC_VIEW_V:
                    mcToHitsMapV[pMC].emplace_back(pCaloHit);
                    break;
                case TPC_VIEW_W:
                    mcToHitsMapW[pMC].emplace_back(pCaloHit);
                    break;
                default:
                    break;
            }
        }
        catch (const StatusCodeException &)
        {
        }
    }

    // Get event descriptor
    MCParticleList primariesList(primaries.begin(), primaries.end());
    const InteractionDescriptor descriptor{LArInteractionTypeHelper::GetInteractionDescriptor(primariesList)};

    int nTrueTracksU{0}, nTrueTracksV{0}, nTrueTracksW{0}, nTrueShowersU{0}, nTrueShowersV{0}, nTrueShowersW{0};
    for (const auto &[pMC, caloHits] : mcToHitsMap)
    {
        const int pdg{std::abs(pMC->GetParticleId())};
        bool isMichel{false};
        if (pdg == E_MINUS)
        {
            const MCParticleList &mcParents{pMC->GetParentList()};
            if (mcParents.size() == 1)
            {
                if (std::abs(mcParents.front()->GetParticleId() == MU_MINUS))
                    isMichel = true;
            }
        }
        const bool isTrack{!(pdg == E_MINUS || pdg == PHOTON) || (pdg == E_MINUS && isMichel)};
        if (isTrack)
            ++nTrueTracksU;
        else
            ++nTrueShowersU;
        if (isTrack)
            ++nTrueTracksV;
        else
            ++nTrueShowersV;
        if (isTrack)
            ++nTrueTracksW;
        else
            ++nTrueShowersW;
    }

    // Mapping reconstructed particles -> reconstruction associated Hits
    PfoList allConnectedPfos;
    LArPfoHelper::GetAllConnectedPfos(*pPfoList, allConnectedPfos);

    std::map<const ParticleFlowObject *, int> pfoToIdMap;
    int i{1};
    for (const ParticleFlowObject *const pPfo : allConnectedPfos)
    {
        if (!LArPfoHelper::IsNeutrino(pPfo))
        {
            pfoToIdMap[pPfo] = i;
            ++i;
        }
        else
            pfoToIdMap[pPfo] = 0;
    }
    
    PfoList finalStatePfos;
    for (const ParticleFlowObject *const pPfo : allConnectedPfos)
    {
        if (LArPfoHelper::IsFinalState(pPfo))
            finalStatePfos.push_back(pPfo);
    }
    int nFinalStateParticles{static_cast<int>(finalStatePfos.size())};
    LArMCParticleHelper::PfoContributionMap pfoToHitsMap;
    LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(allConnectedPfos, mcToHitsMap, pfoToHitsMap, false);
    
    // Matching step
    LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCHitSharingMap;
    LArMCParticleHelper::MCParticleToPfoHitSharingMap mcToPfoHitSharingMap;

    LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(pfoToHitsMap, {mcToHitsMap}, pfoToMCHitSharingMap, mcToPfoHitSharingMap);
    
    // Write one tree entry for each Pfo
    for (const Pfo *const pPfo : allConnectedPfos)
    {
        if (LArPfoHelper::IsNeutrino(pPfo))
            continue;
        const PfoList &parentList{pPfo->GetParentPfoList()};
        const ParticleFlowObject *pParentPfo{parentList.size() == 1 ? parentList.front() : nullptr};
        const int parentId{pParentPfo ? pfoToIdMap[pParentPfo] : -1};
        IntVector childPfos;
        const PfoList &childList{pPfo->GetDaughterPfoList()}; 
        for (const ParticleFlowObject *pChild : childList)
            childPfos.emplace_back(pfoToIdMap[pChild]);

        CaloHitList allHitsInPfo;
        LArPfoHelper::GetAllCaloHits(pPfo, allHitsInPfo);
        const int nHitsInPfoTotal(allHitsInPfo.size()),
            nHitsInPfoU(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, allHitsInPfo)),
            nHitsInPfoV(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, allHitsInPfo)),
            nHitsInPfoW(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, allHitsInPfo));
        CaloHitList uHitsInPfo, vHitsInPfo, wHitsInPfo, hitsInPfo3D;
        LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_U, uHitsInPfo);
        LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_V, vHitsInPfo);
        LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, wHitsInPfo);
        LArPfoHelper::GetCaloHits(pPfo, TPC_3D, hitsInPfo3D);
        
        int nHitsInBestMCParticleTotal(-1),
            nHitsInBestMCParticleU(-1),
            nHitsInBestMCParticleV(-1),
            nHitsInBestMCParticleW(-1),
            bestMCParticlePdgCode(0),
            bestMCParticleTier{-1};
        double bestMCParticleEnergy(-999999),
               bestMCParticleMomentumX(-999999),
               bestMCParticleMomentumY(-999999),
               bestMCParticleMomentumZ(-999999);
        int nHitsSharedWithBestMCParticleTotal(-1);
        const MCParticle *pBestMC{nullptr};
        
        const LArMCParticleHelper::MCParticleToSharedHitsVector &mcParticleToSharedHitsVector(pfoToMCHitSharingMap.at(pPfo));
        
        for (const LArMCParticleHelper::MCParticleCaloHitListPair &mcParticleCaloHitListPair : mcParticleToSharedHitsVector)
        {
            const pandora::MCParticle *const pAssociatedMCParticle(mcParticleCaloHitListPair.first);
            const CaloHitList &allMCHits(mcToHitsMap.at(pAssociatedMCParticle));
            const CaloHitList &associatedMCHits(mcParticleCaloHitListPair.second);
            
            CaloHitList associatedMCHitsW;
            for (const CaloHit *const pCaloHit : associatedMCHits)
            {
                if (TPC_VIEW_W == pCaloHit->GetHitType())
                    associatedMCHitsW.push_back(pCaloHit);
            }
            
            if (static_cast<int>(associatedMCHits.size()) > nHitsSharedWithBestMCParticleTotal)
            {
                // This is the current best matched MCParticle, to be stored
                nHitsSharedWithBestMCParticleTotal = associatedMCHits.size();
                nHitsInBestMCParticleTotal = allMCHits.size();
                nHitsInBestMCParticleU = LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, allMCHits);
                nHitsInBestMCParticleV = LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, allMCHits);
                nHitsInBestMCParticleW = LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, allMCHits);
                bestMCParticlePdgCode = pAssociatedMCParticle->GetParticleId();
                bestMCParticleEnergy = pAssociatedMCParticle->GetEnergy();
                bestMCParticleMomentumX = pAssociatedMCParticle->GetMomentum().GetX();
                bestMCParticleMomentumY = pAssociatedMCParticle->GetMomentum().GetY();
                bestMCParticleMomentumZ = pAssociatedMCParticle->GetMomentum().GetZ();
                bestMCParticleTier = LArMCParticleHelper::GetHierarchyTier(pAssociatedMCParticle);
                pBestMC = pAssociatedMCParticle;
            }
        }
        if (bestMCParticlePdgCode == 0)
        {
            std::cout << "Found 0 match PFO with " << pfoToMCHitSharingMap.at(pPfo).size() << std::endl;
        }
        bool isTrack(true);
        if (std::abs(bestMCParticlePdgCode) == E_MINUS || std::abs(bestMCParticlePdgCode) == PHOTON)
        {
            bool isMichel{false};
            if (std::abs(bestMCParticlePdgCode) == E_MINUS)
            {
                const MCParticleList &mcParents{pBestMC->GetParentList()};
                if (mcParents.size() == 1)
                {
                    if (std::abs(mcParents.front()->GetParticleId() == MU_MINUS))
                        isMichel = true;
                }
            }
            isTrack = isMichel;
        }
        // Cache values here
        const int isPandoraTrack(pPfo->GetParticleId() == MU_MINUS ? 1 : 0);
        const int pandoraTier{LArPfoHelper::GetHierarchyTier(pPfo)};
        const int isTrueTrack(isTrack ? 1 : 0);
        const float purity((nHitsInPfoTotal > 0) ? static_cast<float>(nHitsSharedWithBestMCParticleTotal) / static_cast<float>(nHitsInPfoTotal) : 0.f);
        const float completeness((nHitsInBestMCParticleTotal > 0) ? static_cast<float>(nHitsSharedWithBestMCParticleTotal)/static_cast<float>(nHitsInBestMCParticleTotal) : 0.f);
        FloatVector hitDriftPositionsU, hitWirePositionsU, hitEnergiesU,
                    hitDriftPositionsV, hitWirePositionsV, hitEnergiesV,
                    hitDriftPositionsW, hitWirePositionsW, hitEnergiesW,
                    hit3DPositionsX, hit3DPositionsY, hit3DPositionsZ;

        for (const CaloHit *const pCaloHit : uHitsInPfo)
        {
            hitDriftPositionsU.push_back(pCaloHit->GetPositionVector().GetX());
            hitWirePositionsU.push_back(pCaloHit->GetPositionVector().GetZ());
            hitEnergiesU.push_back(pCaloHit->GetInputEnergy());
        }
        for (const CaloHit *const pCaloHit : vHitsInPfo)
        {
            hitDriftPositionsV.push_back(pCaloHit->GetPositionVector().GetX());
            hitWirePositionsV.push_back(pCaloHit->GetPositionVector().GetZ());
            hitEnergiesV.push_back(pCaloHit->GetInputEnergy());
        }
        for (const CaloHit *const pCaloHit : wHitsInPfo)
        {
            hitDriftPositionsW.push_back(pCaloHit->GetPositionVector().GetX());
            hitWirePositionsW.push_back(pCaloHit->GetPositionVector().GetZ());
            hitEnergiesW.push_back(pCaloHit->GetInputEnergy());
        }
        for (const CaloHit *const pCaloHit : hitsInPfo3D)
        {
            hit3DPositionsX.push_back(pCaloHit->GetPositionVector().GetX());
            hit3DPositionsY.push_back(pCaloHit->GetPositionVector().GetY());
            hit3DPositionsZ.push_back(pCaloHit->GetPositionVector().GetZ());
        }
        //GetVertex
        const pandora::Vertex* vertex{nullptr};
        try
        {
            vertex = LArPfoHelper::GetVertex(pPfo);
        }
        catch (StatusCodeException &)
        {
        }
        CartesianVector vertexPosition(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
        if (vertex)
            vertexPosition = vertex->GetPosition();
        float vertexPositionX = vertexPosition.GetX();
        float vertexPositionY = vertexPosition.GetY();
        float vertexPositionZ = vertexPosition.GetZ();

        //Project vertex on U,V views
        CartesianVector vertexPositionU = vertex ? LArGeometryHelper::ProjectPosition(this->GetPandora(),vertexPosition,TPC_VIEW_U) : vertexPosition;
        CartesianVector vertexPositionV = vertex ? LArGeometryHelper::ProjectPosition(this->GetPandora(),vertexPosition,TPC_VIEW_V) : vertexPosition;
        CartesianVector vertexPositionW = vertex ? LArGeometryHelper::ProjectPosition(this->GetPandora(),vertexPosition,TPC_VIEW_W) : vertexPosition;

        float vertexDriftPosition=vertexPositionU.GetX();
        float vertexWirePositionU=vertexPositionU.GetZ();
        float vertexWirePositionV=vertexPositionV.GetZ();
        float vertexWirePositionW=vertexPositionW.GetZ();

        //Get interaction vertex
        const pandora::Vertex* intVertex = LArPfoHelper::GetVertex(LArPfoHelper::GetParentNeutrino(pPfo));
        CartesianVector intVertexPosition = intVertex->GetPosition();
        CartesianVector intVertexPositionU = LArGeometryHelper::ProjectPosition(this->GetPandora(), intVertexPosition, TPC_VIEW_U);
        CartesianVector intVertexPositionV = LArGeometryHelper::ProjectPosition(this->GetPandora(), intVertexPosition, TPC_VIEW_V);
        CartesianVector intVertexPositionW = LArGeometryHelper::ProjectPosition(this->GetPandora(), intVertexPosition, TPC_VIEW_W);

        const int isCC{descriptor.IsCC()};
        const int isNuE{descriptor.IsElectronNeutrino()};
        const int isNuMu{descriptor.IsMuonNeutrino()};
        // Write to tree here
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "eventId", m_eventId));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isCC", isCC));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isNuE", isNuE));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isNuMu", isNuMu));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nFinalStateParticles", nFinalStateParticles));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "pfoId", pfoToIdMap[pPfo]));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "pfoTier", pandoraTier));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "parent", parentId));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "children", &childPfos));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsInPfoTotal", nHitsInPfoTotal));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsInPfoU", nHitsInPfoU));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsInPfoV", nHitsInPfoV));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsInPfoW", nHitsInPfoW));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsInBestMCParticleTotal", nHitsInBestMCParticleTotal));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsInBestMCParticleU", nHitsInBestMCParticleU));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsInBestMCParticleV", nHitsInBestMCParticleV));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsInBestMCParticleW", nHitsInBestMCParticleW));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMCParticlePdgCode", bestMCParticlePdgCode));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMCParticleTier", bestMCParticleTier));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMCParticleEnergy", bestMCParticleEnergy));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMCParticleMomentumX", bestMCParticleMomentumX));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMCParticleMomentumY", bestMCParticleMomentumY));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMCParticleMomentumZ", bestMCParticleMomentumZ));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "purity", purity));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "completeness", completeness));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isPandoraTrack", isPandoraTrack));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isTrueTrack", isTrueTrack));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "xU", &hitDriftPositionsU));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "u", &hitWirePositionsU));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "adcU", &hitEnergiesU));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "xV", &hitDriftPositionsV));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "v", &hitWirePositionsV));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "adcV", &hitEnergiesV));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "xW", &hitDriftPositionsW));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "w", &hitWirePositionsW));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "adcW", &hitEnergiesW));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "x3d", &hit3DPositionsX));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "y3d", &hit3DPositionsY));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "z3d", &hit3DPositionsZ));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "vtxX3d", vertexPositionX));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "vtxY3d", vertexPositionY));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "vtxZ3d", vertexPositionZ));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "vtxX", vertexDriftPosition));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "vtxU", vertexWirePositionU));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "vtxV", vertexWirePositionV));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "vtxW", vertexWirePositionW));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nuVtxX3d", intVertexPosition.GetX()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nuVtxY3d", intVertexPosition.GetY()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nuVtxZ3d", intVertexPosition.GetZ()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nuVtxX", intVertexPosition.GetX()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nuVtxU", intVertexPositionU.GetZ()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nuVtxV", intVertexPositionV.GetZ()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nuVtxW", intVertexPositionW.GetZ()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trueNuVtxX3d", trueNeutrinoVertexX));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trueNuVtxY3d", trueNeutrinoVertexY));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trueNuVtxZ3d", trueNeutrinoVertexZ));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trueNuVtxX", trueNeutrinoVertexX));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trueNuVtxU", trueNeutrinoVertexU));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trueNuVtxV", trueNeutrinoVertexV));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trueNuVtxW", trueNeutrinoVertexW));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trueNuEnergy", trueNeutrinoEnergy));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTrueTracksU", nTrueTracksU));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTrueTracksV", nTrueTracksV));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTrueTracksW", nTrueTracksW));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTrueShowersU", nTrueShowersU));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTrueShowersV", nTrueShowersV));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nTrueShowersW", nTrueShowersW));
        PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
    }
    std::cout << std::endl;
    
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode UGProjectAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "WriteToTree", m_writeToTree));
    
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "InputPfoListName", m_inputPfoListName));
    
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "CaloHitListName", m_caloHitListName));
    
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "MCParticleListName", m_mcParticleListName));
    
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "TreeName", m_treeName));
    
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "FileName", m_fileName));
    
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
