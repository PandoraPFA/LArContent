/**
 *  @file   larpandoracontent/LArThreeReco/LArEventBuilding/PfoCharacterisationAlgorithm.cc
 * 
 *  @brief  Implementation of the pfo characterisation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingShowerFitResult.h"

#include "larpandoracontent/LArThreeDReco/LArEventBuilding/PfoCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

PfoCharacterisationAlgorithm::PfoCharacterisationAlgorithm() :
    m_writeToTree(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

PfoCharacterisationAlgorithm::~PfoCharacterisationAlgorithm()
{
    if (m_writeToTree)
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoCharacterisationAlgorithm::Run()
{
    PfoList tracksToShowers, showersToTracks;

    for (const std::string &pfoListName : m_inputPfoListNames)
    {
        const PfoList *pPfoList(nullptr);
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, pfoListName, pPfoList));

        if (!pPfoList || pPfoList->empty())
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "PfoCharacterisationAlgorithm: unable to find pfo list " << pfoListName << std::endl;

            continue;
        }

        for (const ParticleFlowObject *const pPfo : *pPfoList)
        {
            if (this->IsClearTrack(pPfo))
            {
                PandoraContentApi::ParticleFlowObject::Metadata metadata;
                metadata.m_particleId = MU_MINUS;
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, metadata));

                if (m_showerPfoListName == pfoListName)
                    showersToTracks.push_back(pPfo);
            }
            else
            {
                PandoraContentApi::ParticleFlowObject::Metadata metadata;
                metadata.m_particleId = E_MINUS;
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, metadata));

                if (m_trackPfoListName == pfoListName)
                    tracksToShowers.push_back(pPfo);
            }
        }
    }

    if (!tracksToShowers.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_trackPfoListName, m_showerPfoListName, tracksToShowers));

    if (!showersToTracks.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_showerPfoListName, m_trackPfoListName, showersToTracks));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool PfoCharacterisationAlgorithm::IsClearTrack(const ParticleFlowObject *const pPfo) const
{
    CaloHitList caloHitList;
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_U, caloHitList);
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_V, caloHitList);
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, caloHitList);

    MCParticleWeightMap mcParticleWeightMap;

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        for (const MCParticleWeightMap::value_type &mapEntry : pCaloHit->GetMCParticleWeightMap())
            mcParticleWeightMap[mapEntry.first] += mapEntry.second;
    }

    float bestWeight(0.f);
    float bestMCParticleLength(-1.f);
    const MCParticle *pBestMCParticle(nullptr);

    MCParticleList mcParticleList;
    for (const auto &mapEntry : mcParticleWeightMap) mcParticleList.push_back(mapEntry.first);
    mcParticleList.sort(LArMCParticleHelper::SortByMomentum);

    for (const MCParticle *const pMCParticle : mcParticleList)
    {
        const float weight(mcParticleWeightMap.at(pMCParticle));

        if (weight > bestWeight)
        {
            pBestMCParticle = pMCParticle;
            bestWeight = weight;
        }
    }

    bool isTrueTrack(false);

    if (pBestMCParticle)
    {
        const MCParticle *const pPrimaryMCParticle(LArMCParticleHelper::GetPrimaryMCParticle(pBestMCParticle));
        const int absParticleId(std::abs(pPrimaryMCParticle->GetParticleId()));
        isTrueTrack = ((MU_MINUS == absParticleId) || (PROTON == absParticleId) || (PI_PLUS == absParticleId));
        bestMCParticleLength = (pBestMCParticle->GetEndpoint() - pBestMCParticle->GetVertex()).GetMagnitude();
    }

    // Tree variables here
    //--------------------------------------------------------------------------------------------------------------------------------------
    const int trueTrack(isTrueTrack ? 1 : 0);
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trueTrack", trueTrack));

    const int inputParticleId(pPfo->GetParticleId());
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "inputParticleId", inputParticleId));

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMCParticleLength", bestMCParticleLength));

    int nHitsU(0), nHitsV(0), nHitsW(0);
    int nClustersU(0), nClustersV(0), nClustersW(0);
    float vertexDistanceU(-1.f), vertexDistanceV(-1.f), vertexDistanceW(-1.f);

    FloatVector xPositionsU, xPositionsV, xPositionsW, zPositionsU, zPositionsV, zPositionsW;
    float mipEnergyU(0.f), mipEnergyV(0.f), mipEnergyW(0.f);

    float straightLineLengthU(-1.f), straightLineLengthV(-1.f), straightLineLengthW(-1.f);
    float integratedPathLengthU(-1.f), integratedPathLengthV(-1.f), integratedPathLengthW(-1.f);

    float straightLineLength10U(-1.f), straightLineLength10V(-1.f), straightLineLength10W(-1.f);
    float integratedPathLength10U(-1.f), integratedPathLength10V(-1.f), integratedPathLength10W(-1.f);

    float showerFitWidthU(-1.f), showerFitWidthV(-1.f), showerFitWidthW(-1.f);
    float showerFitGapLengthU(-1.f), showerFitGapLengthV(-1.f), showerFitGapLengthW(-1.f);

    int nPointsOfContactU(0), nPointsOfContactV(0), nPointsOfContactW(0);
    int nHitsInBranchesU(0), nHitsInBranchesV(0), nHitsInBranchesW(0);

    ClusterList twoDClusterList;
    LArPfoHelper::GetTwoDClusterList(pPfo, twoDClusterList);

    for (const Cluster *const pCluster : twoDClusterList)
    {
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));

        if (!((TPC_VIEW_U == hitType) || (TPC_VIEW_V == hitType) || (TPC_VIEW_W == hitType)))
            continue;

        const std::string hitTypeLabel((TPC_VIEW_U == hitType) ? "U" : (TPC_VIEW_V == hitType) ? "V" : "W");

        // nHits
        int &nHits((TPC_VIEW_U == hitType) ? nHitsU : (TPC_VIEW_V == hitType) ? nHitsV : nHitsW);
        nHits += pCluster->GetNCaloHits();
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHits" + hitTypeLabel, nHits));

        // nClusters
        int &nClusters((TPC_VIEW_U == hitType) ? nClustersU : (TPC_VIEW_V == hitType) ? nClustersV : nClustersW);
        ++nClusters;
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nClusters" + hitTypeLabel, nClusters));

        // vertexDistance
        float &vertexDistance((TPC_VIEW_U == hitType) ? vertexDistanceU : (TPC_VIEW_V == hitType) ? vertexDistanceV : vertexDistanceW);

        const VertexList *pVertexList = nullptr;
        PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetCurrentList(*this, pVertexList));
        const Vertex *const pSelectedVertex((pVertexList && (pVertexList->size() == 1) && (VERTEX_3D == (*(pVertexList->begin()))->GetVertexType())) ? *(pVertexList->begin()) : nullptr);

        if (pSelectedVertex)
        {
            const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pSelectedVertex->GetPosition(), LArClusterHelper::GetClusterHitType(pCluster)));
            vertexDistance = LArClusterHelper::GetClosestDistance(vertexPosition2D, pCluster);
        }

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "vertexDistance" + hitTypeLabel, vertexDistance));

        // hit positions and energy
        FloatVector &xPositions((TPC_VIEW_U == hitType) ? xPositionsU : (TPC_VIEW_V == hitType) ? xPositionsV : xPositionsW);
        FloatVector &zPositions((TPC_VIEW_U == hitType) ? zPositionsU : (TPC_VIEW_V == hitType) ? zPositionsV : zPositionsW);
        float &mipEnergy((TPC_VIEW_U == hitType) ? mipEnergyU : (TPC_VIEW_V == hitType) ? mipEnergyV : mipEnergyW);

        CaloHitList clusterHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(clusterHitList);

        for (const CaloHit *const pCaloHit : clusterHitList)
        {
            xPositions.push_back(pCaloHit->GetPositionVector().GetX());
            zPositions.push_back(pCaloHit->GetPositionVector().GetZ());
            mipEnergy += pCaloHit->GetMipEquivalentEnergy();
        }

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "xPositions" + hitTypeLabel, &xPositions));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "zPositions" + hitTypeLabel, &zPositions));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mipEnergy" + hitTypeLabel, mipEnergy));

        // straight line length and integrated pathlength
        float &straightLineLength((TPC_VIEW_U == hitType) ? straightLineLengthU : (TPC_VIEW_V == hitType) ? straightLineLengthV : straightLineLengthW);
        float &integratedPathLength((TPC_VIEW_U == hitType) ? integratedPathLengthU : (TPC_VIEW_V == hitType) ? integratedPathLengthV : integratedPathLengthW);

        try
        {
            const TwoDSlidingFitResult slidingFitResult(pCluster, 5, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
            const CartesianVector globalMinLayerPosition(slidingFitResult.GetGlobalMinLayerPosition());
            const CartesianVector globalMaxLayerPosition(slidingFitResult.GetGlobalMaxLayerPosition());
            straightLineLength = (globalMaxLayerPosition - globalMinLayerPosition).GetMagnitude();

            integratedPathLength = 0.f;
            CartesianVector previousFitPosition(globalMinLayerPosition);
            const LayerFitResultMap &layerFitResultMap(slidingFitResult.GetLayerFitResultMap());

            for (const auto &mapEntry : layerFitResultMap)
            {
                CartesianVector thisFitPosition(0.f, 0.f, 0.f);
                slidingFitResult.GetGlobalPosition(mapEntry.second.GetL(), mapEntry.second.GetFitT(), thisFitPosition);
                integratedPathLength += (thisFitPosition - previousFitPosition).GetMagnitude();
                previousFitPosition = thisFitPosition;
            }
        }
        catch (const StatusCodeException &)
        {
        }

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "straightLineLength" + hitTypeLabel, straightLineLength));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "integratedPathLength" + hitTypeLabel, integratedPathLength));

        float &straightLineLength10((TPC_VIEW_U == hitType) ? straightLineLength10U : (TPC_VIEW_V == hitType) ? straightLineLength10V : straightLineLength10W);
        float &integratedPathLength10((TPC_VIEW_U == hitType) ? integratedPathLength10U : (TPC_VIEW_V == hitType) ? integratedPathLength10V : integratedPathLength10W);

        try
        {
            const TwoDSlidingFitResult slidingFitResult(pCluster, 10, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
            const CartesianVector globalMinLayerPosition(slidingFitResult.GetGlobalMinLayerPosition());
            const CartesianVector globalMaxLayerPosition(slidingFitResult.GetGlobalMaxLayerPosition());
            straightLineLength10 = (globalMaxLayerPosition - globalMinLayerPosition).GetMagnitude();

            integratedPathLength10 = 0.f;
            CartesianVector previousFitPosition(globalMinLayerPosition);
            const LayerFitResultMap &layerFitResultMap(slidingFitResult.GetLayerFitResultMap());

            for (const auto &mapEntry : layerFitResultMap)
            {
                CartesianVector thisFitPosition(0.f, 0.f, 0.f);
                slidingFitResult.GetGlobalPosition(mapEntry.second.GetL(), mapEntry.second.GetFitT(), thisFitPosition);
                integratedPathLength10 += (thisFitPosition - previousFitPosition).GetMagnitude();
                previousFitPosition = thisFitPosition;
            }
        }
        catch (const StatusCodeException &)
        {
        }

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "straightLineLength10" + hitTypeLabel, straightLineLength10));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "integratedPathLength10" + hitTypeLabel, integratedPathLength10));

        // shower fit width and gap length
        float &showerFitWidth((TPC_VIEW_U == hitType) ? showerFitWidthU : (TPC_VIEW_V == hitType) ? showerFitWidthV : showerFitWidthW);
        float &showerFitGapLength((TPC_VIEW_U == hitType) ? showerFitGapLengthU : (TPC_VIEW_V == hitType) ? showerFitGapLengthV : showerFitGapLengthW);

        try
        {
            const TwoDSlidingShowerFitResult showerFitResult(pCluster, 10, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
            const LayerFitResultMap &layerFitResultMapS(showerFitResult.GetShowerFitResult().GetLayerFitResultMap());
            const LayerFitResultMap &layerFitResultMapP(showerFitResult.GetPositiveEdgeFitResult().GetLayerFitResultMap());
            const LayerFitResultMap &layerFitResultMapN(showerFitResult.GetNegativeEdgeFitResult().GetLayerFitResultMap());

            if (layerFitResultMapS.size() > 1)
            {
                CartesianVector globalMinLayerPositionOnAxis(0.f, 0.f, 0.f), globalMaxLayerPositionOnAxis(0.f, 0.f, 0.f);
                showerFitResult.GetShowerFitResult().GetGlobalPosition(layerFitResultMapS.begin()->second.GetL(), 0.f, globalMinLayerPositionOnAxis);
                showerFitResult.GetShowerFitResult().GetGlobalPosition(layerFitResultMapS.rbegin()->second.GetL(), 0.f, globalMaxLayerPositionOnAxis);
                const float straightLinePathLength((globalMaxLayerPositionOnAxis - globalMinLayerPositionOnAxis).GetMagnitude());

                if (straightLinePathLength > std::numeric_limits<float>::epsilon())
                {
                    showerFitWidth = 0.f;
                    showerFitGapLength = 0.f;
                    CartesianVector previousLayerPosition(globalMinLayerPositionOnAxis);

                    for (LayerFitResultMap::const_iterator iterS = layerFitResultMapS.begin(); iterS != layerFitResultMapS.end(); ++iterS)
                    {
                        CartesianVector thisLayerPosition(0.f, 0.f, 0.f);
                        showerFitResult.GetShowerFitResult().GetGlobalPosition(iterS->second.GetL(), 0.f, thisLayerPosition);
                        const float thisGapLength((thisLayerPosition - previousLayerPosition).GetMagnitude());

                        if (thisGapLength > showerFitGapLength)
                        {
                            const float minZ(std::min(thisLayerPosition.GetZ(), previousLayerPosition.GetZ()));
                            const float maxZ(std::max(thisLayerPosition.GetZ(), previousLayerPosition.GetZ()));

                            if ((maxZ - minZ) < std::numeric_limits<float>::epsilon())
                                throw StatusCodeException(STATUS_CODE_FAILURE);

                            const float gapZ(LArGeometryHelper::CalculateGapDeltaZ(this->GetPandora(), minZ, maxZ, LArClusterHelper::GetClusterHitType(pCluster)));
                            const float correctedGapLength(thisGapLength * (1.f - gapZ / (maxZ - minZ)));

                            if (correctedGapLength > showerFitGapLength)
                                showerFitGapLength = correctedGapLength;
                        }

                        previousLayerPosition = thisLayerPosition;

                        LayerFitResultMap::const_iterator iterP = layerFitResultMapP.find(iterS->first);
                        LayerFitResultMap::const_iterator iterN = layerFitResultMapN.find(iterS->first);

                        if ((layerFitResultMapP.end() == iterP) || (layerFitResultMapN.end() == iterN))
                            continue;

                        showerFitWidth += std::fabs(iterP->second.GetFitT() - iterN->second.GetFitT());
                    }
                }
            }
        }
        catch (StatusCodeException &)
        {
        }

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "showerFitWidth" + hitTypeLabel, showerFitWidth));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "showerFitGapLength" + hitTypeLabel, showerFitGapLength));

        // nPointsOfContact, nHitsInPrimaryBranches
        int &nPointsOfContact((TPC_VIEW_U == hitType) ? nPointsOfContactU : (TPC_VIEW_V == hitType) ? nPointsOfContactV : nPointsOfContactW);
        int &nHitsInBranches((TPC_VIEW_U == hitType) ? nHitsInBranchesU : (TPC_VIEW_V == hitType) ? nHitsInBranchesV : nHitsInBranchesW);

        const std::string &clusterListName((TPC_VIEW_U == hitType) ? m_clusterListNameU : (TPC_VIEW_V == hitType) ? m_clusterListNameV : m_clusterListNameW);
        const ClusterList *pClusterList(nullptr);

        if (STATUS_CODE_SUCCESS == PandoraContentApi::GetList(*this, clusterListName, pClusterList))
        {
            ClusterVector candidateClusters;
            for (const Cluster *const pCandidateCluster : *pClusterList)
            {
                if ((pCandidateCluster != pCluster) && (pCandidateCluster->GetNCaloHits() > 5))
                    candidateClusters.push_back(pCandidateCluster);
            }

            ClusterUsageMap forwardUsageMap, backwardUsageMap;
            this->FindAssociatedClusters(pCluster, candidateClusters, forwardUsageMap, backwardUsageMap);

            SeedAssociationList seedAssociationList;
            this->IdentifyClusterMerges(ClusterVector(1, pCluster), backwardUsageMap, seedAssociationList);

            for (const Cluster *const pBranchCluster : seedAssociationList.at(pCluster))
            {
                ++nPointsOfContact;
                nHitsInBranches += pBranchCluster->GetNCaloHits();
            }
        }

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPointsOfContact" + hitTypeLabel, nPointsOfContact));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsInBranches" + hitTypeLabel, nHitsInBranches));
    }

    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));

    //--------------------------------------------------------------------------------------------------------------------------------------
    // End tree variables calculations

    return isTrueTrack;
}

//------------------------------------------------------------------------------------------------------------------------------------------

PfoCharacterisationAlgorithm::AssociationType PfoCharacterisationAlgorithm::AreClustersAssociated(const Cluster *const pClusterSeed, const Cluster *const pCluster) const
{
    const float m_nearbyClusterDistance(2.5f);
    const float m_remoteClusterDistance(10.f);

    // Calculate distances of association
    const float sOuter(LArClusterHelper::GetClosestDistance(pClusterSeed->GetCentroid(pClusterSeed->GetOuterPseudoLayer()), pCluster));
    const float cOuter(LArClusterHelper::GetClosestDistance(pCluster->GetCentroid(pCluster->GetOuterPseudoLayer()), pClusterSeed));
    const float sInner(LArClusterHelper::GetClosestDistance(pClusterSeed->GetCentroid(pClusterSeed->GetInnerPseudoLayer()), pCluster));
    const float cInner(LArClusterHelper::GetClosestDistance(pCluster->GetCentroid(pCluster->GetInnerPseudoLayer()), pClusterSeed));

    // Association check 1(a), look for enclosed clusters
    if ((cOuter < m_nearbyClusterDistance && cInner < m_nearbyClusterDistance) &&
        (sInner > m_nearbyClusterDistance) &&
        (sOuter > m_nearbyClusterDistance))
    {
        return STRONG;
    }

    // Association check 1(b), look for overlapping clusters
    if ((cInner < m_nearbyClusterDistance && sOuter < m_nearbyClusterDistance) &&
        (sInner > m_nearbyClusterDistance) &&
        (cOuter > m_nearbyClusterDistance))
    {
        return STRONG;
    }

    if ((cOuter < m_nearbyClusterDistance && sInner < m_nearbyClusterDistance) &&
        (sOuter > m_nearbyClusterDistance) &&
        (cInner > m_nearbyClusterDistance))
    {
        return STRONG;
    }

    // Association check 2, look for branching clusters
    if (((sInner > m_remoteClusterDistance)) &&
        ((sOuter > m_remoteClusterDistance)) &&
        ((cInner < m_nearbyClusterDistance) || (cOuter < m_nearbyClusterDistance)))
    {
        return STANDARD;
    }

    // Association check 3, look any distance below threshold
    if ((sOuter < m_nearbyClusterDistance) || (cOuter < m_nearbyClusterDistance) || (sInner < m_nearbyClusterDistance) || (cInner < m_nearbyClusterDistance))
        return SINGLE_ORDER;

    return NONE;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrackPfoListName", m_trackPfoListName));
    m_inputPfoListNames.push_back(m_trackPfoListName);

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ShowerPfoListName", m_showerPfoListName));
    m_inputPfoListNames.push_back(m_showerPfoListName);

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListNameU", m_clusterListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListNameV", m_clusterListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListNameW", m_clusterListNameW));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "WriteToTree", m_writeToTree));

    if (m_writeToTree)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputTree", m_treeName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputFile", m_fileName));
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
