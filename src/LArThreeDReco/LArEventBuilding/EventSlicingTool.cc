/**
 *  @file   LArContent/src/LArThreeDReco/LArEventBuilding/EventSlicingTool.cc
 * 
 *  @brief  Implementation of the event slicing tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"
#include "LArHelpers/LArPfoHelper.h"
#include "LArHelpers/LArPointingClusterHelper.h"

#include "LArObjects/LArThreeDSlidingFitResult.h"

#include "LArThreeDReco/LArEventBuilding/EventSlicingTool.h"

using namespace pandora;

namespace lar_content
{

EventSlicingTool::EventSlicingTool() :
    m_halfWindowLayers(20),
    m_minVertexLongitudinalDistance(-7.5f),
    m_maxVertexLongitudinalDistance(60.f),
    m_maxVertexTransverseDistance(10.5f),
    m_vertexAngularAllowance(9.f),
    m_maxClosestApproach(15.f),
    m_maxInterceptDistance(60.f),
    m_maxHitSeparationSquared(30.f * 30.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventSlicingTool::Slice(NeutrinoParentAlgorithm *const pAlgorithm, const HitTypeToNameMap &caloHitListNames,
    const HitTypeToNameMap &clusterListNames, SliceList &sliceList)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << this->GetType() << std::endl;

    ClusterToPfoMap clusterToPfoMap;

    ClusterList trackClusters3D;
    this->GetThreeDClusters(pAlgorithm, m_trackPfoListName, trackClusters3D, clusterToPfoMap);

    ClusterList showerClusters3D;
    this->GetThreeDClusters(pAlgorithm, m_showerPfoListName, showerClusters3D, clusterToPfoMap);

    ClusterSliceList clusterSliceList;
    this->GetClusterSliceList(trackClusters3D, showerClusters3D, clusterSliceList);

    if (clusterSliceList.size() < 2)
    {
        this->CopyAllHitsToSingleSlice(pAlgorithm, caloHitListNames, sliceList);
    }
    else
    {
        ClusterToSliceIndexMap clusterToSliceIndexMap;        
        this->CreateSlices(clusterSliceList, sliceList, clusterToSliceIndexMap);

        ClusterList assignedClusters;
        this->CopyPfoHitsToSlices(clusterToSliceIndexMap, clusterToPfoMap, sliceList, assignedClusters);

        ClusterList remainingClusters;
        this->GetRemainingClusters(pAlgorithm, clusterListNames, assignedClusters, remainingClusters);

        this->AssignRemainingHitsToSlices(remainingClusters, clusterToSliceIndexMap, sliceList);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventSlicingTool::GetThreeDClusters(const Algorithm *const pAlgorithm, const std::string &pfoListName, ClusterList &clusters3D,
    ClusterToPfoMap &clusterToPfoMap) const
{
    const PfoList *pPfoList(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*pAlgorithm, pfoListName, pPfoList));

    if (!pPfoList || pPfoList->empty())
    {
        if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
            std::cout << "EventSlicingTool: unable to find pfo list " << pfoListName << std::endl;

        return;
    }

    for (const ParticleFlowObject *const pPfo : *pPfoList)
    {
        ClusterList pfoClusters3D;
        LArPfoHelper::GetThreeDClusterList(pPfo, pfoClusters3D);

        for (const Cluster *const pCluster3D : pfoClusters3D)
        {
            if (!clusterToPfoMap.insert(ClusterToPfoMap::value_type(pCluster3D, pPfo)).second)
                throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
        }

        clusters3D.insert(pfoClusters3D.begin(), pfoClusters3D.end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventSlicingTool::GetClusterSliceList(const ClusterList &trackClusters3D, const ClusterList &showerClusters3D,
    ClusterSliceList &clusterSliceList) const
{
    ClusterVector sortedTrackClusters3D(trackClusters3D.begin(), trackClusters3D.end());
    std::sort(sortedTrackClusters3D.begin(), sortedTrackClusters3D.end(), LArClusterHelper::SortByNHits);

    ThreeDSlidingFitResultMap slidingFitResultMap;
    const float layerPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));

    for (const Cluster *const pCluster3D : sortedTrackClusters3D)
    {
        try {slidingFitResultMap.insert(ThreeDSlidingFitResultMap::value_type(pCluster3D, ThreeDSlidingFitResult(pCluster3D, m_halfWindowLayers, layerPitch)));}
        catch (StatusCodeException &) {std::cout << "EventSlicingTool: ThreeDSlidingFitResult failure for track cluster." << std::endl;}
    }

    ClusterVector sortedShowerClusters3D(showerClusters3D.begin(), showerClusters3D.end());
    std::sort(sortedShowerClusters3D.begin(), sortedShowerClusters3D.end(), LArClusterHelper::SortByNHits);

    ClusterVector sortedClusters3D;
    sortedClusters3D.insert(sortedClusters3D.end(), sortedTrackClusters3D.begin(), sortedTrackClusters3D.end());
    sortedClusters3D.insert(sortedClusters3D.end(), sortedShowerClusters3D.begin(), sortedShowerClusters3D.end());

    ClusterList usedClusters;

    for (const Cluster *const pCluster3D : sortedClusters3D)
    {
        if (usedClusters.count(pCluster3D))
            continue;

        clusterSliceList.push_back(ClusterVector(1, pCluster3D));
        ClusterVector &clusterSlice(clusterSliceList.back());
        usedClusters.insert(pCluster3D);

        while (this->AddNextPointing(sortedTrackClusters3D, slidingFitResultMap, clusterSlice, usedClusters) ||
            this->AddNextProximity(sortedClusters3D, clusterSlice, usedClusters)) { /* Deliberately empty */ }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EventSlicingTool::AddNextPointing(const ClusterVector &trackCandidates, const ThreeDSlidingFitResultMap &slidingFitResultMap,
    ClusterVector &clusterSlice, ClusterList &usedClusters) const
{
    for (const Cluster *const pCandidateCluster : trackCandidates)
    {
        if (usedClusters.count(pCandidateCluster))
            continue;

        ThreeDSlidingFitResultMap::const_iterator candidateIter = slidingFitResultMap.find(pCandidateCluster);

        if (slidingFitResultMap.end() == candidateIter)
            continue;

        const LArPointingCluster candidatePointingCluster(candidateIter->second);

        for (const Cluster *const pClusterInSlice : clusterSlice)
        {
            ThreeDSlidingFitResultMap::const_iterator inSliceIter = slidingFitResultMap.find(pClusterInSlice);

            if (slidingFitResultMap.end() == inSliceIter)
                continue;

            const LArPointingCluster inSlicePointingCluster(inSliceIter->second);

            if (this->CheckClosestApproach(inSlicePointingCluster, candidatePointingCluster) ||
                this->IsEmission(inSlicePointingCluster, candidatePointingCluster) ||
                this->IsNode(inSlicePointingCluster, candidatePointingCluster))
            {
                if (!usedClusters.insert(pCandidateCluster).second)
                    throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);

                // ATTN Must return here, as have invalidated clusterSlice iterators
                clusterSlice.push_back(pCandidateCluster);
                return true;
            }
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EventSlicingTool::CheckClosestApproach(const LArPointingCluster &cluster1, const LArPointingCluster &cluster2) const
{
    return (this->CheckClosestApproach(cluster1.GetInnerVertex(), cluster2.GetInnerVertex()) ||
        this->CheckClosestApproach(cluster1.GetOuterVertex(), cluster2.GetInnerVertex()) ||
        this->CheckClosestApproach(cluster1.GetInnerVertex(), cluster2.GetOuterVertex()) ||
        this->CheckClosestApproach(cluster1.GetOuterVertex(), cluster2.GetOuterVertex()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EventSlicingTool::CheckClosestApproach(const LArPointingCluster::Vertex &vertex1, const LArPointingCluster::Vertex &vertex2) const
{
    CartesianVector intersectionPoint(0.f, 0.f, 0.f);
    float displacement1(std::numeric_limits<float>::max()), displacement2(std::numeric_limits<float>::max());
    LArPointingClusterHelper::GetIntersection(vertex1, vertex2, intersectionPoint, displacement1, displacement2);

    const CartesianVector approach1(vertex1.GetPosition() + vertex1.GetDirection() * displacement1);
    const CartesianVector approach2(vertex2.GetPosition() + vertex2.GetDirection() * displacement2);
    const float closestApproach((approach1 - approach2).GetMagnitude());

    return ((closestApproach < m_maxClosestApproach) && (std::fabs(displacement1) < m_maxInterceptDistance) && (std::fabs(displacement2) < m_maxInterceptDistance));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EventSlicingTool::IsNode(const LArPointingCluster &cluster1, const LArPointingCluster &cluster2) const
{
    return (LArPointingClusterHelper::IsNode(cluster1.GetInnerVertex().GetPosition(), cluster2.GetInnerVertex(), m_minVertexLongitudinalDistance, m_maxVertexTransverseDistance) ||
        LArPointingClusterHelper::IsNode(cluster1.GetOuterVertex().GetPosition(), cluster2.GetInnerVertex(), m_minVertexLongitudinalDistance, m_maxVertexTransverseDistance) ||
        LArPointingClusterHelper::IsNode(cluster1.GetInnerVertex().GetPosition(), cluster2.GetOuterVertex(), m_minVertexLongitudinalDistance, m_maxVertexTransverseDistance) ||
        LArPointingClusterHelper::IsNode(cluster1.GetOuterVertex().GetPosition(), cluster2.GetOuterVertex(), m_minVertexLongitudinalDistance, m_maxVertexTransverseDistance) ||
        LArPointingClusterHelper::IsNode(cluster2.GetInnerVertex().GetPosition(), cluster1.GetInnerVertex(), m_minVertexLongitudinalDistance, m_maxVertexTransverseDistance) ||
        LArPointingClusterHelper::IsNode(cluster2.GetOuterVertex().GetPosition(), cluster1.GetInnerVertex(), m_minVertexLongitudinalDistance, m_maxVertexTransverseDistance) ||
        LArPointingClusterHelper::IsNode(cluster2.GetInnerVertex().GetPosition(), cluster1.GetOuterVertex(), m_minVertexLongitudinalDistance, m_maxVertexTransverseDistance) ||
        LArPointingClusterHelper::IsNode(cluster2.GetOuterVertex().GetPosition(), cluster1.GetOuterVertex(), m_minVertexLongitudinalDistance, m_maxVertexTransverseDistance));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EventSlicingTool::IsEmission(const LArPointingCluster &cluster1, const LArPointingCluster &cluster2) const
{
    return (LArPointingClusterHelper::IsEmission(cluster1.GetInnerVertex().GetPosition(), cluster2.GetInnerVertex(), m_minVertexLongitudinalDistance, m_maxVertexLongitudinalDistance, m_maxVertexTransverseDistance, m_vertexAngularAllowance) ||
        LArPointingClusterHelper::IsEmission(cluster1.GetOuterVertex().GetPosition(), cluster2.GetInnerVertex(), m_minVertexLongitudinalDistance, m_maxVertexLongitudinalDistance, m_maxVertexTransverseDistance, m_vertexAngularAllowance) ||
        LArPointingClusterHelper::IsEmission(cluster1.GetInnerVertex().GetPosition(), cluster2.GetOuterVertex(), m_minVertexLongitudinalDistance, m_maxVertexLongitudinalDistance, m_maxVertexTransverseDistance, m_vertexAngularAllowance) ||
        LArPointingClusterHelper::IsEmission(cluster1.GetOuterVertex().GetPosition(), cluster2.GetOuterVertex(), m_minVertexLongitudinalDistance, m_maxVertexLongitudinalDistance, m_maxVertexTransverseDistance, m_vertexAngularAllowance) ||
        LArPointingClusterHelper::IsEmission(cluster2.GetInnerVertex().GetPosition(), cluster1.GetInnerVertex(), m_minVertexLongitudinalDistance, m_maxVertexLongitudinalDistance, m_maxVertexTransverseDistance, m_vertexAngularAllowance) ||
        LArPointingClusterHelper::IsEmission(cluster2.GetOuterVertex().GetPosition(), cluster1.GetInnerVertex(), m_minVertexLongitudinalDistance, m_maxVertexLongitudinalDistance, m_maxVertexTransverseDistance, m_vertexAngularAllowance) ||
        LArPointingClusterHelper::IsEmission(cluster2.GetInnerVertex().GetPosition(), cluster1.GetOuterVertex(), m_minVertexLongitudinalDistance, m_maxVertexLongitudinalDistance, m_maxVertexTransverseDistance, m_vertexAngularAllowance) ||
        LArPointingClusterHelper::IsEmission(cluster2.GetOuterVertex().GetPosition(), cluster1.GetOuterVertex(), m_minVertexLongitudinalDistance, m_maxVertexLongitudinalDistance, m_maxVertexTransverseDistance, m_vertexAngularAllowance));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EventSlicingTool::AddNextProximity(const ClusterVector &clusterCandidates, ClusterVector &clusterSlice, ClusterList &usedClusters) const
{
    for (const Cluster *const pCandidateCluster : clusterCandidates)
    {
        if (usedClusters.count(pCandidateCluster))
            continue;

        for (const Cluster *const pClusterInSlice : clusterSlice)
        {
            if (this->CheckHitSeparation(pCandidateCluster, pClusterInSlice))
            {
                if (!usedClusters.insert(pCandidateCluster).second)
                    throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);

                // ATTN Must return here, as have invalidated clusterSlice iterators
                clusterSlice.push_back(pCandidateCluster);
                return true;
            }
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EventSlicingTool::CheckHitSeparation(const Cluster *const pCluster1, const Cluster *const pCluster2) const
{
    for (const auto &orderedList1 : pCluster1->GetOrderedCaloHitList())
    {
        for (const CaloHit *const pCaloHit1 : *(orderedList1.second))
        {
            const CartesianVector &positionVector1(pCaloHit1->GetPositionVector());

            for (const auto &orderedList2 : pCluster2->GetOrderedCaloHitList())
            {
                for (const CaloHit *const pCaloHit2 : *(orderedList2.second))
                {
                    const CartesianVector &positionVector2(pCaloHit2->GetPositionVector());

                    if ((positionVector1 - positionVector2).GetMagnitudeSquared() < m_maxHitSeparationSquared)
                        return true;
                }
            }
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventSlicingTool::CopyAllHitsToSingleSlice(const Algorithm *const pAlgorithm, const HitTypeToNameMap &caloHitListNames,
    SliceList &sliceList) const
{
    if (!sliceList.empty())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const CaloHitList *pCaloHitListU(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*pAlgorithm,
        caloHitListNames.at(TPC_VIEW_U), pCaloHitListU));

    const CaloHitList *pCaloHitListV(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*pAlgorithm,
        caloHitListNames.at(TPC_VIEW_V), pCaloHitListV));

    const CaloHitList *pCaloHitListW(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*pAlgorithm,
        caloHitListNames.at(TPC_VIEW_W), pCaloHitListW));

    if (pCaloHitListU || pCaloHitListV || pCaloHitListW)
    {
        sliceList.push_back(NeutrinoParentAlgorithm::Slice());
        NeutrinoParentAlgorithm::Slice &slice(sliceList.at(0));

        if (pCaloHitListU) slice.m_caloHitListU = *pCaloHitListU;
        if (pCaloHitListV) slice.m_caloHitListV = *pCaloHitListV;
        if (pCaloHitListW) slice.m_caloHitListW = *pCaloHitListW;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventSlicingTool::CreateSlices(const ClusterSliceList &clusterSliceList, SliceList &sliceList, ClusterToSliceIndexMap &clusterToSliceIndexMap) const
{
    unsigned int index(0);

    for (const ClusterVector &clusterList : clusterSliceList)
    {
        for (const Cluster *const pCluster3D : clusterList)
        {
            if (TPC_3D != LArClusterHelper::GetClusterHitType(pCluster3D))
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

            if (!clusterToSliceIndexMap.insert(ClusterToSliceIndexMap::value_type(pCluster3D, index)).second)
                throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
        }

        sliceList.push_back(NeutrinoParentAlgorithm::Slice());
        ++index;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventSlicingTool::CopyPfoHitsToSlices(const ClusterToSliceIndexMap &clusterToSliceIndexMap, const ClusterToPfoMap &clusterToPfoMap,
    SliceList &sliceList, ClusterList &assignedClusters) const
{
    for (const ClusterToSliceIndexMap::value_type &mapValue : clusterToSliceIndexMap)
    {
        const Cluster *const pCluster3D(mapValue.first);
        const unsigned int index(mapValue.second);

        const Pfo *const pPfo(clusterToPfoMap.at(pCluster3D));
        NeutrinoParentAlgorithm::Slice &slice(sliceList.at(index));

        ClusterList clusters2D;
        LArPfoHelper::GetTwoDClusterList(pPfo, clusters2D);

        for (const Cluster *const pCluster2D : clusters2D)
        {
            const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster2D));

            if ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType))
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

            CaloHitList &targetList((TPC_VIEW_U == hitType) ? slice.m_caloHitListU : (TPC_VIEW_V == hitType) ? slice.m_caloHitListV : slice.m_caloHitListW);

            pCluster2D->GetOrderedCaloHitList().GetCaloHitList(targetList);
            targetList.insert(pCluster2D->GetIsolatedCaloHitList().begin(), pCluster2D->GetIsolatedCaloHitList().end());

            if (!assignedClusters.insert(pCluster2D).second)
                throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventSlicingTool::GetRemainingClusters(const Algorithm *const pAlgorithm, const HitTypeToNameMap &clusterListNames,
    const ClusterList &assignedClusters, ClusterList &remainingClusters) const
{
    this->GetRemainingClusters(pAlgorithm, clusterListNames.at(TPC_VIEW_U), assignedClusters, remainingClusters);
    this->GetRemainingClusters(pAlgorithm, clusterListNames.at(TPC_VIEW_V), assignedClusters, remainingClusters);
    this->GetRemainingClusters(pAlgorithm, clusterListNames.at(TPC_VIEW_W), assignedClusters, remainingClusters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventSlicingTool::GetRemainingClusters(const Algorithm *const pAlgorithm, const std::string &clusterListName,
    const ClusterList &assignedClusters, ClusterList &remainingClusters) const
{
    const ClusterList *pClusterList(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*pAlgorithm, clusterListName, pClusterList));

    if (!pClusterList || pClusterList->empty())
    {
        if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
            std::cout << "EventSlicingTool: unable to find cluster list " << clusterListName << std::endl;

        return;
    }

    for (const Cluster *const pCluster2D : *pClusterList)
    {
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster2D));

        if ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        if (assignedClusters.count(pCluster2D))
            continue;

        if (!remainingClusters.insert(pCluster2D).second)
            throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventSlicingTool::AssignRemainingHitsToSlices(const ClusterList &remainingClusters, const ClusterToSliceIndexMap &clusterToSliceIndexMap,
    SliceList &sliceList) const
{
    for (const Cluster *const pCluster2D : remainingClusters)
    {
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster2D));

        if ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        const Cluster *pBestCluster3D(nullptr);

        for (const ClusterToSliceIndexMap::value_type &mapValue : clusterToSliceIndexMap)
        {
            const Cluster *const pCluster3D(mapValue.first);

            // TODO associate pCluster2D to pCluster3D using some metrics
            // Could work using 3D clusters, or could use existing contents of slice associated with pCluster3D
            if (false)
                pBestCluster3D = pCluster3D;
        }

        if (pBestCluster3D)
        {
            const unsigned int index(clusterToSliceIndexMap.at(pBestCluster3D));
            NeutrinoParentAlgorithm::Slice &slice(sliceList.at(index));
            CaloHitList &targetList((TPC_VIEW_U == hitType) ? slice.m_caloHitListU : (TPC_VIEW_V == hitType) ? slice.m_caloHitListV : slice.m_caloHitListW);

            pCluster2D->GetOrderedCaloHitList().GetCaloHitList(targetList);
            targetList.insert(pCluster2D->GetIsolatedCaloHitList().begin(), pCluster2D->GetIsolatedCaloHitList().end());
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventSlicingTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "TrackPfoListName", m_trackPfoListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowerPfoListName", m_showerPfoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitHalfWindow", m_halfWindowLayers));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinVertexLongitudinalDistance", m_minVertexLongitudinalDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxVertexLongitudinalDistance", m_maxVertexLongitudinalDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxVertexTransverseDistance", m_maxVertexTransverseDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VertexAngularAllowance", m_vertexAngularAllowance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxClosestApproach", m_maxClosestApproach));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxInterceptDistance", m_maxInterceptDistance));

    float maxHitSeparation = std::sqrt(m_maxHitSeparationSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxHitSeparation", maxHitSeparation));
    m_maxHitSeparationSquared = maxHitSeparation * maxHitSeparation;

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
