/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/OvershootSplittingAlgorithm.cc
 *
 *  @brief  Implementation of the overshoot splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"

#include "larpandoracontent/LArObjects/LArPointingCluster.h"

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/OvershootSplittingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

OvershootSplittingAlgorithm::OvershootSplittingAlgorithm() :
    TwoDSlidingFitMultiSplitAlgorithm(),
    m_minClusterLength(5.f),
    m_maxClusterSeparation(5.f),
    m_minVertexDisplacement(2.f),
    m_maxIntersectDisplacement(1.5f),
    m_minSplitDisplacement(10.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OvershootSplittingAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pCluster = *iter;

        if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLength * m_minClusterLength)
            continue;

        clusterVector.push_back(pCluster);
    }

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OvershootSplittingAlgorithm::FindBestSplitPositions(const TwoDSlidingFitResultMap &slidingFitResultMap, ClusterPositionMap &clusterSplittingMap) const
{
    // Use sliding fit results to build a list of intersection points
    ClusterPositionMap clusterIntersectionMap;
    this->BuildIntersectionMap(slidingFitResultMap, clusterIntersectionMap);

    // Sort intersection points according to their position along the sliding fit
    ClusterPositionMap sortedIntersectionMap;
    this->BuildSortedIntersectionMap(slidingFitResultMap, clusterIntersectionMap, sortedIntersectionMap);

    // Use intersection points to decide where/if to split cluster
    this->PopulateSplitPositionMap(sortedIntersectionMap, clusterSplittingMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OvershootSplittingAlgorithm::BuildIntersectionMap(const TwoDSlidingFitResultMap &slidingFitResultMap, ClusterPositionMap &clusterIntersectionMap) const
{
    ClusterList clusterList;
    for (const auto &mapEntry : slidingFitResultMap)
        clusterList.push_back(mapEntry.first);
    clusterList.sort(LArClusterHelper::SortByNHits);

    for (const Cluster *const pCluster1 : clusterList)
    {
        const TwoDSlidingFitResult &slidingFitResult1(slidingFitResultMap.at(pCluster1));

        for (const Cluster *const pCluster2 : clusterList)
        {
            if (pCluster1 == pCluster2)
                continue;

            const TwoDSlidingFitResult &slidingFitResult2(slidingFitResultMap.at(pCluster2));

            try
            {
                const LArPointingCluster pointingCluster(slidingFitResult2);

                // Project pointing cluster onto target cluster
                const CartesianVector innerPosition(pointingCluster.GetInnerVertex().GetPosition());
                const CartesianVector outerPosition(pointingCluster.GetOuterVertex().GetPosition());
                const float innerDisplacement(LArClusterHelper::GetClosestDistance(innerPosition, pCluster1));
                const float outerDisplacement(LArClusterHelper::GetClosestDistance(outerPosition, pCluster1));
                const bool useInner((innerDisplacement < outerDisplacement) ? true : false);

                const LArPointingCluster::Vertex &clusterVertex = (useInner ? pointingCluster.GetInnerVertex() : pointingCluster.GetOuterVertex());

                float rL2(0.f), rT2(0.f);
                CartesianVector intersectPosition2(0.f, 0.f, 0.f);

                try
                {
                    LArPointingClusterHelper::GetIntersection(clusterVertex, pCluster1, intersectPosition2, rL2, rT2);
                }
                catch (const StatusCodeException &)
                {
                    continue;
                }

                if (rL2 < -m_maxIntersectDisplacement || rL2 > m_maxClusterSeparation)
                    continue;

                // Find projected position and direction on target cluster
                float rL1(0.f), rT1(0.f);
                CartesianVector projectedPosition1(0.f, 0.f, 0.f), projectedDirection1(0.f, 0.f, 0.f);
                slidingFitResult1.GetLocalPosition(intersectPosition2, rL1, rT1);

                const StatusCode statusCodePosition(slidingFitResult1.GetGlobalFitPosition(rL1, projectedPosition1));
                if (STATUS_CODE_SUCCESS != statusCodePosition)
                    throw pandora::StatusCodeException(statusCodePosition);

                const StatusCode statusCodeDirection(slidingFitResult1.GetGlobalFitDirection(rL1, projectedDirection1));
                if (STATUS_CODE_SUCCESS != statusCodeDirection)
                    throw pandora::StatusCodeException(statusCodeDirection);

                const CartesianVector projectedPosition2(clusterVertex.GetPosition());
                const CartesianVector projectedDirection2(clusterVertex.GetDirection());

                // Find intersection of pointing cluster and target cluster
                float firstDisplacement(0.f), secondDisplacement(0.f);
                CartesianVector intersectPosition1(0.f, 0.f, 0.f);

                try
                {
                    LArPointingClusterHelper::GetIntersection(projectedPosition1, projectedDirection1, projectedPosition2,
                        projectedDirection2, intersectPosition1, firstDisplacement, secondDisplacement);
                }
                catch (const StatusCodeException &)
                {
                    continue;
                }

                // Store intersections if they're sufficiently far along the cluster trajectory
                const float closestDisplacement1(LArClusterHelper::GetClosestDistance(intersectPosition1, pCluster1));
                const float closestDisplacement2(LArClusterHelper::GetClosestDistance(intersectPosition1, pCluster2));

                if (std::max(closestDisplacement1, closestDisplacement2) > m_maxClusterSeparation)
                    continue;

                const CartesianVector minPosition(slidingFitResult1.GetGlobalMinLayerPosition());
                const CartesianVector maxPosition(slidingFitResult1.GetGlobalMaxLayerPosition());
                const float lengthSquared((maxPosition - minPosition).GetMagnitudeSquared());

                const float minDisplacementSquared((minPosition - intersectPosition1).GetMagnitudeSquared());
                const float maxDisplacementSquared((maxPosition - intersectPosition1).GetMagnitudeSquared());

                if (std::min(minDisplacementSquared, maxDisplacementSquared) < (m_minVertexDisplacement * m_minVertexDisplacement) ||
                    std::max(minDisplacementSquared, maxDisplacementSquared) > lengthSquared)
                    continue;

                clusterIntersectionMap[pCluster1].push_back(intersectPosition1);
            }
            catch (StatusCodeException &statusCodeException)
            {
                if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                    throw statusCodeException;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OvershootSplittingAlgorithm::BuildSortedIntersectionMap(const TwoDSlidingFitResultMap &slidingFitResultMap,
    const ClusterPositionMap &clusterIntersectionMap, ClusterPositionMap &sortedIntersectionMap) const
{
    ClusterList clusterList;
    for (const auto &mapEntry : clusterIntersectionMap)
        clusterList.push_back(mapEntry.first);
    clusterList.sort(LArClusterHelper::SortByNHits);

    for (const Cluster *const pCluster : clusterList)
    {
        const CartesianPointVector &inputPositionVector(clusterIntersectionMap.at(pCluster));

        if (inputPositionVector.empty())
            continue;

        TwoDSlidingFitResultMap::const_iterator sIter = slidingFitResultMap.find(pCluster);
        if (slidingFitResultMap.end() == sIter)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        const TwoDSlidingFitResult &slidingFitResult = sIter->second;

        MyTrajectoryPointList trajectoryPointList;
        for (CartesianPointVector::const_iterator pIter = inputPositionVector.begin(), pIterEnd = inputPositionVector.end(); pIter != pIterEnd; ++pIter)
        {
            const CartesianVector &position = *pIter;
            float rL(0.f), rT(0.f);
            slidingFitResult.GetLocalPosition(position, rL, rT);
            trajectoryPointList.push_back(MyTrajectoryPoint(rL, position));
        }

        std::sort(trajectoryPointList.begin(), trajectoryPointList.end(), OvershootSplittingAlgorithm::SortByHitProjection);

        if (trajectoryPointList.empty())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        for (MyTrajectoryPointList::const_iterator qIter = trajectoryPointList.begin(), qIterEnd = trajectoryPointList.end(); qIter != qIterEnd; ++qIter)
        {
            const CartesianVector &clusterPosition = qIter->second;
            sortedIntersectionMap[pCluster].push_back(clusterPosition);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void OvershootSplittingAlgorithm::PopulateSplitPositionMap(const ClusterPositionMap &clusterIntersectionMap, ClusterPositionMap &clusterSplittingMap) const
{
    ClusterList clusterList;
    for (const auto &mapEntry : clusterIntersectionMap)
        clusterList.push_back(mapEntry.first);
    clusterList.sort(LArClusterHelper::SortByNHits);

    for (const Cluster *const pCluster : clusterList)
    {
        const CartesianPointVector &inputPositionVector(clusterIntersectionMap.at(pCluster));

        if (inputPositionVector.empty())
            continue;

        // Select pairs of positions within a given separation, and calculate their average position
        MyTrajectoryPointList candidatePositionList;

        bool foundPrevPosition(false);
        CartesianVector prevPosition(0.f, 0.f, 0.f);

        for (CartesianPointVector::const_iterator pIter = inputPositionVector.begin(), pIterEnd = inputPositionVector.end(); pIter != pIterEnd; ++pIter)
        {
            const CartesianVector &nextPosition = *pIter;

            if (foundPrevPosition)
            {
                const CartesianVector averagePosition((nextPosition + prevPosition) * 0.5f);
                const float displacementSquared((nextPosition - prevPosition).GetMagnitudeSquared());

                if (displacementSquared < m_maxIntersectDisplacement * m_maxIntersectDisplacement)
                    candidatePositionList.push_back(MyTrajectoryPoint(displacementSquared, averagePosition));
            }

            prevPosition = nextPosition;
            foundPrevPosition = true;
        }

        if (candidatePositionList.empty())
            continue;

        std::sort(candidatePositionList.begin(), candidatePositionList.end(), OvershootSplittingAlgorithm::SortByHitProjection);

        // Use the average positions of the closest pairs of points as the split position
        bool foundPrevCandidate(false);
        CartesianVector prevCandidate(0.f, 0.f, 0.f);

        for (MyTrajectoryPointList::const_iterator pIter = candidatePositionList.begin(), pIterEnd = candidatePositionList.end();
             pIter != pIterEnd; ++pIter)
        {
            const CartesianVector &nextCandidate = pIter->second;

            if (foundPrevCandidate)
            {
                if ((nextCandidate - prevCandidate).GetMagnitudeSquared() < m_minSplitDisplacement * m_minSplitDisplacement)
                    continue;
            }

            clusterSplittingMap[pCluster].push_back(nextCandidate);
            prevCandidate = nextCandidate;
            foundPrevCandidate = true;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool OvershootSplittingAlgorithm::SortByHitProjection(const MyTrajectoryPoint &lhs, const MyTrajectoryPoint &rhs)
{
    if (lhs.first != rhs.first)
        return (lhs.first < rhs.first);

    return (lhs.second.GetMagnitudeSquared() > rhs.second.GetMagnitudeSquared());
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode OvershootSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterLength", m_minClusterLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxClusterSeparation", m_maxClusterSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinVertexDisplacement", m_minVertexDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxIntersectDisplacement", m_maxIntersectDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinSplitDisplacement", m_minSplitDisplacement));

    return TwoDSlidingFitMultiSplitAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
