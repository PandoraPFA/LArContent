/**
 *  @file   larpandoracontent/LArTwoDReco/LArCosmicRay/ExtensionThroughShowerAlgorithm.cc
 *
 *  @brief  Implementation of the extension through shower class.
 *
 *  $Log: $
 */
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArTwoDReco/LArCosmicRay/ExtensionThroughShowerAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ExtensionThroughShowerAlgorithm::ExtensionThroughShowerAlgorithm() :
    m_maxTrackDistanceToShowerBranch(5.f),
    m_maxShowerBranchTransverseDistance(3.f),
    m_maxShowerBranchLongitudinalDistance(20.f),
    m_thresholdShowerClusterCount(5)    
{
}

//------------------------------------------------------------------------------------------------------------------------------------------    
    
bool ExtensionThroughShowerAlgorithm::FindBestClusterAssociation(const ClusterVector &clusterVector, const SlidingFitResultMapPair &slidingFitResultMapPair,
    const ClusterList *const pClusterList, const bool isHigherXBoundary, ClusterEndpointAssociation &clusterAssociation)
{
    const float nearestTPCBoundaryX(isHigherXBoundary ? m_tpcMaxXEdge : m_tpcMinXEdge);
    
    // ATTN: This assumes that clusterVector has been sorted furthest cluster to TPC boundary -> closest
    for (const Cluster *const pCurrentCluster : clusterVector)
    {
        const TwoDSlidingFitResult &microSlidingFitResult(slidingFitResultMapPair.first->at(pCurrentCluster));
        const TwoDSlidingFitResult &macroSlidingFitResult(slidingFitResultMapPair.second->at(pCurrentCluster));

        const bool isEndUpstream = (std::fabs(microSlidingFitResult.GetGlobalMinLayerPosition().GetX() - nearestTPCBoundaryX) <
                                    std::fabs(microSlidingFitResult.GetGlobalMaxLayerPosition().GetX() - nearestTPCBoundaryX));
            
        CartesianVector clusterMergePoint(0.f, 0.f, 0.f), clusterMergeDirection(0.f, 0.f, 0.f);
        if (!GetClusterMergingCoordinates(microSlidingFitResult, macroSlidingFitResult, macroSlidingFitResult, isEndUpstream, clusterMergePoint, clusterMergeDirection))
            continue;

        // Reject clusters that do not cross TPC boundary
        if (std::fabs(clusterMergeDirection.GetX()) < std::numeric_limits<float>::epsilon())
            continue;

        // Reject clusters that are not infront of a shower region
        if(!this->IsClusterEndpointInFrontOfShower(microSlidingFitResult, clusterMergePoint, clusterMergeDirection, isEndUpstream, pClusterList))
            continue;
                
        // ATTN: Temporarily set the other merge point to define extrapolate hits search region
        const float predictedGradient(clusterMergeDirection.GetZ() / clusterMergeDirection.GetX());
        const float predictedIntercept(clusterMergePoint.GetZ() - (predictedGradient * clusterMergePoint.GetX()));
        const CartesianVector extrapolatedHitsEndpoint(nearestTPCBoundaryX, 0.f, predictedIntercept + (predictedGradient * nearestTPCBoundaryX));

        if (isEndUpstream ? clusterMergePoint.GetZ() < extrapolatedHitsEndpoint.GetZ() : clusterMergePoint.GetZ() > extrapolatedHitsEndpoint.GetZ())
            continue;
                
        clusterAssociation = isEndUpstream ?
            ClusterEndpointAssociation(extrapolatedHitsEndpoint, clusterMergeDirection, clusterMergePoint, clusterMergeDirection * (-1.f), pCurrentCluster, true) :
            ClusterEndpointAssociation(clusterMergePoint, clusterMergeDirection, extrapolatedHitsEndpoint, clusterMergeDirection * (-1.f), pCurrentCluster, false);

        return true;
    }

    return false;    
}

//------------------------------------------------------------------------------------------------------------------------------------------  
    
bool ExtensionThroughShowerAlgorithm::IsClusterEndpointInFrontOfShower(const TwoDSlidingFitResult &microSlidingFitResult, const CartesianVector &clusterMergePoint,
    const CartesianVector &clusterMergeDirection, const bool isEndUpstream, const ClusterList *const pClusterList)
{
    const Cluster *const pCurrentCluster(microSlidingFitResult.GetCluster());
    const unsigned int currentInnerPseudoLayer(pCurrentCluster->GetInnerPseudoLayer()), currentOuterPseudoLayer(pCurrentCluster->GetOuterPseudoLayer());
    
    // Find the track endpoint accounting for possible clustering errors
    const int directionFactor(isEndUpstream ? -1 : 1);
    const CartesianVector &endpointPosition(isEndUpstream ? microSlidingFitResult.GetGlobalMinLayerPosition() : microSlidingFitResult.GetGlobalMaxLayerPosition());
    const float separationDistance((endpointPosition - clusterMergePoint).GetMagnitude());
    const CartesianVector extrapolatedEndpoint(clusterMergePoint + (clusterMergeDirection * directionFactor * separationDistance));

    unsigned int showerClusterCount(0);
    for (const Cluster *const pTestCluster : *pClusterList)
    {
        if (pTestCluster == pCurrentCluster)
            continue;
                
        const unsigned int testInnerPseudoLayer(pTestCluster->GetInnerPseudoLayer()), testOuterPseudoLayer(pTestCluster->GetOuterPseudoLayer());
        const CartesianVector &innerPoint(pTestCluster->GetCentroid(testInnerPseudoLayer)), &outerPoint(pTestCluster->GetCentroid(testOuterPseudoLayer));

        // ATTN: If pCurrent cluster is contained within pTestCluster it is likely to be part of the shower
        if ((currentInnerPseudoLayer > testInnerPseudoLayer) && (currentOuterPseudoLayer < testOuterPseudoLayer))
        {
            if (LArClusterHelper::GetClosestDistance(pCurrentCluster, pTestCluster) < m_maxTrackDistanceToShowerBranch)
                return false;
        }

        if (showerClusterCount != m_thresholdShowerClusterCount)
        {
            const float innerTransverseDistance(clusterMergeDirection.GetCrossProduct(innerPoint - extrapolatedEndpoint).GetMagnitude());
            const float outerTransverseDistance(clusterMergeDirection.GetCrossProduct(outerPoint - extrapolatedEndpoint).GetMagnitude());    
            const float innerLongitudinalDistance(clusterMergeDirection.GetDotProduct(innerPoint - extrapolatedEndpoint));
            const float outerLongitudinalDistance(clusterMergeDirection.GetDotProduct(outerPoint - extrapolatedEndpoint));

            // ATTN: Only consider showers that are fully forwards of the extrapolatedEndpoint
            const bool isForward((innerLongitudinalDistance > 0.f) && (outerLongitudinalDistance > 0.f));
            if (!isForward)
                continue;
                        
            const float closestTransverse(std::min(innerTransverseDistance, outerTransverseDistance));
            const float closestLongitudinal(innerTransverseDistance < outerTransverseDistance ? innerLongitudinalDistance : outerLongitudinalDistance);

            if ((closestTransverse < m_maxShowerBranchTransverseDistance) && (closestLongitudinal < m_maxShowerBranchLongitudinalDistance))
                ++showerClusterCount;
        }
    }

    return (showerClusterCount == m_thresholdShowerClusterCount ? true : false);
}    

//------------------------------------------------------------------------------------------------------------------------------------------  
 
StatusCode ExtensionThroughShowerAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTrackDistanceToShowerBranch", m_maxTrackDistanceToShowerBranch));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxShowerBranchTransverseDistance", m_maxShowerBranchTransverseDistance));    

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxShowerBranchLongitudinalDistance", m_maxShowerBranchLongitudinalDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ThresholdShowerClusterCount", m_thresholdShowerClusterCount));        
    
    return TrackExtensionRefinementAlgorithm::ReadSettings(xmlHandle);
}

}  // namespace lar_content
