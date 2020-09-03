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
    m_maxShowerBranchTransverseDistance(3.f),
    m_maxShowerBranchLongitudinalDistance(20.f),
    m_thresholdShowerClusterCount(5)    
{
}

//------------------------------------------------------------------------------------------------------------------------------------------    

bool ExtensionThroughShowerAlgorithm::DoesPassCriteria(const TwoDSlidingFitResult &microSlidingFitResult, const CartesianVector &clusterMergeDirection,
    const bool isEndUpstream, const ClusterList *const pClusterList, CartesianVector &clusterMergePoint) const
{
    const Cluster *const pCurrentCluster(microSlidingFitResult.GetCluster());
    
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

    return (showerClusterCount == m_thresholdShowerClusterCount);
}    

//------------------------------------------------------------------------------------------------------------------------------------------  
 
StatusCode ExtensionThroughShowerAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxShowerBranchTransverseDistance", m_maxShowerBranchTransverseDistance));    

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxShowerBranchLongitudinalDistance", m_maxShowerBranchLongitudinalDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ThresholdShowerClusterCount", m_thresholdShowerClusterCount));        
    
    return TrackExtensionRefinementAlgorithm::ReadSettings(xmlHandle);
}

}  // namespace lar_content
