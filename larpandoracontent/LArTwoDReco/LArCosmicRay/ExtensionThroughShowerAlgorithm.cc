/**
 *  @file   larpandoracontent/LArTwoDReco/LArCosmicRay/ExtensionThroughShowerAlgorithm.cc
 *
 *  @brief  Implementation of the extension through shower class.
 *
 *  $Log: $
 */
#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArControlFlow/MultiPandoraApi.h"

#include "larpandoracontent/LArTwoDReco/LArCosmicRay/ExtensionThroughShowerAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArStitchingHelper.h"

using namespace pandora;

namespace lar_content
{
    
ExtensionThroughShowerAlgorithm::ExtensionThroughShowerAlgorithm()
{
    std::cout << "jam" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------    
    
bool ExtensionThroughShowerAlgorithm::FindBestClusterAssociation(const ClusterVector &clusterVector, const SlidingFitResultMapPair &slidingFitResultMapPair,
    ClusterEndpointAssociation &clusterAssociation, const ClusterList *const pClusterList, const bool isHigherXBoundary)
{
    const float nearestTPCBoundaryX(isHigherXBoundary ? m_tpcMaxXEdge : m_tpcMinXEdge);
    
    // THIS ASSUMES THAT WE HAVE SORTED BY DISTANCE FROM TPC BOUNDARY GREATEST->SMALLEST
    for (const Cluster *const pCluster : clusterVector)
    {
        const TwoDSlidingFitResult &microSlidingFitResult(slidingFitResultMapPair.first->at(pCluster));
        const TwoDSlidingFitResult &macroSlidingFitResult(slidingFitResultMapPair.second->at(pCluster));

        const bool isEndUpstream = (std::fabs(microSlidingFitResult.GetGlobalMinLayerPosition().GetX() - nearestTPCBoundaryX) <
                                    std::fabs(microSlidingFitResult.GetGlobalMaxLayerPosition().GetX() - nearestTPCBoundaryX));

        // ATTN: Match the cluster to itself to get the cluster merging points
        const bool isClusterUpstream(!isEndUpstream);
        
        /////////////////
        ClusterList theCluster({pCluster});
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &theCluster, "CONSIDERED CLUSTER", BLACK);
        std::cout << "isEndUpstream: " << isEndUpstream << std::endl;
        ////////////////
            
        CartesianVector clusterMergePoint(0.f, 0.f, 0.f), clusterMergeDirection(0.f, 0.f, 0.f);
        if (!GetClusterMergingCoordinates(microSlidingFitResult, macroSlidingFitResult, macroSlidingFitResult, isClusterUpstream, clusterMergePoint, clusterMergeDirection))
        {
            std::cout << "CANNOT FIND MERGE POSITION" << std::endl;
            PandoraMonitoringApi::ViewEvent(this->GetPandora());
            continue;
        }

        // WILL NOT CROSS TPC BOUNDARY
        if (std::fabs(clusterMergeDirection.GetX()) < std::numeric_limits<float>::epsilon())
        {
            std::cout << "MERGE DIRECTION HAS NO X COMPONENT" << std::endl;
            PandoraMonitoringApi::ViewEvent(this->GetPandora());
            continue;
        }    

        // Extrapolate track, incase any clustering errors
        const int directionFactor(isEndUpstream ? -1.f : 1.f);
        const CartesianVector &endpointPosition(isEndUpstream ? microSlidingFitResult.GetGlobalMinLayerPosition() : microSlidingFitResult.GetGlobalMaxLayerPosition());
        const float separationDistance((endpointPosition - clusterMergePoint).GetMagnitude());
        const CartesianVector extrapolatedPoint(clusterMergePoint + (clusterMergeDirection * directionFactor * separationDistance));

        // Is there a shower near the cluster endpoint?
        ClusterList goodClusters;
        unsigned int showerClusterCount(0);
        unsigned int showerClusterHitCount(0);
        for (const Cluster *const pTestCluster : *pClusterList)
        {
            if (pTestCluster == pCluster)
                continue;
                
            const unsigned int testInnerPseudoLayer(pTestCluster->GetInnerPseudoLayer()), testOuterPseudoLayer(pTestCluster->GetOuterPseudoLayer());
            const CartesianVector &innerPoint(pTestCluster->GetCentroid(testInnerPseudoLayer)), &outerPoint(pTestCluster->GetCentroid(testOuterPseudoLayer));

            const float innerLongitudinalDistance(clusterMergeDirection.GetDotProduct(innerPoint - extrapolatedPoint) * directionFactor);
            const float outerLongitudinalDistance(clusterMergeDirection.GetDotProduct(innerPoint - extrapolatedPoint) * directionFactor);

            // think about getting rid of this? 
            if (!((innerLongitudinalDistance > 0.f) && (outerLongitudinalDistance > 0.f)))
                continue;

            float closestLongitudinal(0.f), closestTransverse(0.f);
            const float innerTransverseDistance(clusterMergeDirection.GetCrossProduct(innerPoint - extrapolatedPoint).GetMagnitude());
            const float outerTransverseDistance(clusterMergeDirection.GetCrossProduct(outerPoint - extrapolatedPoint).GetMagnitude());

            closestTransverse = (innerTransverseDistance < outerTransverseDistance) ? innerTransverseDistance : outerTransverseDistance;
            closestLongitudinal = (innerTransverseDistance < outerTransverseDistance) ? innerLongitudinalDistance : outerLongitudinalDistance;

            if ((closestTransverse < 3.f) && (closestLongitudinal < 20.f))
            {
                ClusterList frog({pTestCluster});
                std::string frogString("L: " + std::to_string(closestLongitudinal) + " T: " + std::to_string(closestTransverse));
                PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &frog, frogString, BLUE);
                goodClusters.push_back(pTestCluster);
                ++showerClusterCount;
                showerClusterHitCount += pTestCluster->GetNCaloHits();
            }
        }

        std::cout << "FOUND " << showerClusterCount << " GOOD CLUSTER(S)" << std::endl;
        std::cout << "WITH " << showerClusterHitCount << " HIT(S)" << std::endl;
        if (goodClusters.size() < 7 && showerClusterHitCount < 60)
        {
            std::cout << "DOES NOT PASS CRITERIA" << std::endl;
            PandoraMonitoringApi::ViewEvent(this->GetPandora());
            continue;
        }
                
        // ATTN: Temporarily set the other merge point to define region to search for in extrapolate hits
        const float predictedGradient(clusterMergeDirection.GetZ() / clusterMergeDirection.GetX());
        const float predictedIntercept(clusterMergePoint.GetZ() - (predictedGradient * clusterMergePoint.GetX()));
        const CartesianVector extrapolatedEndpointPosition(nearestTPCBoundaryX, 0.f, predictedIntercept + (predictedGradient * nearestTPCBoundaryX));

        // IS THIS EVER NEEDED? - I THINK IT CAN GO WRONG IF THE FIT DOESN'T MAKE SENSE
        if (isEndUpstream ? extrapolatedEndpointPosition.GetZ() > clusterMergePoint.GetZ() : extrapolatedEndpointPosition.GetZ() < clusterMergePoint.GetZ())
        {
            std::cout << "EXTRAPOLATED ENDPOINT IS NOT IN FORWARD DIRECTION" << std::endl;
            PandoraMonitoringApi::ViewEvent(this->GetPandora());
            continue;
        }
                
        clusterAssociation = isEndUpstream ?
            ClusterEndpointAssociation(extrapolatedEndpointPosition, clusterMergeDirection, clusterMergePoint, clusterMergeDirection * (-1.f), pCluster, true) :
            ClusterEndpointAssociation(clusterMergePoint, clusterMergeDirection, extrapolatedEndpointPosition, clusterMergeDirection * (-1.f), pCluster, false);

        ClusterList mainCluster({clusterAssociation.GetMainTrackCluster()});
        const CartesianVector &upstream(clusterAssociation.GetUpstreamMergePoint());
        const CartesianVector &downstream(clusterAssociation.GetDownstreamMergePoint());
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &upstream, "UPSTREAM", VIOLET, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &downstream, "DOWNSTREAM", RED, 2);
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &mainCluster, "CLUSTER", BLACK);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        
        return true;
    }
                
    return false;
    
}

//------------------------------------------------------------------------------------------------------------------------------------------  
    
StatusCode ExtensionThroughShowerAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    return TrackExtensionRefinementAlgorithm::ReadSettings(xmlHandle);
}

}  // namespace lar_content
