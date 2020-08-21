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
    //std::cout << "jam" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------    
    
bool ExtensionThroughShowerAlgorithm::FindBestClusterAssociation(const ClusterVector &clusterVector, const SlidingFitResultMapPair &slidingFitResultMapPair,
    const ClusterList *const pClusterList, const bool isHigherXBoundary, ClusterEndpointAssociation &clusterAssociation)
{
    const float nearestTPCBoundaryX(isHigherXBoundary ? m_tpcMaxXEdge : m_tpcMinXEdge);
    
    // THIS ASSUMES THAT WE HAVE SORTED BY DISTANCE FROM TPC BOUNDARY GREATEST->SMALLEST
    for (const Cluster *const pCluster : clusterVector)
    {
        const TwoDSlidingFitResult &microSlidingFitResult(slidingFitResultMapPair.first->at(pCluster));
        const TwoDSlidingFitResult &macroSlidingFitResult(slidingFitResultMapPair.second->at(pCluster));

        const bool isEndUpstream = (std::fabs(microSlidingFitResult.GetGlobalMinLayerPosition().GetX() - nearestTPCBoundaryX) <
                                    std::fabs(microSlidingFitResult.GetGlobalMaxLayerPosition().GetX() - nearestTPCBoundaryX));
        /////////////////
        ClusterList theCluster({pCluster});
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &theCluster, "CONSIDERED CLUSTER", BLACK);
        std::cout << "isEndUpstream: " << isEndUpstream << std::endl;
        ////////////////
            
        CartesianVector clusterMergePoint(0.f, 0.f, 0.f), clusterMergeDirection(0.f, 0.f, 0.f);
        if (!GetClusterMergingCoordinates(microSlidingFitResult, macroSlidingFitResult, macroSlidingFitResult, isEndUpstream, clusterMergePoint, clusterMergeDirection))
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

        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedPoint, "extrapolatedPoint", RED, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &endpointPosition, "endpoint", BLUE, 2);

        // Is there a shower near the cluster endpoint?
        unsigned int showerClusterCount(0);
        //unsigned int showerClusterHitCount(0);


        const unsigned int currentInnerPseudoLayer(pCluster->GetInnerPseudoLayer()), currentOuterPseudoLayer(pCluster->GetOuterPseudoLayer());
        
        for (const Cluster *const pTestCluster : *pClusterList)
        {
            if (pTestCluster == pCluster)
                continue;
                
            const unsigned int testInnerPseudoLayer(pTestCluster->GetInnerPseudoLayer()), testOuterPseudoLayer(pTestCluster->GetOuterPseudoLayer());
            const CartesianVector &innerPoint(pTestCluster->GetCentroid(testInnerPseudoLayer)), &outerPoint(pTestCluster->GetCentroid(testOuterPseudoLayer));              

            const float innerLongitudinalDistance(clusterMergeDirection.GetDotProduct(innerPoint - extrapolatedPoint) * directionFactor);
            const float outerLongitudinalDistance(clusterMergeDirection.GetDotProduct(outerPoint - extrapolatedPoint) * directionFactor);
            
            const float innerTransverseDistance(clusterMergeDirection.GetCrossProduct(innerPoint - extrapolatedPoint).GetMagnitude());
            const float outerTransverseDistance(clusterMergeDirection.GetCrossProduct(outerPoint - extrapolatedPoint).GetMagnitude());
            
            if ((currentInnerPseudoLayer > testInnerPseudoLayer) && (currentOuterPseudoLayer < testOuterPseudoLayer))
            {
                
                std::cout << "CONTAINED SEGMENT OF THE TRACK" << std::endl;
                std::cout << "innerTransverseDistance: " << innerTransverseDistance << std::endl;
                std::cout << "outerTransverseDistance: " << outerTransverseDistance << std::endl;

                ClusterList frog({pTestCluster});
                PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &frog, "CONTAINED", VIOLET);
                
                if ((innerTransverseDistance < 30.f) || (outerTransverseDistance < 30.f))
                {
                    std::cout << "contained segement of the track" << std::endl;
                    showerClusterCount = 0; //showerClusterHitCount = 0;
                    break;
                }
            }


            if (!((innerLongitudinalDistance > 0.f) && (outerLongitudinalDistance > 0.f)))
                continue;
                        
            float closestLongitudinal(0.f), closestTransverse(0.f);
            closestTransverse = (innerTransverseDistance < outerTransverseDistance) ? innerTransverseDistance : outerTransverseDistance;
            closestLongitudinal = (innerTransverseDistance < outerTransverseDistance) ? innerLongitudinalDistance : outerLongitudinalDistance;

            
            if ((closestTransverse < 3.f) && (closestLongitudinal < 20.f))
            {
                
                ClusterList frog({pTestCluster});
                std::string frogString("L: " + std::to_string(closestLongitudinal) + " T: " + std::to_string(closestTransverse));
                PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &frog, frogString, BLUE);    
                
                ++showerClusterCount;
                //showerClusterHitCount += pTestCluster->GetNCaloHits();
            }
        }

        std::cout << "FOUND " << showerClusterCount << " GOOD CLUSTER(S)" << std::endl;
        //std::cout << "WITH " << showerClusterHitCount << " HIT(S)" << std::endl;
        if ((showerClusterCount < 5)) //|| (showerClusterHitCount < 50))
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

    std::cout << "DID NOT FIND AN ASSOCIATION" << std::endl;
    return false;
    
}

//------------------------------------------------------------------------------------------------------------------------------------------  
    
StatusCode ExtensionThroughShowerAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    return TrackExtensionRefinementAlgorithm::ReadSettings(xmlHandle);
}

}  // namespace lar_content
