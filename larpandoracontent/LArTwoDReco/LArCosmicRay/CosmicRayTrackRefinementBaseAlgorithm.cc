/**
 *  @file   larpandoracontent/LArTwoDReco/LArCosmicRay/CosmicRayTrackRefinementBaseAlgorithm.cc
 *
 *  @brief  Implementation of the cosmic ray track refinement base class.
 *
 *  $Log: $
 */
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArTwoDReco/LArCosmicRay/CosmicRayTrackRefinementBaseAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

using namespace pandora;

namespace lar_content
{

CosmicRayTrackRefinementBaseAlgorithm::CosmicRayTrackRefinementBaseAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------    

void CosmicRayTrackRefinementBaseAlgorithm::InitialiseSlidingFitResultMaps(const ClusterVector &clusterVector, TwoDSlidingFitResultMap &microSlidingFitResultMap,
    TwoDSlidingFitResultMap &macroSlidingFitResultMap) const
{
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const int m_slidingFitWindow(10);
    const int m_macroSlidingFitWindow(1000);
    
    for (const Cluster *const pCluster : clusterVector)
    {
        try
        {
            TwoDSlidingFitResult microSlidingFitResult(pCluster, m_slidingFitWindow, slidingFitPitch);
            TwoDSlidingFitResult macroSlidingFitResult(pCluster, m_macroSlidingFitWindow, slidingFitPitch);

            (void) microSlidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pCluster, microSlidingFitResult));
            (void) macroSlidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pCluster, macroSlidingFitResult));
        }
        catch (StatusCodeException &) {}
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CosmicRayTrackRefinementBaseAlgorithm::GetClusterMergingCoordinates(const TwoDSlidingFitResult &clusterMicroFitResult, const TwoDSlidingFitResult &clusterMacroFitResult,
    const TwoDSlidingFitResult &associatedMacroFitResult, const bool isUpstream, CartesianVector &clusterMergePosition, CartesianVector &clusterMergeDirection) const
{
    CartesianVector clusterAverageDirection(0.f, 0.f, 0.f), associatedAverageDirection(0.f, 0.f, 0.f);
    clusterMacroFitResult.GetGlobalDirection(clusterMacroFitResult.GetLayerFitResultMap().begin()->second.GetGradient(), clusterAverageDirection);
    associatedMacroFitResult.GetGlobalDirection(associatedMacroFitResult.GetLayerFitResultMap().begin()->second.GetGradient(), associatedAverageDirection);    

    //ISOBEL CHANGE THIS
    float m_stableRegionClusterFraction(0.05);
    float m_mergePointMinCosAngleDeviation(0.995);
    
    const LayerFitResultMap &clusterMicroLayerFitResultMap(clusterMicroFitResult.GetLayerFitResultMap());
    const int startLayer(isUpstream ? clusterMicroFitResult.GetMaxLayer() : clusterMicroFitResult.GetMinLayer());
    const int endLayer(isUpstream ? clusterMicroFitResult.GetMinLayer() : clusterMicroFitResult.GetMaxLayer());
    const int loopTerminationLayer(endLayer + (isUpstream ? -1 : 1));
    const int step(isUpstream ? -1 : 1);

    std::cout << "DO NOT FORGET THAT YOU HAVE THE DIRECTIONS FACING ONE ANOTHER" << std::endl;
    std::cout << "DO NOT FORGET THAT YOU HAVE CHANGED THIS SO THAT WE USE THE LOCAL GRADIENT AND NOT THE AVERAGE" << std::endl; 
    
    // ATTN: Search for stable region for which the local layer gradient agrees well with associated cluster global gradient
    unsigned int gradientStabilityWindow(std::ceil(clusterMicroLayerFitResultMap.size() *  m_stableRegionClusterFraction));
    unsigned int goodLayerCount(0);
    
    for (int i = startLayer; i != loopTerminationLayer; i += step)
    {
        const auto microIter(clusterMicroLayerFitResultMap.find(i));

        if (microIter == clusterMicroLayerFitResultMap.end())
            continue;
        
        CartesianVector microDirection(0.f, 0.f, 0.f);
        clusterMicroFitResult.GetGlobalDirection(microIter->second.GetGradient(), microDirection);
    
        const float cosDirectionOpeningAngle(microDirection.GetCosOpeningAngle(associatedAverageDirection));
        if (cosDirectionOpeningAngle > m_mergePointMinCosAngleDeviation)
        {
            // ATTN: Set merge point and direction as that at the first layer in the stable region
            if (goodLayerCount == 0)
            {
                // ATTN: Cluster direction vectors must point to one another
                //clusterMergeDirection = clusterAverageDirection * (isUpstream ? 1.f : -1.f);
                //clusterMergeDirection = microDirection;
                clusterMergeDirection = clusterAverageDirection;
                clusterMicroFitResult.GetGlobalFitPosition(microIter->second.GetL(), clusterMergePosition);
            }
            
            ++goodLayerCount;
        }
        else
        {
            goodLayerCount = 0;
        }

        if (goodLayerCount > gradientStabilityWindow)
            break;
        
        // ATTN: Abort merging process have not found a stable region 
        if (i == endLayer)
            return false;                                         
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------    

} // namespace lar_content
