/**
 *  @file   larpandoracontent/LArTwoDReco/LArCosmicRay/CosmicRayEndpointCorrectionAlgorithm.cc 
 *
 *  @brief  Implementation of the cosmic ray endpoint correction class
 *
 *  $Log: $
 */
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArTwoDReco/LArCosmicRay/CosmicRayEndpointCorrectionAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

using namespace pandora;

namespace lar_content
{
    
CosmicRayEndpointCorrectionAlgorithm::CosmicRayEndpointCorrectionAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------ 
    
StatusCode CosmicRayEndpointCorrectionAlgorithm::Run()
{
    const ClusterList *pClusterList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));
    
    ClusterVector clusterVector;
    this->SelectCleanClusters(pClusterList, clusterVector);

    TwoDSlidingFitResultMap microSlidingFitResultMap, macroSlidingFitResultMap;
    this->InitialiseSlidingFitResultMaps(clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap);

    const LArTPC &pLArTPC(this->GetPandora().GetGeometry()->GetLArTPC());

    const float tpcEdge1(pLArTPC.GetCenterX() + (pLArTPC.GetWidthX() / 2.f));
    const float tpcEdge2(pLArTPC.GetCenterX() - (pLArTPC.GetWidthX() / 2.f));

    const float minZ(pLArTPC.GetCenterZ() - (pLArTPC.GetWidthZ() / 2.f));
    const float maxZ(pLArTPC.GetCenterZ() + (pLArTPC.GetWidthZ() / 2.f));

    const float bufferRegion(10.f);

    CartesianVector tpcBoundary1_1(tpcEdge1 - bufferRegion, 0.f, minZ), tpcBoundary1_2(tpcEdge1 - bufferRegion, 0.f, maxZ);
    CartesianVector tpcBoundary2_1(tpcEdge2 + bufferRegion, 0.f, minZ), tpcBoundary2_2(tpcEdge2 + bufferRegion, 0.f, maxZ);


    for (const Cluster *const pCluster : clusterVector)
    {

    PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &tpcBoundary1_1, &tpcBoundary1_2, "EDGE", RED, 2, 2);
    PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &tpcBoundary2_1, &tpcBoundary2_2, "EDGE", RED, 2, 2);
        
        TwoDSlidingFitResultMap::const_iterator microFitResult(microSlidingFitResultMap.find(pCluster));
        TwoDSlidingFitResultMap::const_iterator macroFitResult(macroSlidingFitResultMap.find(pCluster));

        if (microFitResult == microSlidingFitResultMap.end() || macroFitResult == macroSlidingFitResultMap.end())
            continue;        

        const CartesianVector &upstreamPosition(microFitResult->second.GetGlobalMinLayerPosition()), &downstreamPosition(microFitResult->second.GetGlobalMaxLayerPosition());
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &upstreamPosition, "UPSTREAM", RED, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &downstreamPosition, "DOWNSTREAM", RED, 2);
        
        const float upstreamX(microFitResult->second.GetGlobalMinLayerPosition().GetX()), downstreamX(microFitResult->second.GetGlobalMaxLayerPosition().GetX());

        bool upstreamModification(false);
        if (upstreamX < (tpcEdge2 + bufferRegion) || upstreamX > (tpcEdge1 - bufferRegion))
        {
            upstreamModification = true;
            std::cout << "UPSTREAM POINT CLOSE ENOUGH TO TPC BOUNDARY" << std::endl;
        }
        
        bool downstreamModification(false);
        if (downstreamX > (tpcEdge1 - bufferRegion) || downstreamX < (tpcEdge2 + bufferRegion))
        {
            downstreamModification = true;
            std::cout << "DOWNSTREAM POINT CLOSE ENOUGH TO TPC BOUNDARY" << std::endl;
        }

        CartesianVector upstreamMergePoint(0.f, 0.f, 0.f), upstreamMergeDirection(0.f, 0.f, 0.f), downstreamMergePoint(0.f, 0.f, 0.f), downstreamMergeDirection(0.f, 0.f, 0.f);
        if (!GetClusterMergingCoordinates(microFitResult->second, macroFitResult->second, macroFitResult->second, false, upstreamMergePoint, upstreamMergeDirection))
            continue;
        if (!GetClusterMergingCoordinates(microFitResult->second, macroFitResult->second, macroFitResult->second, true, downstreamMergePoint, downstreamMergeDirection))
            continue;

        ClusterList theCluster({pCluster});
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &theCluster, "THE CLUSTER", BLUE);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &upstreamMergePoint, "UPSTREAM MERGE POINT", BLACK, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &downstreamMergePoint, "DOWNSTREAM MERGE POINT", BLACK, 2);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());

        
        if (downstreamModification)
        {
            std::cout << "SEPARATION DISTANCE: " << downstreamMergePoint.GetDistanceSquared(downstreamPosition) << std::endl;
            CartesianVector straightLinePrediction(downstreamPosition.GetX() - downstreamMergePoint.GetX(), 0.f, downstreamPosition.GetZ() - downstreamMergePoint.GetZ());
            
            if (straightLinePrediction.GetMagnitude() > std::numeric_limits<float>::epsilon())
            {
            straightLinePrediction = straightLinePrediction.GetUnitVector();

            std::cout << "DR OPENING ANGLE TO CR: " << downstreamMergeDirection.GetOpeningAngle(straightLinePrediction)*180/3.14 << std::endl;

            bool isUpstream(true);
            const LayerFitResultMap &clusterMicroLayerFitResultMap(microFitResult->second.GetLayerFitResultMap());
            const int startLayer(isUpstream ? microFitResult->second.GetMaxLayer() : microFitResult->second.GetMinLayer());
            const int endLayer(isUpstream ? microFitResult->second.GetMinLayer() : microFitResult->second.GetMaxLayer());
            const int loopTerminationLayer(endLayer + (isUpstream ? -1 : 1));
            const int step(isUpstream ? -1 : 1);

            int count(0);
            float openingAngleSum(0.f);
            for (int i = startLayer; i != loopTerminationLayer; i += step)
            {
                const auto microIter(clusterMicroLayerFitResultMap.find(i));

                if (microIter == clusterMicroLayerFitResultMap.end())
                    continue;

                float splitL(0.f), splitT(0.f);
                microFitResult->second.GetLocalPosition(downstreamMergePoint, splitL, splitT);
                
                if (microIter->second.GetL() > splitL && count < 3)
                {
                    ++count;
                    CartesianVector microDirection(0.f, 0.f, 0.f);
                    microFitResult->second.GetGlobalDirection(microIter->second.GetGradient(), microDirection);
                    openingAngleSum += (microDirection.GetOpeningAngle(downstreamMergeDirection)*180/3.14);

                    CartesianVector position(0.f, 0.f, 0.f);
                    microFitResult->second.GetGlobalFitPosition(microIter->second.GetL(), position);
                    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &position, "POSITION", BLUE, 2);
                }
            }

            std::cout << "OPENING ANGLE AVERAGE: " << openingAngleSum / static_cast<float>(count) << std::endl;
            }
            PandoraMonitoringApi::ViewEvent(this->GetPandora());
        }

        if (upstreamModification)
        {
            std::cout << "SEPARATION DISTANCE: " << upstreamMergePoint.GetDistanceSquared(upstreamPosition)*180/3.14 << std::endl;
            CartesianVector straightLinePrediction(upstreamPosition.GetX() - upstreamMergePoint.GetX(), 0.f, upstreamPosition.GetZ() - upstreamMergePoint.GetZ());
            
            if (straightLinePrediction.GetMagnitude() > std::numeric_limits<float>::epsilon())
            {
            
            straightLinePrediction = straightLinePrediction.GetUnitVector();
            straightLinePrediction *= (-1.f);

            std::cout << "DR OPENING ANGLE TO CR: " << upstreamMergeDirection.GetOpeningAngle(straightLinePrediction) << std::endl;
            
            bool isUpstream(false);
            const LayerFitResultMap &clusterMicroLayerFitResultMap(microFitResult->second.GetLayerFitResultMap());
            const int startLayer(isUpstream ? microFitResult->second.GetMaxLayer() : microFitResult->second.GetMinLayer());
            const int endLayer(isUpstream ? microFitResult->second.GetMinLayer() : microFitResult->second.GetMaxLayer());
            const int loopTerminationLayer(endLayer + (isUpstream ? -1 : 1));
            const int step(isUpstream ? -1 : 1);

            int count(0);
            float openingAngleSum(0.f);            
            for (int i = startLayer; i != loopTerminationLayer; i += step)
            {
                const auto microIter(clusterMicroLayerFitResultMap.find(i));

                if (microIter == clusterMicroLayerFitResultMap.end())
                    continue;

                float splitL(0.f), splitT(0.f);
                microFitResult->second.GetLocalPosition(upstreamMergePoint, splitL, splitT);
                
                if (microIter->second.GetL() < splitL && count < 5)
                {
                    ++count;
                    CartesianVector microDirection(0.f, 0.f, 0.f);
                    microFitResult->second.GetGlobalDirection(microIter->second.GetGradient(), microDirection);
                    openingAngleSum += (microDirection.GetOpeningAngle(upstreamMergeDirection)*180/3.14);

                    CartesianVector position(0.f, 0.f, 0.f);
                    microFitResult->second.GetGlobalFitPosition(microIter->second.GetL(), position);
                    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &position, "POSITION", BLUE, 2);
                }
            }

            std::cout << "OPENING ANGLE AVERAGE: " << openingAngleSum / static_cast<float>(count) << std::endl;            
            }
            PandoraMonitoringApi::ViewEvent(this->GetPandora());
        }
        
        
        // search for region along theorised CR that has the largest angle deviation from the mergin point direction
        // if found, remove hits??
        
    }
    
    return STATUS_CODE_SUCCESS;
}    

//------------------------------------------------------------------------------------------------------------------------------------------    

void CosmicRayEndpointCorrectionAlgorithm::SelectCleanClusters(const ClusterList *pClusterList, ClusterVector &clusterVector) const
{
    for (const Cluster *const pCluster : *pClusterList)
    {
        if (pCluster->GetNCaloHits() < 100)
            continue;

        clusterVector.push_back(pCluster);
    }
    
    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------    
    


//------------------------------------------------------------------------------------------------------------------------------------------    

StatusCode CosmicRayEndpointCorrectionAlgorithm::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    std::cout << "frog" << std::endl;
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
