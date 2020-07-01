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
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

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
        const CartesianVector &averageMinLayerDirection(macroFitResult->second.GetGlobalMinLayerDirection()), &averageMaxLayerDirection(macroFitResult->second.GetGlobalMaxLayerDirection());

        std::cout << "MIN LAYER DIRECTION (AVERAGE): " << averageMinLayerDirection << std::endl;
        std::cout << "MAX LAYER DIRECTION (AVERAGE): " << averageMaxLayerDirection << std::endl;
        
        //const float upstreamX(microFitResult->second.GetGlobalMinLayerPosition().GetX()), downstreamX(microFitResult->second.GetGlobalMaxLayerPosition().GetX());

        //bool upstreamModification(true);
        /*
        if (upstreamX < (tpcEdge2 + bufferRegion) || upstreamX > (tpcEdge1 - bufferRegion))
        {
            upstreamModification = true;
            std::cout << "UPSTREAM POINT CLOSE ENOUGH TO TPC BOUNDARY" << std::endl;
        }
        */
        //bool downstreamModification(true);
        /*
        if (downstreamX > (tpcEdge1 - bufferRegion) || downstreamX < (tpcEdge2 + bufferRegion))
        {
            downstreamModification = true;
            std::cout << "DOWNSTREAM POINT CLOSE ENOUGH TO TPC BOUNDARY" << std::endl;
        }
        */

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


        this->NewIsCosmicRay(downstreamPosition, downstreamMergePoint, downstreamMergeDirection, pCluster, true);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        this->NewIsCosmicRay(upstreamPosition, upstreamMergePoint, upstreamMergeDirection, pCluster, false);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        
        /*
        this->IsCosmicRay(downstreamPosition, downstreamMergePoint, downstreamMergeDirection, true, microFitResult->second, averageMinLayerDirection, pCluster);
                PandoraMonitoringApi::ViewEvent(this->GetPandora());
                this->IsCosmicRay(upstreamPosition, upstreamMergePoint, upstreamMergeDirection, false, microFitResult->second, averageMinLayerDirection, pCluster);
                PandoraMonitoringApi::ViewEvent(this->GetPandora());
        */
        
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


    void CosmicRayEndpointCorrectionAlgorithm::NewIsCosmicRay(const CartesianVector &clusterEndpoint, const CartesianVector &clusterMergePoint, const CartesianVector &clusterMergeDirection, const Cluster *const pCluster, const bool isUpstream) const
{
    float straightLinePrediction(clusterEndpoint.GetX() - clusterMergePoint.GetX(), 0.f, clusterEndpoint.GetZ() - clusterMergePoint.GetZ());

    CartesianPointVector hitSubset;
    if (straightLinePrediction.GetMagnitude() > std::numeric_limits<float>::epsilon())
    {
        const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
        for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
        {
            for (const CaloHit *const pCaloHit : *mapEntry.second) 
            {
                const CartesianVector &position(pCaloHit->GetPositionVector());

                if ((isUpstream && (position.GetZ() > clusterMergePoint.GetZ())) || (!isUpstream && (position.GetZ() < clusterMergePoint.GetZ())))
                {
                    hitSubset.push_back(position);                
                }
            }
        }

        try
        {
            const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
            const int m_slidingFitWindow(5);
            const TwoDSlidingFitResult subsetFit(&hitSubset, m_slidingFitWindow, slidingFitPitch);

            const LayerFitResultMap &clusterMicroLayerFitResultMap(subsetFit.GetLayerFitResultMap());
            const int startLayer(subsetFit.GetMinLayer());
            const int endLayer(subsetFit.GetMaxLayer());
            const int loopTerminationLayer(endLayer + 1);
            const int step(1);

            int count(0);
            int angleCount(0);
            float angleSum(0.f);

            CartesianVector previousDirection(0.f,0.f,0.f);
            
            for (int i = startLayer; i != loopTerminationLayer; i += step)
            {
                const auto microIter(clusterMicroLayerFitResultMap.find(i));

                if (microIter == clusterMicroLayerFitResultMap.end())
                    continue;

                CartesianVector microDirection(0.f, 0.f, 0.f);
                subsetFit.GetGlobalDirection(microIter->second.GetGradient(), microDirection);
                //std::cout << "ANGLE DEVIATION: " << (microDirection.GetOpeningAngle(clusterMergeDirection)*180/3.14) << std::endl;

                if ((microDirection.GetOpeningAngle(clusterMergeDirection)*180/3.14) > 5)
                {
                    angleSum += (microDirection.GetOpeningAngle(clusterMergeDirection)*180/3.14);
                    ++angleCount;
                }
                
                ++count;
                
                if(count != 1)
                {
                    //std::cout << "ANGLE FROM LAST ANGLE: " << microDirection.GetOpeningAngle(previousDirection) << std::endl;
                }

                previousDirection = microDirection;

                CartesianVector position(0.f, 0.f, 0.f);
                subsetFit.GetGlobalFitPosition(microIter->second.GetL(), position);
                PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &position, "POSITION", BLACK, 2);
            }

            if (angleCount > 0)
                std::cout << "AVERAGE OPENING ANGLE: " << angleSum / angleCount << std::endl;
        }
        catch (StatusCodeException &)
        {
            std::cout << "CANNOT MAKE A FIT" << std::endl;
            return;
        }


    }

}




    
void CosmicRayEndpointCorrectionAlgorithm::IsCosmicRay(const CartesianVector &clusterEndpoint, const CartesianVector &clusterMergePoint, const CartesianVector &clusterMergeDirection, const bool isUpstream, const TwoDSlidingFitResult &microFitResult, const CartesianVector &averageDirection, const Cluster *const pCluster) const
{

    std::cout << "AVERAGE DIRECTION: " << averageDirection << std::endl;
    std::cout << "SEPARATION DISTANCE: " << std::sqrt(clusterMergePoint.GetDistanceSquared(clusterEndpoint)) << std::endl;
    
    CartesianVector straightLinePrediction(clusterEndpoint.GetX() - clusterMergePoint.GetX(), 0.f, clusterEndpoint.GetZ() - clusterMergePoint.GetZ());
    //float straightLinePredictionGradient(straightLinePrediction.GetZ() / straightLinePrediction.GetX());
    //float intercept(clusterMergePoint.GetZ() - (straightLinePredictionGradient * clusterMergePoint.GetX()));


    const CartesianVector  &endpointDirection(isUpstream ? microFitResult.GetGlobalMaxLayerPosition() : microFitResult.GetGlobalMinLayerPosition());

    if (straightLinePrediction.GetMagnitude() > std::numeric_limits<float>::epsilon())
    {

        //ORDERED BY MIN LAYER
        float hitSeparationDistance(0.f);
        CartesianVector firstHitPosition(0.f,0.f,0.f), previousHitPosition(0.f,0.f,0.f);
        int hitCount(0);
        const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
        for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
        {
            for (const CaloHit *const pCaloHit : *mapEntry.second) 
            {
                const CartesianVector &position(pCaloHit->GetPositionVector());

                if ((isUpstream && (position.GetZ() > clusterMergePoint.GetZ())) || (!isUpstream && (position.GetZ() < clusterMergePoint.GetZ())))
                {
                    ++hitCount;
                    
                    if (hitCount == 1)
                        firstHitPosition = position;
                    
                    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &position, "hit position", DARKGREEN, 2);

                    if (hitCount != 1)
                        hitSeparationDistance += (position - previousHitPosition).GetMagnitude();
                
                    previousHitPosition = position;
                
                }
            }
        }

        CartesianVector hitStraightLinePredition(isUpstream ? CartesianVector(previousHitPosition.GetX() - firstHitPosition.GetX(), 0.f, previousHitPosition.GetZ() - firstHitPosition.GetZ()) :
                                                                    CartesianVector(firstHitPosition.GetX() - previousHitPosition.GetX(), 0.f, firstHitPosition.GetZ() - previousHitPosition.GetZ()));

        hitStraightLinePredition = hitStraightLinePredition.GetUnitVector();

        const CartesianVector hitExtrapolationPoint(isUpstream ? firstHitPosition + (hitStraightLinePredition * hitSeparationDistance) : previousHitPosition + (hitStraightLinePredition * hitSeparationDistance));
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitExtrapolationPoint, "HIT EXTRAP POINT", VIOLET, 2);

        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        
        straightLinePrediction = straightLinePrediction.GetUnitVector();
        straightLinePrediction *= isUpstream ? 1.f : -1.f;
        std::cout << "DR OPENING ANGLE TO CR: " << clusterMergeDirection.GetOpeningAngle(straightLinePrediction)*180/3.14 << std::endl; //I THINK THIS SHOULD BE AN AVERAGE DIRECTION...

        const LayerFitResultMap &clusterMicroLayerFitResultMap(microFitResult.GetLayerFitResultMap());
        const int startLayer(isUpstream ? microFitResult.GetMaxLayer() : microFitResult.GetMinLayer());
        const int endLayer(isUpstream ? microFitResult.GetMinLayer() : microFitResult.GetMaxLayer());
        const int loopTerminationLayer(endLayer + (isUpstream ? -1 : 1));
        const int step(isUpstream ? -1 : 1);

        CartesianVector previousDirection(0.f,0.f,0.f), previousPosition(0.f,0.f,0.f);

        int count(0);
        float openingAngleSum(0.f);
        float separationDistance(0.0);
        float averageDistanceFromLine(0.f);
        float straightLineOpeningAngleSum(0.f);

        float splitL(0.f), splitT(0.f);
        microFitResult.GetLocalPosition(clusterMergePoint, splitL, splitT);

        for (int i = startLayer; i != loopTerminationLayer; i += step)
        {
            const auto microIter(clusterMicroLayerFitResultMap.find(i));

            if (microIter == clusterMicroLayerFitResultMap.end())
                continue;

            if ((isUpstream ? microIter->second.GetL() > splitL : microIter->second.GetL() < splitL) && count < 100)
            {

                    ++count;
                
                    CartesianVector microDirection(0.f, 0.f, 0.f);
                    microFitResult.GetGlobalDirection(microIter->second.GetGradient(), microDirection);
                    
                    //float openingAngle(microDirection.GetOpeningAngle(downstreamMergeDirection)*180/3.14);
                    //std::cout << "OPENING ANGLE: " << openingAngle << std::endl;
                    //if (openingAngle < 10.f)
                    //continue;
                    
                    openingAngleSum += (microDirection.GetOpeningAngle(clusterMergeDirection)*180/3.14);
                    straightLineOpeningAngleSum += (straightLinePrediction.GetOpeningAngle(microDirection)*180/3.14);
                    

                    // LOOK AT THE ANGLE FOR THE STRAIGHT LINE FIT???????? (might tell us how curvy it is)

                    CartesianVector position(0.f, 0.f, 0.f);
                    microFitResult.GetGlobalFitPosition(microIter->second.GetL(), position);
                    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &position, "POSITION", BLUE, 2);

                    //bool isAbove(position.GetZ() > ((straightLinePredictionGradient * position.GetZ()) + intercept));
                    //std::cout << (isAbove ? "ABOVE" : "BELOW") << std::endl;

                    averageDistanceFromLine += straightLinePrediction.GetCrossProduct(position - clusterMergePoint).GetMagnitude();

                    if(count != 1)
                    {
                        std::cout << "OPENING ANGLE FROM LAST LAYER: " << microDirection.GetOpeningAngle(previousDirection)*180/3.14 << std::endl;
                        std::cout << "separationDistance" << separationDistance << std::endl;
                        separationDistance += (position - previousPosition).GetMagnitude();
                        std::cout << "separationDistance" << separationDistance << std::endl;
                        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &position, "POSITION", VIOLET, 2);
                        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &previousPosition, "PREVIOUS POSITION", BLACK, 2);
                        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
                    }

                    previousDirection = microDirection;
                    previousPosition = position;
            }
        }

            separationDistance += std::sqrt(previousPosition.GetDistanceSquared(clusterMergePoint));


                /*
                if (microIter->second.GetL() < splitL)
                {
                    CartesianVector microDirection(0.f, 0.f, 0.f);
                    microFitResult.GetGlobalDirection(microIter->second.GetGradient(), microDirection);                    
                    CartesianVector position(0.f, 0.f, 0.f);
                    microFitResult.GetGlobalFitPosition(microIter->second.GetL(), position);
                    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &position, "POSITION", RED, 2);
                    std::cout << "AFTER THE SPLIT POINT" << std::endl;
                    std::cout << "OPENING ANGLE FROM LAST LAYER: " << microDirection.GetOpeningAngle(previousDirection)*180/3.14 << std::endl;
                    previousDirection = microDirection;

                }
                */
                

            if (count > 0)
            {
                float impactL(0.f), impactT(0.f);
                CartesianVector extrapolatedPoint(clusterMergePoint + (isUpstream ? (straightLinePrediction * separationDistance) : (straightLinePrediction * separationDistance * (-1.f))));
                PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedPoint, "EXTRAPOLATED MERGE POINT", GREEN, 2);
            
                LArPointingClusterHelper::GetImpactParameters(extrapolatedPoint, straightLinePrediction, clusterEndpoint, impactL, impactT);
                std::cout << "impactL: " << impactL << std::endl;
                std::cout << "impactT: " << impactT << std::endl;

                std::cout << "OVERAL ANGLE CHANGE: " << endpointDirection.GetOpeningAngle(clusterMergeDirection)*180/3.14 << std::endl;
                std::cout << "averageDistanceFromLine" << averageDistanceFromLine / static_cast<float>(count) << std::endl;         
                std::cout << "OPENING ANGLE AVERAGE: " << openingAngleSum / static_cast<float>(count) << std::endl;
                std::cout << "straightLineOpeningAngleSum" << straightLineOpeningAngleSum / static_cast<float>(count) << std::endl;
            }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------    

StatusCode CosmicRayEndpointCorrectionAlgorithm::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    std::cout << "frog" << std::endl;
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
