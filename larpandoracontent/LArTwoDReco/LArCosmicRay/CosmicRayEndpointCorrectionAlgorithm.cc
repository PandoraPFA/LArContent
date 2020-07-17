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

    for (const Cluster *const pCluster : clusterVector)
    {
        TwoDSlidingFitResultMap::const_iterator microFitResult(microSlidingFitResultMap.find(pCluster));
        TwoDSlidingFitResultMap::const_iterator macroFitResult(macroSlidingFitResultMap.find(pCluster));

        if (microFitResult == microSlidingFitResultMap.end() || macroFitResult == macroSlidingFitResultMap.end())
            continue;        

        const CartesianVector &upstreamPosition(microFitResult->second.GetGlobalMinLayerPosition()), &downstreamPosition(microFitResult->second.GetGlobalMaxLayerPosition());
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &upstreamPosition, "UPSTREAM", RED, 2);
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &downstreamPosition, "DOWNSTREAM", RED, 2);

        CartesianVector upstreamMergePoint(0.f, 0.f, 0.f), upstreamMergeDirection(0.f, 0.f, 0.f), downstreamMergePoint(0.f, 0.f, 0.f), downstreamMergeDirection(0.f, 0.f, 0.f);
        if (!GetClusterMergingCoordinates(microFitResult->second, macroFitResult->second, macroFitResult->second, false, upstreamMergePoint, upstreamMergeDirection))
            continue;
        if (!GetClusterMergingCoordinates(microFitResult->second, macroFitResult->second, macroFitResult->second, true, downstreamMergePoint, downstreamMergeDirection))
            continue;

        /*
        ClusterList theCluster({pCluster});
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &theCluster, "THE CLUSTER", BLUE);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &upstreamMergePoint, "UPSTREAM MERGE POINT", BLACK, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &downstreamMergePoint, "DOWNSTREAM MERGE POINT", BLACK, 2);
        */

        if (this->NewIsCosmicRay(downstreamPosition, downstreamMergePoint, downstreamMergeDirection, pCluster, true))
        {
            std::cout << "WILL REMOVE" << std::endl;
            PandoraMonitoringApi::ViewEvent(this->GetPandora());
        }
        else
        {
            std::cout << "WILL NOT REMOVE" << std::endl;
        }
        
        if (this->NewIsCosmicRay(upstreamPosition, upstreamMergePoint, upstreamMergeDirection, pCluster, false))
        {
            std::cout << "WILL REMOVE" << std::endl;
            PandoraMonitoringApi::ViewEvent(this->GetPandora());
        }
        else
        {
            std::cout << "WILL NOT REMOVE" << std::endl;
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
        if (pCluster->GetNCaloHits() < 50)
            continue;

        clusterVector.push_back(pCluster);
    }
    
    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------


bool CosmicRayEndpointCorrectionAlgorithm::NewIsCosmicRay(const CartesianVector &clusterEndpoint, const CartesianVector &clusterMergePoint,
    const CartesianVector &clusterMergeDirection, const Cluster *const pCluster, const bool isUpstream) const
{

    ClusterList theCluster({pCluster});
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &theCluster, "THE CLUSTER", BLUE);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &clusterEndpoint, "ENDPOINT", BLACK, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &clusterMergePoint, "CLUSTER MERGE POINT", RED, 2);
    
    // Whether the CR endpoint is near the TPC boundary
    const LArTPC &pLArTPC(this->GetPandora().GetGeometry()->GetLArTPC());
    const float tpcHighXEdge(pLArTPC.GetCenterX() + (pLArTPC.GetWidthX() / 2.f)), tpcLowXEdge(pLArTPC.GetCenterX() - (pLArTPC.GetWidthX() / 2.f));
    const float allowanceWidth(15.f), allowanceHighXEdge(tpcHighXEdge - allowanceWidth), allowanceLowXEdge(tpcLowXEdge + allowanceWidth);

    const bool isClusterEndpointInBoundary((clusterEndpoint.GetX() < allowanceLowXEdge) || (clusterEndpoint.GetX() > allowanceHighXEdge));
    const bool isClusterMergePointInBoundary((clusterMergePoint.GetX() < allowanceLowXEdge) || (clusterMergePoint.GetX() > allowanceHighXEdge));

    std::cout << "Distance from boundary: " << std::min(std::fabs(clusterMergePoint.GetX() - tpcHighXEdge), std::fabs(clusterMergePoint.GetX() - tpcLowXEdge)) << std::endl;

    //ISOBEL: NEED TO MAKE SURE NOT THE OUTER EDGE OF THE DETECTOR
    
    if (!(isClusterEndpointInBoundary || isClusterMergePointInBoundary))
    {
        std::cout << "NOT CLOSE ENOUGH TO THE BOUNDARY" << std::endl;
        PandoraMonitoringApi::ViewEvent(this->GetPandora());        
        return false;
    }

    // Investigate the separation of the endpoints points
    const float endpointSeparation((clusterEndpoint - clusterMergePoint).GetMagnitude());
    std::cout << "Endpoint Separation: " << endpointSeparation << std::endl;

    // ISOBEL: MAKE THIS LARGER?
    if (endpointSeparation < std::numeric_limits<float>::epsilon())
    {
        std::cout << "MERGE POINT AND ENDPOINT ARE THE SAME" << std::endl;
        PandoraMonitoringApi::ViewEvent(this->GetPandora());        
        return false;
    }

    // Investigate the cluster subset
    CartesianPointVector hitSubset;
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
    for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
    {
        for (const CaloHit *const pCaloHit : *mapEntry.second) 
        {
            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());

            if ((isUpstream && (hitPosition.GetZ() > clusterMergePoint.GetZ())) || (!isUpstream && (hitPosition.GetZ() < clusterMergePoint.GetZ())))
                hitSubset.push_back(hitPosition);          
        }
    }

    TwoDSlidingFitResultList subsetFitVector;
    try
    {
        const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const int m_slidingFitWindow(5);
        const TwoDSlidingFitResult subsetFit(&hitSubset, m_slidingFitWindow, slidingFitPitch);

        subsetFitVector.push_back(subsetFit);
    }
    catch (const StatusCodeException &)
    {
        std::cout << "CANNOT MAKE A FIT" << std::endl;
        PandoraMonitoringApi::ViewEvent(this->GetPandora());        
        return false;
    }

    // INVESTIGATE THE STRAIGHTNESS OF THE CLUSTER SUBSET
    const CartesianVector &minPosition(subsetFitVector.front().GetGlobalMinLayerPosition()), &maxPosition(subsetFitVector.front().GetGlobalMaxLayerPosition());
    CartesianVector straightLinePrediction(maxPosition - minPosition);
    straightLinePrediction = straightLinePrediction.GetUnitVector();
    
    float distanceFromLine(0.f);
    unsigned int distanceCount(0);
    for (const CartesianVector &hitPosition : hitSubset)
    {
        if (hitPosition.GetZ() < minPosition.GetZ() || hitPosition.GetZ() > maxPosition.GetZ())
            continue;

        distanceFromLine += (straightLinePrediction.GetCrossProduct(hitPosition - minPosition).GetMagnitude());

        ++distanceCount;
    }

    if (distanceCount == 0)
    {
        std::cout << "DISTANCE COUNT TOO SMALL" << std::endl;
        PandoraMonitoringApi::ViewEvent(this->GetPandora());        
        return false;
    }

    const float averageDistanceFromLine(distanceFromLine / distanceCount);
    std::cout << "averageDistanceFromLine: " << averageDistanceFromLine << std::endl;    
    
    if (averageDistanceFromLine < 0.3)
    {
        std::cout << "TOO STRAIGHT" << std::endl;
        PandoraMonitoringApi::ViewEvent(this->GetPandora());        
        return false;
    }

    
    // INVESTIGATE ENDPOINT Z SEPARATION
    const float predictedGradient(clusterMergeDirection.GetZ() / clusterMergeDirection.GetX());
    const float predictedIntercept(clusterMergePoint.GetZ() - (predictedGradient * clusterMergePoint.GetX()));
    
    const CartesianVector extrapolatedPoint(clusterEndpoint.GetX(), 0.f, predictedIntercept + (predictedGradient * clusterEndpoint.GetX()));
    const float deltaZ(std::fabs(clusterEndpoint.GetZ() - (predictedIntercept + (predictedGradient * clusterEndpoint.GetX()))));
    const float extrapolatedSeparation((extrapolatedPoint - clusterMergePoint).GetMagnitude());

    const CartesianVector start(clusterMergePoint + (clusterMergeDirection*20));
    const CartesianVector end(clusterMergePoint - (clusterMergeDirection*20));
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedPoint, "JAM", GREEN, 2);
    PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &start, &end, "JAM LINE", GREEN, 2, 2);
    
    const float scaledDeltaZ(deltaZ / extrapolatedSeparation); 
    std::cout << "deltaZ: " << deltaZ << std::endl;
    std::cout << "scaled deltaZ: " << scaledDeltaZ << std::endl;

    if (scaledDeltaZ < 0.25)
    {
        std::cout << "SCALED DELTA Z IS TOO HIGH" << std::endl;
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return false;
    }

    // INVESTIGATE THE DIRECTION CHANGE
    const TwoDSlidingFitResult subsetFit(subsetFitVector.front());
    const LayerFitResultMap &clusterMicroLayerFitResultMap(subsetFit.GetLayerFitResultMap());
    const int startLayer(isUpstream ? subsetFit.GetMinLayer() : subsetFit.GetMaxLayer());
    const int endLayer(isUpstream ? subsetFit.GetMaxLayer() : subsetFit.GetMinLayer());
    const int loopTerminationLayer(isUpstream ? endLayer + 1 : endLayer - 1);
    const int step(isUpstream ? 1 : -1);

    bool reachedFirstCurve(false);
    bool isClockwise(false);

    float previousOpeningAngle(0.f);

    unsigned int numberOfWobbles(0);
    unsigned int maxNumberOfWobbles(3);

    bool metCriteria(false);
    for (int i = startLayer; i != loopTerminationLayer; i += step)
    {
        const auto microIter(clusterMicroLayerFitResultMap.find(i));

        if (microIter == clusterMicroLayerFitResultMap.end())
            continue;

        CartesianVector microDirection(0.f, 0.f, 0.f);
        subsetFit.GetGlobalDirection(microIter->second.GetGradient(), microDirection);
        float microOpeningAngle(microDirection.GetOpeningAngle(clusterMergeDirection) * 180 / 3.14);

        /////////////////////
        CartesianVector microPosition(0.f, 0.f, 0.f);
        subsetFit.GetGlobalPosition(microIter->second.GetL(), microIter->second.GetFitT(), microPosition);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &microPosition, "MICRO POSITION", BLACK, 2);
        /////////////////////
        
        bool isAbove(microDirection.GetZ() > (predictedGradient*microDirection.GetX()));
        if (!isAbove)
            microOpeningAngle *= (-1.f);

        std::cout << "ANGLE DEVIATION: " << microOpeningAngle << std::endl;

        float layerAngleDeviation(0.f);
        if (reachedFirstCurve)
        {
            layerAngleDeviation = (std::fabs(microOpeningAngle - previousOpeningAngle));
            std::cout << "<-------------------- CHANGE FROM LAST LAYER: " << layerAngleDeviation << std::endl;
        }
        
        if (((microOpeningAngle > 10.f) || (microOpeningAngle < -10.f)) && !reachedFirstCurve)
        {
            reachedFirstCurve = true;
            isClockwise = (microOpeningAngle < 0.f);
            previousOpeningAngle = microOpeningAngle;
            continue;
        }

        if (reachedFirstCurve)
        {
            if (((isClockwise && (microOpeningAngle > previousOpeningAngle)) || (!isClockwise && (microOpeningAngle < previousOpeningAngle))) && (layerAngleDeviation > 1 ))
            {
                ++numberOfWobbles;

                if (numberOfWobbles > maxNumberOfWobbles)
                {
                    std::cout << "TOO MANY WOBBLES - RESET" << std::endl;
                    PandoraMonitoringApi::ViewEvent(this->GetPandora());
                    return false;
                }
            }

            if ((isClockwise && (microOpeningAngle < -25.f)) || (!isClockwise && (microOpeningAngle > 25.f)))
            {
                std::cout << "MEET ANGLE CRITERIA" << std::endl;
                metCriteria = true;
                break;
            }
        }

        previousOpeningAngle = microOpeningAngle;
    }


    if(!metCriteria)
    {
        std::cout << "DIDN'T FIND ANGLE CRITERIA" << std::endl;
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return false;
    }
    
    return true;

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
