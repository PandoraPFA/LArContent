/**
 *  @file   larpandoracontent/LArTwoDReco/LArCosmicRay/CosmicRayEndpointCorrectionAlgorithm.cc 
 *
 *  @brief  Implementation of the cosmic ray endpoint correction class
 *
 *  $Log: $
 */
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArControlFlow/MultiPandoraApi.h"

#include "larpandoracontent/LArTwoDReco/LArCosmicRay/CosmicRayEndpointCorrectionAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArStitchingHelper.h"

using namespace pandora;

namespace lar_content
{
    
CosmicRayEndpointCorrectionAlgorithm::CosmicRayEndpointCorrectionAlgorithm() :
    m_minCaloHits(50),
    m_maxDistanceFromTPC(15.f),
    m_curveThreshold(0.3f),
    m_minScaledZOffset(2.f), //0.25
    m_thresholdAngleDeviation(10.f),
    m_thresholdAngleDeviationBetweenLayers(1.f),
    m_maxAnomalousPoints(2),
    m_thresholdMaxAngleDeviation(21.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------    
    
void CosmicRayEndpointCorrectionAlgorithm::FindBestClusterAssociation(const ClusterVector &clusterVector, const SlidingFitResultMapPair &slidingFitResultMapPair,
    ClusterAssociationVector &clusterAssociationVector)
{
    // Collect LArTPC information
    const Pandora *primaryPandoraInstance;
    try
    {
        primaryPandoraInstance = MultiPandoraApi::GetPrimaryPandoraInstance(&this->GetPandora());
    }
    catch (const StatusCodeException &)
    {
        primaryPandoraInstance = &this->GetPandora();
    }
            
    float detectorMinXEdge(std::numeric_limits<float>::max());
    float detectorMaxXEdge(-std::numeric_limits<float>::max());

    const LArTPC *pLArTPC(&this->GetPandora().GetGeometry()->GetLArTPC());
    const LArTPCMap &larTPCMap(primaryPandoraInstance->GetGeometry()->GetLArTPCMap());
    
    for (const LArTPCMap::value_type &mapEntry : larTPCMap)
    {
        const LArTPC *const pSubLArTPC(mapEntry.second);
        detectorMinXEdge = std::min(detectorMinXEdge, pSubLArTPC->GetCenterX() - 0.5f * pSubLArTPC->GetWidthX());
        detectorMaxXEdge = std::max(detectorMaxXEdge, pSubLArTPC->GetCenterX() + 0.5f * pSubLArTPC->GetWidthX());

        if (std::fabs(pSubLArTPC->GetCenterX() - pLArTPC->GetCenterX()) < std::numeric_limits<float>::epsilon())
            pLArTPC = pSubLArTPC;
    }
    
    //ATTN: Method assumes clusters ordered by their hits 
    for (const Cluster *const pCluster : clusterVector)
    {
        const TwoDSlidingFitResult &microSlidingFitResult(slidingFitResultMapPair.first->at(pCluster));
        const TwoDSlidingFitResult &macroSlidingFitResult(slidingFitResultMapPair.second->at(pCluster));
        
        for (unsigned int isEndUpstream = 0; isEndUpstream < 2; ++isEndUpstream)
        {
            /////////////////
            ClusterList theCluster({pCluster});
            PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &theCluster, "CONSIDERED CLUSTER", BLACK);
            ////////////////
            
            std::cout << "isEndUpstream: " << isEndUpstream << std::endl;
            
            const bool isClusterUpstream(!isEndUpstream);
            
            CartesianVector clusterMergePoint(0.f, 0.f, 0.f), clusterMergeDirection(0.f, 0.f, 0.f);
            if (!GetClusterMergingCoordinates(microSlidingFitResult, macroSlidingFitResult, macroSlidingFitResult, isClusterUpstream, clusterMergePoint, clusterMergeDirection))
            {
                std::cout << "CANNOT FIND MERGE POSITION" << std::endl;
                PandoraMonitoringApi::ViewEvent(this->GetPandora());
                continue;
            }
            
            const CartesianVector &endpointPosition(isEndUpstream ? microSlidingFitResult.GetGlobalMinLayerPosition() : microSlidingFitResult.GetGlobalMaxLayerPosition());
            const float endpointSeparation((endpointPosition - clusterMergePoint).GetMagnitude());
            
            /////////////////
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &endpointPosition, "ENDPOINT", BLUE, 2);
            std::cout << "Endpoint Separation: " << endpointSeparation << std::endl;
            /////////////////    

            if (endpointSeparation < std::numeric_limits<float>::epsilon())
            {
                std::cout << "MERGE POINT AND ENDPOINT ARE THE SAME" << std::endl;
                PandoraMonitoringApi::ViewEvent(this->GetPandora());
                continue;
            }

            // INVESTIGATE ENDPOINT Z SEPARATION
            const float predictedGradient(clusterMergeDirection.GetZ() / clusterMergeDirection.GetX()), predictedIntercept(clusterMergePoint.GetZ() - (predictedGradient * clusterMergePoint.GetX()));
            const CartesianVector fitEndpointPosition(endpointPosition.GetX(), 0.f, predictedIntercept + (predictedGradient * endpointPosition.GetX()));

            ///////////////////
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &fitEndpointPosition, "FIT ENDPOINT", ORANGE, 2);
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &clusterMergePoint, "clusterMergePoint", VIOLET, 2);
            ///////////////////
            
            const float deltaZ(std::fabs(endpointPosition.GetZ() - fitEndpointPosition.GetZ()));
            //const float extrapolatedSeparation((fitEndpointPosition - clusterMergePoint).GetMagnitude());
            const float scaledDeltaZ(deltaZ / endpointSeparation);
    
            /////////////////
            const CartesianVector start(clusterMergePoint + (clusterMergeDirection*20));
            const CartesianVector end(clusterMergePoint - (clusterMergeDirection*20));
            PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &start, &end, "JAM LINE", GREEN, 4, 2);
            std::cout << "deltaZ: " << deltaZ << std::endl;
            std::cout << "scaled deltaZ: " << scaledDeltaZ << std::endl;    
            /////////////////

            if (deltaZ < m_minScaledZOffset)
            {
                std::cout << "SCALED DELTA Z IS TOO LOW" << std::endl;
                PandoraMonitoringApi::ViewEvent(this->GetPandora());
                continue;
            }
    
            // INVESTIGATE PROXIMITY TO TPC BOUNDARY
           
            
            bool isNearestBoundaryHigherX((endpointPosition.GetX() - pLArTPC->GetCenterX()) > 0.f);
            const float nearestTPCBoundaryX(pLArTPC->GetCenterX() + (pLArTPC->GetWidthX() * 0.5f * (isNearestBoundaryHigherX ? 1.f : -1.f)));

            if ((std::fabs(nearestTPCBoundaryX - detectorMinXEdge) < std::numeric_limits<float>::epsilon()) ||
                (std::fabs(nearestTPCBoundaryX - detectorMaxXEdge) < std::numeric_limits<float>::epsilon()))
            {
                std::cout << "ENDPOINT IS NEAR DETECTOR EDGE" << std::endl;
                PandoraMonitoringApi::ViewEvent(this->GetPandora());
                continue;
            }

            const float allowanceBoundaryX(nearestTPCBoundaryX + (m_maxDistanceFromTPC * (isNearestBoundaryHigherX ? -1.f : 1.f)));
            if (isNearestBoundaryHigherX && !((clusterMergePoint.GetX() > allowanceBoundaryX) || (endpointPosition.GetX() > allowanceBoundaryX)))
            {
                std::cout << "NOT CLOSE ENOUGH TO THE BOUNDARY" << std::endl;
                PandoraMonitoringApi::ViewEvent(this->GetPandora());
                continue;
            }

            if (!isNearestBoundaryHigherX && !((clusterMergePoint.GetX() < allowanceBoundaryX) || (endpointPosition.GetX() < allowanceBoundaryX)))
            {
                std::cout << "NOT CLOSE ENOUGH TO THE BOUNDARY" << std::endl;
                PandoraMonitoringApi::ViewEvent(this->GetPandora());
                continue;
            }            
                
            if(this->IsDeltaRay(pCluster, clusterMergePoint, clusterMergeDirection, isEndUpstream, clusterMergePoint))
            {
                const LArTPC &neighboughTPC(LArStitchingHelper::FindClosestTPC(*primaryPandoraInstance, *pLArTPC, isNearestBoundaryHigherX));

                std::cout << "ADDRESS neighbour: " << &neighboughTPC << std::endl;
                std::cout << "ADDRESS: " << &pLArTPC << std::endl;
                
                const float gapSizeX(std::fabs(neighboughTPC.GetCenterX() - pLArTPC->GetCenterX()) - (neighboughTPC.GetWidthX() * 0.5f) - (pLArTPC->GetWidthX() * 0.5f));

                std::cout << "THIS TPC CENTRE: " << pLArTPC->GetCenterX() << std::endl;
                std::cout << "NEIGHBOUR TPC CENTRE: " << neighboughTPC.GetCenterX() << std::endl;
                
                const float endpointX(nearestTPCBoundaryX + (gapSizeX * 0.5f * (isNearestBoundaryHigherX ? 1.f : -1.f)));
                const CartesianVector extrapolatedEndpointPosition(endpointX, 0.f, predictedIntercept + (predictedGradient * endpointX));
                //const CartesianVector extrapolatedEndpointPosition(clusterMergePoint + (clusterMergeDirection * endpointSeparation * (isEndUpstream ? -1.f : 1.f)));
                
                ClusterEndpointAssociation clusterEndpointAssociation(isEndUpstream ?
                    ClusterEndpointAssociation(extrapolatedEndpointPosition, clusterMergeDirection, clusterMergePoint, clusterMergeDirection * (-1.f), pCluster, true) :
                    ClusterEndpointAssociation(clusterMergePoint, clusterMergeDirection, extrapolatedEndpointPosition, clusterMergeDirection * (-1.f), pCluster, false));
                
                clusterAssociationVector.push_back(clusterEndpointAssociation);
            }

            PandoraMonitoringApi::ViewEvent(this->GetPandora());
        }

        if (!clusterAssociationVector.empty())
        {
            /////////////////
            std::cout << "clusterAssociationVector size: " << clusterAssociationVector.size() << std::endl; 
            for (const ClusterEndpointAssociation &clusterEndpointAssociation : clusterAssociationVector)
            {
                ClusterList mainCluster({clusterEndpointAssociation.GetMainTrackCluster()});
                const CartesianVector &upstream(clusterEndpointAssociation.GetUpstreamMergePoint());
                const CartesianVector &downstream(clusterEndpointAssociation.GetDownstreamMergePoint());

                PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &upstream, "UPSTREAM", VIOLET, 2);
                PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &downstream, "DOWNSTREAM", RED, 2);
                PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &mainCluster, "CLUSTER", BLACK);
                PandoraMonitoringApi::ViewEvent(this->GetPandora());
                /////////////////
            }
            
            break;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayEndpointCorrectionAlgorithm::GetExtrapolatedCaloHits(const ClusterEndpointAssociation &clusterAssociation, const ClusterList *const pClusterList, ClusterToCaloHitListMap &clusterToCaloHitListMap) const
{
    // Fill subset information to try and minimise loops
    const CartesianVector &upstreamPoint(clusterAssociation.GetUpstreamMergePoint()), &downstreamPoint(clusterAssociation.GetDownstreamMergePoint());

    ClusterVector clustersInRegion;
    ClusterToCaloHitListMap hitsInRegion;
    for (const Cluster *const pCluster : *pClusterList)
    {
        const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
        for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
        {
            for (const CaloHit *const pCaloHit : *mapEntry.second)
            {
                const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
                
                if(!IsInLineSegment(upstreamPoint, downstreamPoint, hitPosition))
                    continue;

                if (std::find(clustersInRegion.begin(), clustersInRegion.end(), pCluster) == clustersInRegion.end())
                    clustersInRegion.push_back(pCluster);
                
                hitsInRegion[pCluster].push_back(pCaloHit);
            }
        }
    }

    for (const Cluster *const pCluster : clustersInRegion)
    {
        for (const CaloHit *const pCaloHit : hitsInRegion.at(pCluster))
        {
            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "HIt", BLUE, 2);
        }
    }

    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    
    float m_initialLengthRequired(40.f);

    const CartesianVector clusterMergePosition(clusterAssociation.IsEndUpstream() ? clusterAssociation.GetDownstreamMergePoint() : clusterAssociation.GetUpstreamMergePoint());
    const CartesianVector clusterMergeDirection(clusterAssociation.IsEndUpstream() ? clusterAssociation.GetDownstreamMergeDirection() : clusterAssociation.GetUpstreamMergeDirection());
    const CartesianVector clusterStartPosition(clusterMergePosition + (clusterMergeDirection * (-1.f) * m_initialLengthRequired));

    const float minX(std::min(clusterMergePosition.GetX(), clusterStartPosition.GetX())), maxX(std::max(clusterMergePosition.GetX(), clusterStartPosition.GetX()));
    const float minZ(std::min(clusterMergePosition.GetZ(), clusterStartPosition.GetZ())), maxZ(std::max(clusterMergePosition.GetZ(), clusterStartPosition.GetZ()));

    CartesianPointVector hitPositionVector;
 
    const OrderedCaloHitList &orderedCaloHitList(clusterAssociation.GetMainTrackCluster()->GetOrderedCaloHitList());
    for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
    {
        for (const CaloHit *const pCaloHit : *mapEntry.second)
        {
            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            if ((hitPosition.GetX() < minX) || (hitPosition.GetX() > maxX) || (hitPosition.GetZ() < minZ) || (hitPosition.GetZ() > maxZ))
                continue;

            hitPositionVector.push_back(hitPosition);    
        }
    }

    for (auto entry : hitPositionVector)
    {
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &entry, "PICKED POINT", VIOLET, 2);
    }

    PandoraMonitoringApi::ViewEvent(this->GetPandora());

    unsigned int count(0);
    unsigned int hitsCollected;
    do
    {
        hitsCollected = 0;
        
        try
        {
            const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
            const int m_slidingFitWindow(20);         
            const TwoDSlidingFitResult extrapolatedFit(&hitPositionVector, m_slidingFitWindow, slidingFitPitch);

            CartesianVector startPosition(0.f,0.f,0.f), startDirection(0.f,0.f,0.f);

            if (count == 0)
            {
                startPosition = (clusterAssociation.IsEndUpstream() ? clusterAssociation.GetDownstreamMergePoint() : clusterAssociation.GetUpstreamMergePoint());
                startDirection = (clusterAssociation.IsEndUpstream() ? clusterAssociation.GetDownstreamMergeDirection() : clusterAssociation.GetUpstreamMergeDirection());
            }
            else
            {
                startPosition = clusterAssociation.IsEndUpstream() ? extrapolatedFit.GetGlobalMinLayerPosition() : extrapolatedFit.GetGlobalMaxLayerPosition();
                startDirection = clusterAssociation.IsEndUpstream() ? extrapolatedFit.GetGlobalMinLayerDirection() * (-1.f) : extrapolatedFit.GetGlobalMaxLayerDirection();
            }

            const CartesianVector endPosition(startPosition + (startDirection * 5.f));

            const float gradient((endPosition.GetZ() - startPosition.GetZ()) / (endPosition.GetX() - startPosition.GetX()));

            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &startPosition, "start", RED, 2);
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &endPosition, "end", RED, 2);

            for (const Cluster *const pCluster : clustersInRegion)
            {
                for (const CaloHit *const pCaloHit : hitsInRegion.at(pCluster))
                {
                    const CartesianVector &hitPosition(pCaloHit->GetPositionVector());

                    const ClusterToCaloHitListMap::iterator iter(clusterToCaloHitListMap.find(pCluster));
                    if (iter != clusterToCaloHitListMap.end())
                    {
                        if (std::find(iter->second.begin(), iter->second.end(), pCaloHit) != iter->second.end())
                            continue;
                    }
                    
                    if (!IsInLineSegment(startPosition, endPosition, hitPosition))
                        continue;

                    const float &hitWidth(pCaloHit->GetCellSize1());
                    const CartesianVector highHitEdge(hitPosition.GetX() + (hitWidth * 0.5f), 0, hitPosition.GetZ());
                    const CartesianVector lowHitEdge(hitPosition.GetX() - (hitWidth * 0.5f), 0, hitPosition.GetZ());

                    //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &highHitEdge, "HITENDS", BLACK, 2);
                    //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &lowHitEdge, "HITENDS", BLACK, 2);
                    
                    const float highEdgeDistanceFromLine((endPosition - startPosition).GetCrossProduct(highHitEdge - startPosition).GetMagnitude());
                    const float lowEdgeDistanceFromLine((endPosition - startPosition).GetCrossProduct(lowHitEdge - startPosition).GetMagnitude());

                    // PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "COLLECTED", VIOLET, 2);
                    //std::cout << "highEdgeDistanceFromLine: " << highEdgeDistanceFromLine << std::endl;  
                    //std::cout << "lowEdgeDistanceFromLine: " << lowEdgeDistanceFromLine << std::endl;

                    //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "COLLECTED", GREEN, 2);
                    
                    if ((highEdgeDistanceFromLine > 10.0f) || (lowEdgeDistanceFromLine > 10.f))
                    {
                        //std::cout << "Too wide" << std::endl;
                        //PandoraMonitoringApi::Pause(this->GetPandora());
                        continue;
                    }

                    float xOnLine(((hitPosition.GetZ() - startPosition.GetZ()) / gradient) + startPosition.GetX());
                    CartesianVector positionOnLine(xOnLine, 0, hitPosition.GetZ());
                    //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &positionOnLine, "positionOnLine", BLUE, 2);

                    if (((highHitEdge.GetX() > xOnLine) && (lowHitEdge.GetX() > xOnLine)) ||
                        ((highHitEdge.GetX() < xOnLine) && (lowHitEdge.GetX() < xOnLine)))
                    {
                    
                        if (!((highEdgeDistanceFromLine < 0.5f) || (lowEdgeDistanceFromLine < 0.5f)))
                        {
                            //std::cout << "doesn't straddle and not clse enough" << std::endl;
                            //PandoraMonitoringApi::Pause(this->GetPandora());
                            continue;
                        }
                    }

                    //std::cout << "hit added" << std::endl;

                    //PandoraMonitoringApi::Pause(this->GetPandora());

                    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "COLLECTED", GREEN, 2);
                    
                    hitPositionVector.push_back(hitPosition);
                    ++hitsCollected;
                    clusterToCaloHitListMap[pCluster].push_back(pCaloHit);
                }
            }

            std::cout << "hitsCollected: " << hitsCollected << std::endl;

        }
        catch (const StatusCodeException &)
        {
            std::cout << "FIT FAILED IN EXTRAPOLATED HITS" << std::endl;
            break;
        }

        ++count;
    }
    while (hitsCollected > 0);

    

    PandoraMonitoringApi::ViewEvent(this->GetPandora());
        /*

    for (const Cluster *const pCluster : *pClusterList)
    {
        const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
        for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
        {
            for (const CaloHit *const pCaloHit : *mapEntry.second)
            {
                const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
                if ((hitPosition.GetX() < minX) || (hitPosition.GetX() > maxX) || (hitPosition.GetZ() < upstreamPoint.GetZ()) || (hitPosition.GetZ() > downstreamPoint.GetZ()))
                    continue;

                const float distanceFromLine(clusterAssociation.GetConnectingLineDirection().GetCrossProduct(hitPosition - upstreamPoint).GetMagnitude());
                if (distanceFromLine > m_distanceFromLine)
                    continue;

                clusterToCaloHitListMap[pCluster].push_back(pCaloHit);
            }
        }
    }
    */

    /*
    const CartesianVector &upstreamPoint(clusterAssociation.GetUpstreamMergePoint()), &downstreamPoint(clusterAssociation.GetDownstreamMergePoint());
    const float minX(std::min(upstreamPoint.GetX(), downstreamPoint.GetX())), maxX(std::max(upstreamPoint.GetX(), downstreamPoint.GetX()));
        
    for (const Cluster *const pCluster : *pClusterList)
    {
        const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
        for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
        {
            for (const CaloHit *const pCaloHit : *mapEntry.second)
            {
                const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
                if ((hitPosition.GetX() < minX) || (hitPosition.GetX() > maxX) || (hitPosition.GetZ() < upstreamPoint.GetZ()) || (hitPosition.GetZ() > downstreamPoint.GetZ()))
                    continue;

                const float distanceFromLine(clusterAssociation.GetConnectingLineDirection().GetCrossProduct(hitPosition - upstreamPoint).GetMagnitude());
                if (distanceFromLine > m_distanceFromLine)
                    continue;

                clusterToCaloHitListMap[pCluster].push_back(pCaloHit);
            }
        }
    }
    */
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CosmicRayEndpointCorrectionAlgorithm::IsDeltaRay(const Cluster *const pCluster, const CartesianVector &clusterMergePoint,
    const CartesianVector &clusterMergeDirection, const bool isEndUpstream, CartesianVector &mergePosition) const
{
    // MAKE MORE DETAILED FIT OF THE AMBIGUOUS SECTION
    CartesianPointVector hitSubset;
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
    for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
    {
        for (const CaloHit *const pCaloHit : *mapEntry.second) 
        {
            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());

            if ((isEndUpstream && (hitPosition.GetZ() < clusterMergePoint.GetZ())) || (!isEndUpstream && (hitPosition.GetZ() > clusterMergePoint.GetZ())))
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

    // INVESTIGATE THE DIRECTION CHANGE
    const TwoDSlidingFitResult subsetFit(subsetFitVector.front());
    const LayerFitResultMap &clusterMicroLayerFitResultMap(subsetFit.GetLayerFitResultMap());
    const int startLayer(isEndUpstream ? subsetFit.GetMaxLayer() : subsetFit.GetMinLayer());
    const int endLayer(isEndUpstream ? subsetFit.GetMinLayer() : subsetFit.GetMaxLayer());
    const int loopTerminationLayer(isEndUpstream ? endLayer - 1 : endLayer + 1);
    const int step(isEndUpstream ? -1 : 1);

    //unsigned int anomalousLayerCount(0);    
    //bool reachedFirstCurve(false), isCurveClockwise(false);
    float previousOpeningAngle(std::numeric_limits<float>::max());
    
    for (int i = startLayer; i != loopTerminationLayer; i += step)
    {
        const auto microIter(clusterMicroLayerFitResultMap.find(i));

        if (microIter == clusterMicroLayerFitResultMap.end())
            continue;

        CartesianVector microDirection(0.f, 0.f, 0.f);
        subsetFit.GetGlobalDirection(microIter->second.GetGradient(), microDirection);
        float microOpeningAngle(microDirection.GetOpeningAngle(clusterMergeDirection) * 180 / 3.14);

        // IT IS USING THE CORRECT DEFINITION - AT THIS POINT ALL DIRECTIONS HAVE POSITIVE Z?
        if(microDirection.GetZ() < (clusterMergeDirection.GetZ() * microDirection.GetX() / clusterMergeDirection.GetX()))
            microOpeningAngle *= (-1.f);

        const float layerAngleDeviation(previousOpeningAngle > 180.f ? microOpeningAngle : std::fabs(microOpeningAngle - previousOpeningAngle));
        
        /////////////////////
        CartesianVector microPosition(0.f, 0.f, 0.f);
        subsetFit.GetGlobalPosition(microIter->second.GetL(), microIter->second.GetFitT(), microPosition);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &microPosition, "MICRO POSITION", BLACK, 2);
        std::cout << "ANGLE DEVIATION: " << microOpeningAngle << std::endl;
        std::cout << "<-------------------- CHANGE FROM LAST LAYER: " << layerAngleDeviation << std::endl;
        /////////////////////

        if (std::fabs(microOpeningAngle) > m_thresholdMaxAngleDeviation)
        {
            std::cout << "MEET ANGLE CRITERIA" << std::endl;

            mergePosition = (isEndUpstream ? subsetFit.GetGlobalMaxLayerPosition() : subsetFit.GetGlobalMinLayerPosition());
            return true;
        } 

        /*
        // ISOBEL - DO YOU NEED FINAL THRESHOLD
        // CHECK IT DOESNT CHANGE SIGN WITH SIGNIFICANCE.
        if ((!reachedFirstCurve) && (std::fabs(microOpeningAngle) > m_thresholdAngleDeviation))
        {
            reachedFirstCurve = true;
            isCurveClockwise = (microOpeningAngle > 0.f);
            continue;
        }

        if (reachedFirstCurve)
        {
            if ((isCurveClockwise && (microOpeningAngle < previousOpeningAngle)) || (!isCurveClockwise && (microOpeningAngle > previousOpeningAngle)) ||
                (layerAngleDeviation < m_thresholdAngleDeviationBetweenLayers))
            {
                ++anomalousLayerCount;

                if (anomalousLayerCount > m_maxAnomalousPoints)
                {
                    std::cout << "TOO SMOOTH" << std::endl;
                    PandoraMonitoringApi::ViewEvent(this->GetPandora());
                    //return false;
                    reachedFirstCurve = false;
                    continue;
                }
            }
            else
            {
                anomalousLayerCount = 0;
            }
            
                if (std::fabs(microOpeningAngle) > m_thresholdMaxAngleDeviation)
                {
                    std::cout << "MEET ANGLE CRITERIA" << std::endl;

                    mergePosition = (isEndUpstream ? subsetFit.GetGlobalMaxLayerPosition() : subsetFit.GetGlobalMinLayerPosition());
                    return true;
                }

        
        }
        */

        previousOpeningAngle = microOpeningAngle;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------ 

void CosmicRayEndpointCorrectionAlgorithm::CreateMainTrack(ClusterEndpointAssociation &clusterEndpointAssociation, const ClusterToCaloHitListMap &clusterToCaloHitListMap, const ClusterList *const pClusterList, ClusterVector &clusterVector, SlidingFitResultMapPair &slidingFitResultMapPair) const
{
    const Cluster *pMainTrackCluster(clusterEndpointAssociation.GetMainTrackCluster());
    const CartesianVector &clusterMergePoint(clusterEndpointAssociation.IsEndUpstream() ?
        clusterEndpointAssociation.GetDownstreamMergePoint() : clusterEndpointAssociation.GetUpstreamMergePoint());

    ////////////////
    ClusterList originalTrack({pMainTrackCluster});
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &originalTrack, "ORIGINAL TRACK", BLACK);
    ///////////////
    
    // Determine the shower clusters which contain hits that belong to the main track
    ClusterVector showerClustersToFragment;
    for (auto &entry : clusterToCaloHitListMap)
    {
        if (entry.first != pMainTrackCluster) //COULD TURN THIS INTO A VECTOR - TO EXTEND TO EM SHOWER
            showerClustersToFragment.push_back(entry.first);
    }

    std::sort(showerClustersToFragment.begin(), showerClustersToFragment.end(), LArClusterHelper::SortByNHits);

    ////////////////
    for (auto &entry : showerClustersToFragment)
    {
        ClusterList jam({entry});
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &jam, "SHOWER CLUSTER", VIOLET);
    }
    ////////////////    

    ClusterList remnantClusterList;
    
    pMainTrackCluster = RemoveOffAxisHitsFromTrack(pMainTrackCluster, clusterMergePoint, clusterEndpointAssociation.IsEndUpstream(),
        clusterToCaloHitListMap, remnantClusterList, *slidingFitResultMapPair.first, *slidingFitResultMapPair.second);
   
    ////////////////
    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    ClusterList refinedCluster({pMainTrackCluster});
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &refinedCluster, "REFINED MAIN TRACK", BLACK);
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &remnantClusterList, "REMNANT CLUSTERS", VIOLET);
    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    ////////////////   
   
    for (const Cluster *const pShowerCluster : showerClustersToFragment)
    {
        const CaloHitList &caloHitsToMerge(clusterToCaloHitListMap.at(pShowerCluster));
        this->AddHitsToMainTrack(pMainTrackCluster, pShowerCluster, caloHitsToMerge, clusterEndpointAssociation, remnantClusterList);
    }

    ////////////////       
    ClusterList extendedCluster({pMainTrackCluster});
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &extendedCluster, "REFINED MAIN TRACK", BLACK);
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &remnantClusterList, "REMNANT CLUSTERS", VIOLET);
    PandoraMonitoringApi::ViewEvent(this->GetPandora());    
    ////////////////   
    
    ClusterList createdClusters;
    this->ProcessRemnantClusters(remnantClusterList, pMainTrackCluster, pClusterList, createdClusters);

    ////////////////
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &createdClusters, "CREATED CLUSTERS", RED);
    PandoraMonitoringApi::ViewEvent(this->GetPandora());     
    ////////////////       

    // ATTN: Update containers and the cluster association
    ClusterList modifiedClusters(showerClustersToFragment.begin(), showerClustersToFragment.end());   
    this->UpdateContainers(createdClusters, modifiedClusters, clusterVector, slidingFitResultMapPair);

    //ATTN: Update references to pMainTrackCluster 
    this->UpdateAfterMainTrackModification(pMainTrackCluster, clusterEndpointAssociation, slidingFitResultMapPair);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayEndpointCorrectionAlgorithm::UpdateAfterMainTrackModification(const Cluster *const pMainTrackCluster, ClusterEndpointAssociation &clusterEndpointAssociation, SlidingFitResultMapPair &slidingFitResultMapPair) const
{

    const Cluster *const pDeletedTrackCluster(clusterEndpointAssociation.GetMainTrackCluster());
    
    // REPLACE IN MICRO FIT MAP
    const TwoDSlidingFitResultMap::const_iterator microFitToDeleteIter(slidingFitResultMapPair.first->find(pDeletedTrackCluster));
    if (microFitToDeleteIter == slidingFitResultMapPair.first->end())
    {
        std::cout << "ISOBEL THIS IS BAD" << std::endl;
        throw STATUS_CODE_FAILURE;
    }
    
    slidingFitResultMapPair.first->insert(std::make_pair(pMainTrackCluster, microFitToDeleteIter->second));
    slidingFitResultMapPair.first->erase(microFitToDeleteIter);

    // REPLACE IN MACRO FIT MAP
    const TwoDSlidingFitResultMap::const_iterator macroFitToDeleteIter(slidingFitResultMapPair.second->find(pDeletedTrackCluster));
    if (macroFitToDeleteIter == slidingFitResultMapPair.second->end())
    {
        std::cout << "ISOBEL THIS IS BAD" << std::endl;
        throw STATUS_CODE_FAILURE;
    }
    
    slidingFitResultMapPair.second->insert(std::make_pair(pMainTrackCluster, macroFitToDeleteIter->second));
    slidingFitResultMapPair.second->erase(macroFitToDeleteIter);

    clusterEndpointAssociation.SetMainTrackCluster(pMainTrackCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayEndpointCorrectionAlgorithm::RemoveClusterAssociationFromClusterVector(const ClusterEndpointAssociation &clusterAssociation, ClusterVector &clusterVector) const
{
    const ClusterVector::const_iterator clusterToDeleteIter(std::find(clusterVector.begin(), clusterVector.end(), clusterAssociation.GetMainTrackCluster()));
    if (clusterToDeleteIter != clusterVector.end())
        clusterVector.erase(clusterToDeleteIter);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayEndpointCorrectionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCaloHits", m_minCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxDistanceFromTPC", m_maxDistanceFromTPC));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CurveThreshold", m_curveThreshold));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinScaledZOffset", m_minScaledZOffset));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ThresholdAngleDeviation", m_thresholdAngleDeviation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ThresholdAngleDeviationBetweenLayers", m_thresholdAngleDeviationBetweenLayers));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxAnomalousPoints", m_maxAnomalousPoints));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ThresholdMaxAngleDeviation", m_thresholdMaxAngleDeviation));           

    return CosmicRayTrackRefinementBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content



/*
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
    
    if (averageDistanceFromLine < m_curveThreshold)
    {
        std::cout << "TOO STRAIGHT" << std::endl;
        PandoraMonitoringApi::ViewEvent(this->GetPandora());        
        return false;
    }
*/




///OLD METHOD
/*

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
        
        if (((microOpeningAngle > m_thresholdAngleDeviation) || (microOpeningAngle < -m_thresholdAngleDeviation)) && !reachedFirstCurve)
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

                if (numberOfWobbles > m_maxNumberOfWobbles)
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



 */
