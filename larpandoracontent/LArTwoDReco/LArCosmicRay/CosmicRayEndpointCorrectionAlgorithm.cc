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
    m_minScaledZOffset(2.f), //0.25 m_curveThreshold(0.3f),
    m_thresholdAngleDeviation(10.f),
    m_thresholdAngleDeviationBetweenLayers(1.f),
    m_maxAnomalousPoints(2),
    m_thresholdMaxAngleDeviation(21.f),
    m_deltaRayslidingFitWindow(5),
    m_growingFitInitialLength(40.f),
    m_growingFitSegmentLength(7.f),
    m_furthestDistanceToLine(10.f),
    m_closestDistanceToLine(0.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------    
    
void CosmicRayEndpointCorrectionAlgorithm::FindBestClusterAssociation(const ClusterVector &clusterVector, const SlidingFitResultMapPair &slidingFitResultMapPair,
    ClusterAssociationVector &clusterAssociationVector)
{
    // Collect LArTPC information
    const Pandora *pPrimaryPandoraInstance;
    try
    {
        pPrimaryPandoraInstance = MultiPandoraApi::GetPrimaryPandoraInstance(&this->GetPandora());
    }
    catch (const StatusCodeException &)
    {
        pPrimaryPandoraInstance = &this->GetPandora();
    }
           
    float detectorMinXEdge(std::numeric_limits<float>::max());
    float detectorMaxXEdge(-std::numeric_limits<float>::max());

    const LArTPC *pLArTPC(&this->GetPandora().GetGeometry()->GetLArTPC());
    const LArTPCMap &larTPCMap(pPrimaryPandoraInstance->GetGeometry()->GetLArTPCMap());
    
    for (const LArTPCMap::value_type &mapEntry : larTPCMap)
    {
        const LArTPC *const pSubLArTPC(mapEntry.second);
        detectorMinXEdge = std::min(detectorMinXEdge, pSubLArTPC->GetCenterX() - 0.5f * pSubLArTPC->GetWidthX());
        detectorMaxXEdge = std::max(detectorMaxXEdge, pSubLArTPC->GetCenterX() + 0.5f * pSubLArTPC->GetWidthX());

        //ATTN: Child & parent pandora instance TPCs have different addresses
        if (std::fabs(pSubLArTPC->GetCenterX() - pLArTPC->GetCenterX()) < std::numeric_limits<float>::epsilon())
            pLArTPC = pSubLArTPC;
    }
    
    //ATTN: Find best cluster association - method assumes clusters ordered by their hits
    
    for (const Cluster *const pCluster : clusterVector)
    {
        const TwoDSlidingFitResult &microSlidingFitResult(slidingFitResultMapPair.first->at(pCluster));
        const TwoDSlidingFitResult &macroSlidingFitResult(slidingFitResultMapPair.second->at(pCluster));

        // ATTN: Investigate both ends of the track
        for (unsigned int isEndUpstream = 0; isEndUpstream < 2; ++isEndUpstream)
        {
            /////////////////
            /*
            ClusterList theCluster({pCluster});
            PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &theCluster, "CONSIDERED CLUSTER", BLACK);
            std::cout << "isEndUpstream: " << isEndUpstream << std::endl;
            */
            ////////////////

            // ATTN: Match the cluster to itself to get the cluster merging points
            const bool isClusterUpstream(!isEndUpstream);
            
            CartesianVector clusterMergePoint(0.f, 0.f, 0.f), clusterMergeDirection(0.f, 0.f, 0.f);
            if (!GetClusterMergingCoordinates(microSlidingFitResult, macroSlidingFitResult, macroSlidingFitResult, isClusterUpstream, clusterMergePoint, clusterMergeDirection))
            {
                //std::cout << "CANNOT FIND MERGE POSITION" << std::endl;
                //PandoraMonitoringApi::ViewEvent(this->GetPandora());
                continue;
            }
            
            const CartesianVector &endpointPosition(isEndUpstream ? microSlidingFitResult.GetGlobalMinLayerPosition() : microSlidingFitResult.GetGlobalMaxLayerPosition());
            const float endpointSeparation((endpointPosition - clusterMergePoint).GetMagnitude());
            
            /////////////////
            /*
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &endpointPosition, "ENDPOINT", BLUE, 2);
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &clusterMergePoint, "MERGE POINT", BLUE, 2);
            std::cout << "Endpoint Separation: " << endpointSeparation << std::endl;
            */
            /////////////////    

            if (endpointSeparation < std::numeric_limits<float>::epsilon())
            {
                //std::cout << "MERGE POINT AND ENDPOINT ARE THE SAME" << std::endl;
                //PandoraMonitoringApi::ViewEvent(this->GetPandora());
                continue;
            }

            // INVESTIGATE PROXIMITY TO TPC BOUNDARY
            const bool isNearestBoundaryHigherX((endpointPosition.GetX() - pLArTPC->GetCenterX()) > 0.f);
            const float nearestTPCBoundaryX(pLArTPC->GetCenterX() + (pLArTPC->GetWidthX() * 0.5f * (isNearestBoundaryHigherX ? 1.f : -1.f)));

            if ((std::fabs(nearestTPCBoundaryX - detectorMinXEdge) < std::numeric_limits<float>::epsilon()) ||
                (std::fabs(nearestTPCBoundaryX - detectorMaxXEdge) < std::numeric_limits<float>::epsilon()))
            {
                //std::cout << "ENDPOINT IS NEAR DETECTOR EDGE" << std::endl;
                //PandoraMonitoringApi::ViewEvent(this->GetPandora());
                continue;
            }

            // is it closer to the boundary then the other endpoint
            const CartesianVector &otherEndpointPosition(isEndUpstream ? microSlidingFitResult.GetGlobalMaxLayerPosition() : microSlidingFitResult.GetGlobalMinLayerPosition());

            if (std::fabs(endpointPosition.GetX() - nearestTPCBoundaryX) > std::fabs(otherEndpointPosition.GetX() - nearestTPCBoundaryX))
            {
                //std::cout << "OTHER ENDPOINT IS CLOSE TO THIS BOUNDARY" << std::endl;
                //PandoraMonitoringApi::ViewEvent(this->GetPandora());
                continue;
            }

            // this doesn't take into account the gap.
            const float allowanceBoundaryX(nearestTPCBoundaryX + (m_maxDistanceFromTPC * (isNearestBoundaryHigherX ? -1.f : 1.f)));
            const bool isMergePointInAllowanceRegion(isNearestBoundaryHigherX ? clusterMergePoint.GetX() > allowanceBoundaryX : clusterMergePoint.GetX() < allowanceBoundaryX);
            const bool isEndpointInAllowanceRegion(isNearestBoundaryHigherX ? endpointPosition.GetX() > allowanceBoundaryX : endpointPosition.GetX() < allowanceBoundaryX);

            if (!(isMergePointInAllowanceRegion || isEndpointInAllowanceRegion))
            {
                //std::cout << "NOT CLOSE ENOUGH TO THE BOUNDARY" << std::endl;
                //PandoraMonitoringApi::ViewEvent(this->GetPandora());
                continue;
            }

            // INVESTIGATE ENDPOINT Z SEPARATION
            if (std::fabs(clusterMergeDirection.GetX()) < std::numeric_limits<float>::epsilon())
            {
                //std::cout << "MERGE DIRECTION HAS NO X COMPONENT" << std::endl;
                //PandoraMonitoringApi::ViewEvent(this->GetPandora());
                continue;
            }
            
            const float predictedGradient(clusterMergeDirection.GetZ() / clusterMergeDirection.GetX());
            const float predictedIntercept(clusterMergePoint.GetZ() - (predictedGradient * clusterMergePoint.GetX()));
            const CartesianVector fitEndpointPosition(endpointPosition.GetX(), 0.f, predictedIntercept + (predictedGradient * endpointPosition.GetX()));
            const float deltaZ(std::fabs(endpointPosition.GetZ() - fitEndpointPosition.GetZ()));
            //const float scaledDeltaZ(deltaZ / endpointSeparation);            

            ///////////////////
            /*
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &fitEndpointPosition, "FIT ENDPOINT", ORANGE, 2);
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &clusterMergePoint, "clusterMergePoint", VIOLET, 2);
            const CartesianVector start(clusterMergePoint + (clusterMergeDirection*40));
            const CartesianVector end(clusterMergePoint - (clusterMergeDirection*40));
            PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &start, &end, "CLUSTER DIRECTION", GREEN, 4, 2);            
            std::cout << "deltaZ: " << deltaZ << std::endl;
            std::cout << "scaled deltaZ: " << scaledDeltaZ << std::endl;  
            */
            /////////////////

            if (deltaZ < m_minScaledZOffset)
            {
                //std::cout << "DELTA Z IS TOO LOW" << std::endl;
                //PandoraMonitoringApi::ViewEvent(this->GetPandora());
                continue;
            }
           
            // SEARCH FOR DELTA RAY BEND    
            if(this->IsDeltaRay(pCluster, clusterMergePoint, clusterMergeDirection, isEndUpstream))
            {
                // ATTN: Temporarily set the other merge point to define region to search for in extrapolate hits
                const LArTPC *const pNeighboughTPC(&LArStitchingHelper::FindClosestTPC(*pPrimaryPandoraInstance, *pLArTPC, isNearestBoundaryHigherX));
                
                const float gapSizeX(std::fabs(pNeighboughTPC->GetCenterX() - pLArTPC->GetCenterX()) - (pNeighboughTPC->GetWidthX() * 0.5f) - (pLArTPC->GetWidthX() * 0.5f));
                const float endpointX(nearestTPCBoundaryX + (gapSizeX * 0.5f * (isNearestBoundaryHigherX ? 1.f : -1.f)));
                
                const CartesianVector extrapolatedEndpointPosition(endpointX, 0.f, predictedIntercept + (predictedGradient * endpointX));

                /*
                //Check whether endpoint is in the detector
                const float detectorMinZEdge = pLArTPC->GetCenterZ() - 0.5f * pLArTPC->GetWidthZ();                
                const float detectorMaxZEdge = pLArTPC->GetCenterZ() + 0.5f * pLArTPC->GetWidthZ();

                if ((extrapolatedEndpointPosition.GetZ() < detectorMinZEdge) || (extrapolatedEndpointPosition.GetZ() > detectorMaxZEdge))
                {
                    std::cout << "EXTRAPOLATED ENDPOINT IS NOT IN THE DETECTOR" << std::endl;
                    PandoraMonitoringApi::ViewEvent(this->GetPandora());
                    continue;
                }
                */

                if (isEndUpstream ? extrapolatedEndpointPosition.GetZ() > clusterMergePoint.GetZ() : extrapolatedEndpointPosition.GetZ() < clusterMergePoint.GetZ())
                {
                    //std::cout << "EXTRAPOLATED ENDPOINT IS NOT IN FORWARD DIRECTION" << std::endl;
                    //PandoraMonitoringApi::ViewEvent(this->GetPandora());
                    continue;
                }

                //std::cout << "LArClusterHelper::GetAverageHitSeparation: " << LArClusterHelper::GetAverageHitSeparation(pCluster) << std::endl;
                //std::cout << "Cluster Hits: " <<  pCluster->GetNCaloHits() << std::endl;
                //std::cout << "Cluster Hits / Length: " << pCluster->GetNCaloHits() / LArClusterHelper::GetLength(pCluster) << std::endl;
                //std::cout << "distance from line: " << GetAverageDeviationFromLine(pCluster, clusterMergeDirection, clusterMergePoint) << std::endl;
                
                ClusterEndpointAssociation clusterEndpointAssociation(isEndUpstream ?
                    ClusterEndpointAssociation(extrapolatedEndpointPosition, clusterMergeDirection, clusterMergePoint, clusterMergeDirection * (-1.f), pCluster, true) :
                    ClusterEndpointAssociation(clusterMergePoint, clusterMergeDirection, extrapolatedEndpointPosition, clusterMergeDirection * (-1.f), pCluster, false));
                
                clusterAssociationVector.push_back(clusterEndpointAssociation);
            }

            //PandoraMonitoringApi::ViewEvent(this->GetPandora());
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
            }
            
            /////////////////
            break;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CosmicRayEndpointCorrectionAlgorithm::IsDeltaRay(const Cluster *const pCluster, CartesianVector &clusterMergePoint,
    const CartesianVector &clusterMergeDirection, const bool isEndUpstream) const 
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
        const TwoDSlidingFitResult subsetFit(&hitSubset, m_deltaRayslidingFitWindow, slidingFitPitch);

        subsetFitVector.push_back(subsetFit);
    }
    catch (const StatusCodeException &)
    {
        //std::cout << "CANNOT MAKE A FIT" << std::endl;
        //PandoraMonitoringApi::ViewEvent(this->GetPandora());  
        return false;
    }

    // INVESTIGATE THE DIRECTION CHANGE
    const TwoDSlidingFitResult subsetFit(subsetFitVector.front());
    const LayerFitResultMap &clusterMicroLayerFitResultMap(subsetFit.GetLayerFitResultMap());
    int startLayer(isEndUpstream ? subsetFit.GetMaxLayer() : subsetFit.GetMinLayer());
    int endLayer(isEndUpstream ? subsetFit.GetMinLayer() : subsetFit.GetMaxLayer());
    int loopTerminationLayer(isEndUpstream ? endLayer - 1 : endLayer + 1);
    int step(isEndUpstream ? -1 : 1);
    
    unsigned int anomalousLayerCount(0);
    bool reachedFirstCurve(false);
    float previousOpeningAngle(std::numeric_limits<float>::max());

    // CHECK MERGE POINT END
    for (int i = startLayer; i != loopTerminationLayer; i += step)
    {
        const auto microIter(clusterMicroLayerFitResultMap.find(i));

        if (microIter == clusterMicroLayerFitResultMap.end())
            continue;

        CartesianVector microDirection(0.f, 0.f, 0.f);
        subsetFit.GetGlobalDirection(microIter->second.GetGradient(), microDirection);
        float microOpeningAngle(microDirection.GetOpeningAngle(clusterMergeDirection) * 180 / 3.14);

        if(microDirection.GetZ() < (clusterMergeDirection.GetZ() * microDirection.GetX() / clusterMergeDirection.GetX()))
            microOpeningAngle *= (-1.f);

        if (!reachedFirstCurve)
        {
            if (std::fabs(microOpeningAngle) > m_thresholdAngleDeviation)
                reachedFirstCurve = true;
            
            previousOpeningAngle = microOpeningAngle;
            continue;
        }

        const float layerAngleDeviation(std::fabs(microOpeningAngle - previousOpeningAngle));

        if ((std::fabs(microOpeningAngle) < std::fabs(previousOpeningAngle)) || (microOpeningAngle * previousOpeningAngle < 0.f) ||
            (layerAngleDeviation < m_thresholdAngleDeviationBetweenLayers))
        {
            ++anomalousLayerCount;

            if (anomalousLayerCount > m_maxAnomalousPoints)
            {
                //std::cout << "MERGE POINT END TOO SMOOTH" << std::endl;
                break;
            }
        }
        else
        {
            if (std::fabs(microOpeningAngle) > m_thresholdMaxAngleDeviation)
            {
                std::cout << "MERGE POINT END MEET ANGLE CRITERIA" << std::endl;
                // ATTN: Make cluster merge points more precise - NEED CLUSTER FIT START
                clusterMergePoint = isEndUpstream ? subsetFit.GetGlobalMaxLayerPosition() : subsetFit.GetGlobalMinLayerPosition();
                
                return true;
            }

            anomalousLayerCount = 0;
        }

        previousOpeningAngle = microOpeningAngle;
    }

    // CHECK CLUSTER ENDPOINT END
    CartesianVector endDirection(isEndUpstream ? subsetFit.GetGlobalMinLayerDirection() : subsetFit.GetGlobalMaxLayerDirection());
    float endOpeningAngle(endDirection.GetOpeningAngle(clusterMergeDirection) * 180 / 3.14);
    if (std::fabs(endOpeningAngle) < m_thresholdMaxAngleDeviation)
    {
        /*
        std::cout << "END DEVIATION ANGLE NOT BIG ENOUGH" << std::endl;
        std::cout << "endOpeningAngle: " << endOpeningAngle << std::endl;
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        */
        return false;
    }
    
    startLayer = (isEndUpstream ? subsetFit.GetMinLayer() : subsetFit.GetMaxLayer());
    endLayer = (isEndUpstream ? subsetFit.GetMaxLayer() : subsetFit.GetMinLayer());
    loopTerminationLayer = (isEndUpstream ? endLayer + 1 : endLayer - 1);
    step = (isEndUpstream ? 1 : -1);
    anomalousLayerCount = 0;    
    reachedFirstCurve = false;

    for (int i = startLayer; i != loopTerminationLayer; i += step)
    {
        const auto microIter(clusterMicroLayerFitResultMap.find(i));

        if (microIter == clusterMicroLayerFitResultMap.end())
            continue;

        CartesianVector microDirection(0.f, 0.f, 0.f);
        subsetFit.GetGlobalDirection(microIter->second.GetGradient(), microDirection);
        float microOpeningAngle(microDirection.GetOpeningAngle(clusterMergeDirection) * 180 / 3.14);

        if(microDirection.GetZ() < (clusterMergeDirection.GetZ() * microDirection.GetX() / clusterMergeDirection.GetX()))
            microOpeningAngle *= (-1.f);

        if (!reachedFirstCurve)
        {
            if (std::fabs(microOpeningAngle) < m_thresholdMaxAngleDeviation)
                reachedFirstCurve = true;
            
            previousOpeningAngle = microOpeningAngle;
            continue;
        }

        const float layerAngleDeviation(std::fabs(microOpeningAngle - previousOpeningAngle));

        if ((std::fabs(microOpeningAngle) > std::fabs(previousOpeningAngle)) || (microOpeningAngle * previousOpeningAngle < 0.f) ||
            (layerAngleDeviation < m_thresholdAngleDeviationBetweenLayers))
        {
            ++anomalousLayerCount;

            if (anomalousLayerCount > m_maxAnomalousPoints)
            {
                //std::cout << "ENDPOINT END TOO SMOOTH" << std::endl;
                break;
            }
        }
        else
        {
            if (std::fabs(microOpeningAngle) < m_thresholdAngleDeviation)
            {
                std::cout << "ENDPOINT END MEET ANGLE CRITERIA" << std::endl;
                
                // ATTN: Make cluster merge points more precise - NEED CLUSTER FIT START
                clusterMergePoint = isEndUpstream ? subsetFit.GetGlobalMaxLayerPosition() : subsetFit.GetGlobalMinLayerPosition();
                return true;
            }

            anomalousLayerCount = 0;
        }

        previousOpeningAngle = microOpeningAngle;
    }    

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayEndpointCorrectionAlgorithm::GetExtrapolatedCaloHits(ClusterEndpointAssociation &clusterAssociation, const ClusterList *const pClusterList,
    ClusterToCaloHitListMap &clusterToCaloHitListMap) const
{

    // Look for clusters in the region of interest
    ClusterToCaloHitListMap hitsInRegion;
    for (const Cluster *const pCluster : *pClusterList)
    {
        const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
        for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
        {
            for (const CaloHit *const pCaloHit : *mapEntry.second)
            {
                const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
                
                if(!IsInLineSegment(clusterAssociation.GetUpstreamMergePoint(), clusterAssociation.GetDownstreamMergePoint(), hitPosition))
                    continue;
                
                hitsInRegion[pCluster].push_back(pCaloHit);
            }
        }
    }

    // ATTN: hitsInRegion is ordered (it isnt)
    ClusterVector clustersInRegion;
    for (const ClusterToCaloHitListMap::value_type &entry : hitsInRegion)
        clustersInRegion.push_back(entry.first);

    std::sort(clustersInRegion.begin(), clustersInRegion.end(), LArClusterHelper::SortByNHits);
    
    CartesianVector extrapolatedStartPosition(clusterAssociation.IsEndUpstream() ? clusterAssociation.GetDownstreamMergePoint() : clusterAssociation.GetUpstreamMergePoint());
    CartesianVector extrapolatedDirection(clusterAssociation.IsEndUpstream() ? clusterAssociation.GetDownstreamMergeDirection() : clusterAssociation.GetUpstreamMergeDirection());
    const CartesianVector clusterSubsetBoundary(extrapolatedStartPosition + (extrapolatedDirection * (-1.f) * m_growingFitInitialLength));

    const float minX(std::min(extrapolatedStartPosition.GetX(), clusterSubsetBoundary.GetX())), maxX(std::max(extrapolatedStartPosition.GetX(), clusterSubsetBoundary.GetX()));
    const float minZ(std::min(extrapolatedStartPosition.GetZ(), clusterSubsetBoundary.GetZ())), maxZ(std::max(extrapolatedStartPosition.GetZ(), clusterSubsetBoundary.GetZ()));

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
    
    unsigned int count(0);
    //unsigned int count(1);
    unsigned int hitsCollected(std::numeric_limits<int>::max());
    CartesianVector extrapolatedEndPosition(0.f, 0.f, 0.f);
    while (hitsCollected > 0)
    {
        hitsCollected = 0;
        
        try
        {
            const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
            const TwoDSlidingFitResult extrapolatedFit(&hitPositionVector, m_microSlidingFitWindow, slidingFitPitch);

            if (count > 0)
            {
                extrapolatedStartPosition = clusterAssociation.IsEndUpstream() ? extrapolatedFit.GetGlobalMinLayerPosition() : extrapolatedFit.GetGlobalMaxLayerPosition();
                extrapolatedDirection = clusterAssociation.IsEndUpstream() ? extrapolatedFit.GetGlobalMinLayerDirection() * (-1.f) : extrapolatedFit.GetGlobalMaxLayerDirection();
            }
            
            extrapolatedEndPosition = extrapolatedStartPosition + (extrapolatedDirection * m_growingFitSegmentLength);
            const float gradient((extrapolatedEndPosition.GetZ() - extrapolatedStartPosition.GetZ()) / (extrapolatedEndPosition.GetX() - extrapolatedStartPosition.GetX()));
                        
            ////////////////
            //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedStartPosition, "start", RED, 2);
            //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedEndPosition, "end", RED, 2);
            ////////////////

            for (const Cluster *const pCluster : clustersInRegion)
            {
                const CaloHitList &theCaloHits(hitsInRegion.at(pCluster));
                for (const CaloHit *const pCaloHit : theCaloHits)
                {
                    const CartesianVector &hitPosition(pCaloHit->GetPositionVector());

                    // ATTN: To avoid counting same hit twice
                    const ClusterToCaloHitListMap::iterator iter(clusterToCaloHitListMap.find(pCluster));
                    if (iter != clusterToCaloHitListMap.end())
                    {
                        if (std::find(iter->second.begin(), iter->second.end(), pCaloHit) != iter->second.end())
                            continue;
                    }

                    if (!IsInLineSegment(extrapolatedStartPosition, extrapolatedEndPosition, hitPosition))
                        continue;

                    const float &hitWidth(pCaloHit->GetCellSize1());
                    const CartesianVector hitHighEdge(hitPosition.GetX() + (hitWidth * 0.5f), 0, hitPosition.GetZ());
                    const CartesianVector hitLowEdge(hitPosition.GetX() - (hitWidth * 0.5f), 0, hitPosition.GetZ());
                    
                    const float highEdgeDistanceFromLine((extrapolatedEndPosition - extrapolatedStartPosition).GetCrossProduct(hitHighEdge - extrapolatedStartPosition).GetMagnitude());
                    const float lowEdgeDistanceFromLine((extrapolatedEndPosition - extrapolatedStartPosition).GetCrossProduct(hitLowEdge - extrapolatedStartPosition).GetMagnitude());

                    if ((highEdgeDistanceFromLine > m_furthestDistanceToLine) || (lowEdgeDistanceFromLine > m_furthestDistanceToLine))
                        continue;

                    float xOnLine(((hitPosition.GetZ() - extrapolatedStartPosition.GetZ()) / gradient) + extrapolatedStartPosition.GetX());

                    if (((hitHighEdge.GetX() > xOnLine) && (hitLowEdge.GetX() > xOnLine)) ||
                        ((hitHighEdge.GetX() < xOnLine) && (hitLowEdge.GetX() < xOnLine)))
                    {
                    
                        if (!((highEdgeDistanceFromLine < m_closestDistanceToLine) || (lowEdgeDistanceFromLine < m_closestDistanceToLine)))
                            continue;
                    }

                    ++hitsCollected;
                    hitPositionVector.push_back(hitPosition);
                    clusterToCaloHitListMap[pCluster].push_back(pCaloHit);
                }
            }
        }
        catch (const StatusCodeException &)
        {
            return;
        }

        ++count;
    }

    if (!clusterToCaloHitListMap.empty())
    {
        std::sort(hitPositionVector.begin(), hitPositionVector.end(), LArClusterHelper::SortCoordinatesByPosition);
        clusterAssociation.IsEndUpstream() ? clusterAssociation.SetUpstreamMergePoint(hitPositionVector.front()) : clusterAssociation.SetDownstreamMergePoint(hitPositionVector.back());
    }
    else
    {
        clusterAssociation.IsEndUpstream() ? clusterAssociation.SetUpstreamMergePoint(extrapolatedEndPosition) : clusterAssociation.SetDownstreamMergePoint(extrapolatedEndPosition);
    }

    //////////////////////////
    
    std::cout << "VIEW IN EXTRAPOLATED HIT FUNCTION" << std::endl;
    ClusterList mainCluster({clusterAssociation.GetMainTrackCluster()});
    const CartesianVector &upstream(clusterAssociation.GetUpstreamMergePoint());
    const CartesianVector &downstream(clusterAssociation.GetDownstreamMergePoint());

    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &upstream, "UPSTREAM", BLACK, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &downstream, "DOWNSTREAM", BLACK, 2);
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &mainCluster, "CLUSTER", BLACK);
    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    
    //////////////////////////    
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayEndpointCorrectionAlgorithm::CreateMainTrack(ClusterEndpointAssociation &clusterEndpointAssociation, const ClusterToCaloHitListMap &clusterToCaloHitListMap,
    const ClusterList *const pClusterList, ClusterVector &clusterVector, SlidingFitResultMapPair &slidingFitResultMapPair, ClusterAssociationVector &clusterAssociationVector) const
{
    const Cluster *pMainTrackCluster(clusterEndpointAssociation.GetMainTrackCluster());
    const CartesianVector &clusterMergePoint(clusterEndpointAssociation.IsEndUpstream() ?
        clusterEndpointAssociation.GetDownstreamMergePoint() : clusterEndpointAssociation.GetUpstreamMergePoint());

    ////////////////
    //ClusterList originalTrack({pMainTrackCluster});
    //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &originalTrack, "ORIGINAL TRACK", BLACK);
    ///////////////
    
    // Determine the shower clusters which contain hits that belong to the main track
    ClusterVector showerClustersToFragment;
    for (auto &entry : clusterToCaloHitListMap)
    {
        if (entry.first != pMainTrackCluster)
            showerClustersToFragment.push_back(entry.first);
    }

    std::sort(showerClustersToFragment.begin(), showerClustersToFragment.end(), LArClusterHelper::SortByNHits);

    ////////////////
    /*
    for (auto &entry : showerClustersToFragment)
    {
        ClusterList jam({entry});
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &jam, "SHOWER CLUSTER", VIOLET);
    }
    */
    ////////////////    

    ClusterList remnantClusterList;
    pMainTrackCluster = RemoveOffAxisHitsFromTrack(pMainTrackCluster, clusterMergePoint, clusterEndpointAssociation.IsEndUpstream(),
        clusterToCaloHitListMap, remnantClusterList, *slidingFitResultMapPair.first, *slidingFitResultMapPair.second);
   
    ////////////////
    /*
    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    ClusterList refinedCluster({pMainTrackCluster});
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &refinedCluster, "REFINED MAIN TRACK", BLACK);
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &remnantClusterList, "REMNANT CLUSTERS", VIOLET);
    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    */
    ////////////////   
   
    for (const Cluster *const pShowerCluster : showerClustersToFragment)
    {
        const CaloHitList &caloHitsToMerge(clusterToCaloHitListMap.at(pShowerCluster));
        this->AddHitsToMainTrack(pMainTrackCluster, pShowerCluster, caloHitsToMerge, clusterEndpointAssociation, remnantClusterList);
    }

    ////////////////
    /*
    ClusterList extendedCluster({pMainTrackCluster});
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &extendedCluster, "REFINED MAIN TRACK", BLACK);
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &remnantClusterList, "REMNANT CLUSTERS", VIOLET);
    PandoraMonitoringApi::ViewEvent(this->GetPandora());    
    */
    ////////////////   
    
    ClusterList createdClusters;
    this->ProcessRemnantClusters(remnantClusterList, pMainTrackCluster, pClusterList, createdClusters);

    ////////////////
    /*
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &createdClusters, "CREATED CLUSTERS", RED);
    ClusterList extendedCluster({pMainTrackCluster});
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &extendedCluster, "REFINED MAIN TRACK", BLACK);
    PandoraMonitoringApi::ViewEvent(this->GetPandora());     
    */
    ////////////////       

    //ATTN: Update references to pMainTrackCluster
    //      IT IS REALLY IMPORTANT THAT THIS IS DONE FIRST I.E. DELETE NULL POINTERS FROM MAP BEFORE ANY ARE ADDED IN AS THEY CAN BECOME THE SAME
    this->UpdateAfterMainTrackModification(pMainTrackCluster, clusterEndpointAssociation, slidingFitResultMapPair, clusterAssociationVector);
    
    // ATTN: Update containers and the cluster association
    ClusterList modifiedClusters(showerClustersToFragment.begin(), showerClustersToFragment.end());   
    this->UpdateContainers(createdClusters, modifiedClusters, clusterVector, slidingFitResultMapPair);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayEndpointCorrectionAlgorithm::UpdateAfterMainTrackModification(const Cluster *const pMainTrackCluster, ClusterEndpointAssociation &clusterEndpointAssociation,
    SlidingFitResultMapPair &slidingFitResultMapPair, ClusterAssociationVector &clusterAssociationVector) const
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

    for (ClusterEndpointAssociation &clusterAssociation : clusterAssociationVector)
    {
        if (clusterAssociation.GetMainTrackCluster() == pDeletedTrackCluster)
            clusterAssociation.SetMainTrackCluster(pMainTrackCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayEndpointCorrectionAlgorithm::RemoveClusterAssociationFromClusterVector(const ClusterEndpointAssociation &clusterAssociation, ClusterVector &clusterVector) const
{
    const ClusterVector::const_iterator clusterToDeleteIter(std::find(clusterVector.begin(), clusterVector.end(), clusterAssociation.GetMainTrackCluster()));
    if (clusterToDeleteIter != clusterVector.end())
        clusterVector.erase(clusterToDeleteIter);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CosmicRayEndpointCorrectionAlgorithm::IsExtrapolatedEndpointNearBoundary(const ClusterEndpointAssociation &clusterAssociation, const float boundaryTolerance) const
{
    const CartesianVector extrapolatedPoint(clusterAssociation.IsEndUpstream() ? clusterAssociation.GetUpstreamMergePoint() : clusterAssociation.GetDownstreamMergePoint());
    
    const LArTPC *pLArTPC(&this->GetPandora().GetGeometry()->GetLArTPC());
    const bool isNearestBoundaryHigherX((extrapolatedPoint.GetX() - pLArTPC->GetCenterX()) > 0.f);
    const float nearestTPCBoundaryX(pLArTPC->GetCenterX() + (pLArTPC->GetWidthX() * 0.5f * (isNearestBoundaryHigherX ? 1.f : -1.f)));

    const float distanceFromTPCBoundary(std::fabs(extrapolatedPoint.GetX() - nearestTPCBoundaryX));

    return (distanceFromTPCBoundary < boundaryTolerance);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayEndpointCorrectionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCaloHits", m_minCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxDistanceFromTPC", m_maxDistanceFromTPC));

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

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DeltaRayslidingFitWindow", m_deltaRayslidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "GrowingFitInitialLength", m_growingFitInitialLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "GrowingFitSegmentLength", m_growingFitSegmentLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FurthestDistanceToLine", m_furthestDistanceToLine));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClosestDistanceToLine", m_closestDistanceToLine));    

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





/*

        CartesianVector microPosition(0.f, 0.f, 0.f);
        subsetFit.GetGlobalPosition(microIter->second.GetL(), microIter->second.GetFitT(), microPosition);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &microPosition, "MICRO POSITION", BLACK, 2);
        std::cout << "ANGLE DEVIATION: " << microOpeningAngle << std::endl;
        std::cout << "<-------------------- CHANGE FROM LAST LAYER: " << layerAngleDeviation << std::endl;
*/
