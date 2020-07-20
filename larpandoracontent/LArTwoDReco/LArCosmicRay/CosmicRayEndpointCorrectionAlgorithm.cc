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
    
CosmicRayEndpointCorrectionAlgorithm::CosmicRayEndpointCorrectionAlgorithm() :
    m_minCaloHits(50),
    m_maxDistanceFromTPC(15.f),
    m_curveThreshold(0.3f),
    m_minScaledZOffset(0.25f),
    m_thresholdAngleDeviation(10.f),
    m_thresholdAngleDeviationBetweenLayers(1.f),
    m_maxSmoothPoints(2),
    m_thresholdMaxAngleDeviation(25.f)
{
}

StatusCode CosmicRayEndpointCorrectionAlgorithm::Run()
{
    PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_DEFAULT, -1.f, 1.f, 1.f);
    
    const ClusterList *pClusterList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));
    
    ClusterVector clusterVector;
    this->SelectCleanClusters(pClusterList, clusterVector);

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);
    
    TwoDSlidingFitResultMap microSlidingFitResultMap, macroSlidingFitResultMap;
    this->InitialiseSlidingFitResultMaps(clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap);

    while (!clusterVector.empty())
    {
        const Cluster *pCluster(clusterVector.front());
        
        TwoDSlidingFitResultMap::const_iterator microFitResult(microSlidingFitResultMap.find(pCluster));
        TwoDSlidingFitResultMap::const_iterator macroFitResult(macroSlidingFitResultMap.find(pCluster));

        if (microFitResult == microSlidingFitResultMap.end() || macroFitResult == macroSlidingFitResultMap.end())
            continue;        

        // EXAMINE UPSTREAM POINT OF CLUSTER

        // EXAMINE DOWNSTREAM POINT OF CLUSTER
    }
        

    
    return STATUS_CODE_SUCCESS;
}    
    
//------------------------------------------------------------------------------------------------------------------------------------------    

    
void CosmicRayEndpointCorrectionAlgorithm::ExamineClusterEndpoint(const Cluster *const pCluster, const bool isUpstream, const TwoDSlidingFitResult &microFitResult, const TwoDSlidingFitResult &macroFitResult)
{
    const CartesianVector &endpointPosition(isUpstream ? microFitResult.GetGlobalMaxLayerPosition() : microFitResult.GetGlobalMinLayerPosition());

    //ClusterList theCluster({pCluster});
    //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &theCluster, "CLUSTER", BLACK);
    //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &upstreamPosition, "UPSTREAM", RED, 2);
    //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &downstreamPosition, "DOWNSTREAM", RED, 2);

    CartesianVector clusterMergePoint(0.f, 0.f, 0.f), clusterMergeDirection(0.f, 0.f, 0.f);
    if (!GetClusterMergingCoordinates(microFitResult, macroFitResult, macroFitResult, isUpstream, clusterMergePoint, clusterMergeDirection))
    {
        std::cout << "CANNOT FIND MERGE POSITION" << std::endl;
        PandoraMonitoringApi::ViewEvent(this->GetPandora());        
        return;
    }

    // Check if needs endpoint removal - meetCriteria?
    if (!this->NewIsCosmicRay(endpointPosition, clusterMergePoint, clusterMergeDirection, pCluster, isUpstream))
    {
        std::cout << "DOESN'T SATISY CRITERIA" << std::endl;
        PandoraMonitoringApi::ViewEvent(this->GetPandora());              
        return;
    }

    // REMOVE HITS - CREATE AN ACTUAL CLUSTER ASSOCIATION HERE
    std::cout << "REMOVING HITS..." << std::endl;

    const float predictedGradient(downstreamMergeDirection.GetZ() / downstreamMergeDirection.GetX());
    const float predictedIntercept(downstreamMergePoint.GetZ() - (predictedGradient * downstreamMergePoint.GetX()));
    const CartesianVector extrapolatedPoint(endpointPosition.GetX(), 0.f, predictedIntercept + (predictedGradient * endpointPosition.GetX()));

    ClusterAssociation clusterAssociation;
    if (isUpstream)
    {
        // ISOBEL : CHECK NEGATIVES
        clusterAssociation = ClusterAssociation(pCluster, pCluster, extrapolatedPoint, clusterMergeDirection, clusterMergePoint, (-1) * clusterMergeDirection);
    }
    else
    {
        // ISOBEL : CHECK NEGATIVES
        clusterAssociation = ClusterAssociation(pCluster, pCluster, clusterMergePoint, (-1) * clusterMergeDirection, extrpolatedPoint, clusterMergeDirection);    

    }

    std::cout << "isUpstream: " << isUpstream << std::endl;

    /*
    GetExtrapolatedCaloHits(downstreamMergePoint, extrapolatedPoint, downstreamMergeDirection, pClusterList, downstreamClusterToCaloHitListMap);

            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &downstreamMergePoint, "DOWNSTREAM", BLACK, 2);
            
            // VISUALIZE COLLECTED CALO HITS
            for (auto entry : downstreamClusterToCaloHitListMap)
            {
                const CaloHitList &caloHitList(entry.second);

                for (auto &hit :caloHitList)
                {
                    const CartesianVector &hitPosition(hit->GetPositionVector());
                    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "POSITION", GREEN, 2);
                }
            }
            
            PandoraMonitoringApi::ViewEvent(this->GetPandora());

            ClusterList remnantClusterList;
            if (!downstreamClusterToCaloHitListMap.empty())
            {
                pCluster = RemoveOffAxisHitsFromTrack(pCluster, downstreamMergePoint, true, downstreamClusterToCaloHitListMap,
                    remnantClusterList, microSlidingFitResultMap, macroSlidingFitResultMap, clusterVector);

                ClusterList theNe({pCluster});
                PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &theNe, "NEW CLUSTER", RED);
                PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &remnantClusterList, "REMNANT", BLUE);
                PandoraMonitoringApi::ViewEvent(this->GetPandora());
            }

        }
        else
        {
            std::cout << "WILL NOT REMOVE" << std::endl;
        }

        ClusterToCaloHitListMap upstreamClusterToCaloHitListMap;        
        if (this->NewIsCosmicRay(upstreamPosition, upstreamMergePoint, upstreamMergeDirection, pCluster, false))
        {
            std::cout << "WILL REMOVE" << std::endl;
            PandoraMonitoringApi::ViewEvent(this->GetPandora());

            // REMOVE HITS
            std::cout << "REMOVING HITS..." << std::endl;

            const float predictedGradient(upstreamMergeDirection.GetZ() / upstreamMergeDirection.GetX());
            const float predictedIntercept(upstreamMergePoint.GetZ() - (predictedGradient * upstreamMergePoint.GetX()));
            const CartesianVector extrapolatedPoint(upstreamPosition.GetX(), 0.f, predictedIntercept + (predictedGradient * upstreamPosition.GetX()));
            
            GetExtrapolatedCaloHits(extrapolatedPoint, upstreamMergePoint, upstreamMergeDirection, pClusterList, upstreamClusterToCaloHitListMap);

            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &upstreamMergePoint, "UPSTREAM", BLACK, 2);
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &upstreamMergePoint, "UPSTREAM", BLACK, 2);

            // VISUALIZE COLLECTED CALO HITS
            for (auto entry : upstreamClusterToCaloHitListMap)
            {
                const CaloHitList &caloHitList(entry.second);

                for (auto &hit :caloHitList)
                {
                    const CartesianVector &hitPosition(hit->GetPositionVector());
                    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "POSITION", GREEN, 2);
                }
            }

            PandoraMonitoringApi::ViewEvent(this->GetPandora());


            ClusterList remnantClusterList;
            if (!upstreamClusterToCaloHitListMap.empty())
            {
                const Cluster *const newCluster(RemoveOffAxisHitsFromTrack(pCluster, upstreamMergePoint, false, upstreamClusterToCaloHitListMap,
                    remnantClusterList, microSlidingFitResultMap, macroSlidingFitResultMap, clusterVector));

                ClusterList theNe({newCluster});
                PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &theNe, "NEW CLUSTER", RED);
                PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &remnantClusterList, "REMNANT", BLUE);
                PandoraMonitoringApi::ViewEvent(this->GetPandora());
            }

        }
        else
        {
            std::cout << "WILL NOT REMOVE" << std::endl;
        }            

        UpdateForClusterDeletion(pCluster, clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap);
        
        // search for region along theorised CR that has the largest angle deviation from the mergin point direction
        // if found, remove hits??
        
    }
    */
}

    
//------------------------------------------------------------------------------------------------------------------------------------------ 
    /*
StatusCode CosmicRayEndpointCorrectionAlgorithm::Run()
{
    PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_DEFAULT, -1.f, 1.f, 1.f);
    
    const ClusterList *pClusterList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));
    
    ClusterVector clusterVector;
    this->SelectCleanClusters(pClusterList, clusterVector);

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);
    
    TwoDSlidingFitResultMap microSlidingFitResultMap, macroSlidingFitResultMap;
    this->InitialiseSlidingFitResultMaps(clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap);

    while (!clusterVector.empty())
    {
        const Cluster *pCluster(clusterVector.front());
        
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


        ClusterToCaloHitListMap downstreamClusterToCaloHitListMap;
        if (this->NewIsCosmicRay(downstreamPosition, downstreamMergePoint, downstreamMergeDirection, pCluster, true))
        {
            std::cout << "WILL REMOVE" << std::endl;
            PandoraMonitoringApi::ViewEvent(this->GetPandora());

            // REMOVE HITS - CREATE AN ACTUAL CLUSTER ASSOCIATION HERE
            std::cout << "REMOVING HITS..." << std::endl;

            const float predictedGradient(downstreamMergeDirection.GetZ() / downstreamMergeDirection.GetX());
            const float predictedIntercept(downstreamMergePoint.GetZ() - (predictedGradient * downstreamMergePoint.GetX()));
            const CartesianVector extrapolatedPoint(downstreamPosition.GetX(), 0.f, predictedIntercept + (predictedGradient * downstreamPosition.GetX()));

            GetExtrapolatedCaloHits(downstreamMergePoint, extrapolatedPoint, downstreamMergeDirection, pClusterList, downstreamClusterToCaloHitListMap);

            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &downstreamMergePoint, "DOWNSTREAM", BLACK, 2);
            
            // VISUALIZE COLLECTED CALO HITS
            for (auto entry : downstreamClusterToCaloHitListMap)
            {
                const CaloHitList &caloHitList(entry.second);

                for (auto &hit :caloHitList)
                {
                    const CartesianVector &hitPosition(hit->GetPositionVector());
                    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "POSITION", GREEN, 2);
                }
            }
            
            PandoraMonitoringApi::ViewEvent(this->GetPandora());

            ClusterList remnantClusterList;
            if (!downstreamClusterToCaloHitListMap.empty())
            {
                pCluster = RemoveOffAxisHitsFromTrack(pCluster, downstreamMergePoint, true, downstreamClusterToCaloHitListMap,
                    remnantClusterList, microSlidingFitResultMap, macroSlidingFitResultMap, clusterVector);

                ClusterList theNe({pCluster});
                PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &theNe, "NEW CLUSTER", RED);
                PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &remnantClusterList, "REMNANT", BLUE);
                PandoraMonitoringApi::ViewEvent(this->GetPandora());
            }

        }
        else
        {
            std::cout << "WILL NOT REMOVE" << std::endl;
        }

        ClusterToCaloHitListMap upstreamClusterToCaloHitListMap;        
        if (this->NewIsCosmicRay(upstreamPosition, upstreamMergePoint, upstreamMergeDirection, pCluster, false))
        {
            std::cout << "WILL REMOVE" << std::endl;
            PandoraMonitoringApi::ViewEvent(this->GetPandora());

            // REMOVE HITS
            std::cout << "REMOVING HITS..." << std::endl;

            const float predictedGradient(upstreamMergeDirection.GetZ() / upstreamMergeDirection.GetX());
            const float predictedIntercept(upstreamMergePoint.GetZ() - (predictedGradient * upstreamMergePoint.GetX()));
            const CartesianVector extrapolatedPoint(upstreamPosition.GetX(), 0.f, predictedIntercept + (predictedGradient * upstreamPosition.GetX()));
            
            GetExtrapolatedCaloHits(extrapolatedPoint, upstreamMergePoint, upstreamMergeDirection, pClusterList, upstreamClusterToCaloHitListMap);

            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &upstreamMergePoint, "UPSTREAM", BLACK, 2);
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &upstreamMergePoint, "UPSTREAM", BLACK, 2);

            // VISUALIZE COLLECTED CALO HITS
            for (auto entry : upstreamClusterToCaloHitListMap)
            {
                const CaloHitList &caloHitList(entry.second);

                for (auto &hit :caloHitList)
                {
                    const CartesianVector &hitPosition(hit->GetPositionVector());
                    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "POSITION", GREEN, 2);
                }
            }

            PandoraMonitoringApi::ViewEvent(this->GetPandora());


            ClusterList remnantClusterList;
            if (!upstreamClusterToCaloHitListMap.empty())
            {
                const Cluster *const newCluster(RemoveOffAxisHitsFromTrack(pCluster, upstreamMergePoint, false, upstreamClusterToCaloHitListMap,
                    remnantClusterList, microSlidingFitResultMap, macroSlidingFitResultMap, clusterVector));

                ClusterList theNe({newCluster});
                PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &theNe, "NEW CLUSTER", RED);
                PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &remnantClusterList, "REMNANT", BLUE);
                PandoraMonitoringApi::ViewEvent(this->GetPandora());
            }

        }
        else
        {
            std::cout << "WILL NOT REMOVE" << std::endl;
        }            

        UpdateForClusterDeletion(pCluster, clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap);
        
        // search for region along theorised CR that has the largest angle deviation from the mergin point direction
        // if found, remove hits??
        
    }
    
    return STATUS_CODE_SUCCESS;
}    
    */
//------------------------------------------------------------------------------------------------------------------------------------------    

void CosmicRayEndpointCorrectionAlgorithm::SelectCleanClusters(const ClusterList *pClusterList, ClusterVector &clusterVector) const
{
    for (const Cluster *const pCluster : *pClusterList)
    {
        if (pCluster->GetNCaloHits() < m_minCaloHits)
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
    const float allowanceHighXEdge(tpcHighXEdge - m_maxDistanceFromTPC), allowanceLowXEdge(tpcLowXEdge + m_maxDistanceFromTPC);

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
    
    if (averageDistanceFromLine < m_curveThreshold)
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

    if (scaledDeltaZ < m_minScaledZOffset)
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
    unsigned int tooSmooth(0);
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
        
        if (std::fabs(microOpeningAngle > m_thresholdAngleDeviation) && !reachedFirstCurve)
        {
            reachedFirstCurve = true;
            isClockwise = (microOpeningAngle < 0.f);
            previousOpeningAngle = microOpeningAngle;
            continue;
        }

        if (reachedFirstCurve)
        {
            if ((isClockwise && (microOpeningAngle > previousOpeningAngle)) || (!isClockwise && (microOpeningAngle < previousOpeningAngle)) ||
                (layerAngleDeviation < m_thresholdAngleDeviationBetweenLayers))
            {
                ++tooSmooth;

                if (tooSmooth > m_maxSmoothPoints)
                {
                    std::cout << "TOO SMOOTH" << std::endl;
                    PandoraMonitoringApi::ViewEvent(this->GetPandora());
                    return false;
                }
            }
            else
            {
                tooSmooth = 0;
            
                if (std::fabs(microOpeningAngle) > m_thresholdMaxAngleDeviation)
                {
                    std::cout << "MEET ANGLE CRITERIA" << std::endl;
                    metCriteria = true;
                    break;
                }
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

//------------------------------------------------------------------------------------------------------------------------------------------    

StatusCode CosmicRayEndpointCorrectionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCaloHits", m_minCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxDistanceFromTPC", m_maxDistanceFromTPC));

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
        "MaxSmoothPoints", m_maxSmoothPoints));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ThresholdMaxAngleDeviation", m_thresholdMaxAngleDeviation));           

        
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content


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
