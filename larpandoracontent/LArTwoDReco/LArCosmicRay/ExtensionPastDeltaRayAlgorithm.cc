/**
 *  @file   larpandoracontent/LArTwoDReco/LArCosmicRay/ExtensionPastDeltaRayAlgorithm.cc 
 *
 *  @brief  Implementation of the cosmic ray endpoint correction class
 *
 *  $Log: $
 */
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArTwoDReco/LArCosmicRay/ExtensionPastDeltaRayAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

using namespace pandora;

namespace lar_content
{
    
ExtensionPastDeltaRayAlgorithm::ExtensionPastDeltaRayAlgorithm() :
    m_maxDistanceFromTPC(15.f),
    m_minZOffset(2.f),
    m_deltaRayslidingFitWindow(5),    
    m_thresholdAngleDeviation(10.f),
    m_thresholdMaxAngleDeviation(21.f),     
    m_thresholdAngleDeviationBetweenLayers(1.f),
    m_maxAnomalousPoints(2)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------    
    
bool ExtensionPastDeltaRayAlgorithm::FindBestClusterAssociation(const ClusterVector &clusterVector, const SlidingFitResultMapPair &slidingFitResultMapPair,
    const ClusterList *const /*pClusterList*/, const bool isHigherXBoundary, ClusterEndpointAssociation &clusterAssociation)
{
    const float nearestTPCBoundaryX(isHigherXBoundary ? m_tpcMaxXEdge : m_tpcMinXEdge);
    
    //ATTN: Find best cluster association - method assumes clusters ordered by furthest distance from TPC
    for (const Cluster *const pCluster : clusterVector)
    {
        const TwoDSlidingFitResult &microSlidingFitResult(slidingFitResultMapPair.first->at(pCluster));
        const TwoDSlidingFitResult &macroSlidingFitResult(slidingFitResultMapPair.second->at(pCluster));


        const bool isEndUpstream = (std::fabs(microSlidingFitResult.GetGlobalMinLayerPosition().GetX() - nearestTPCBoundaryX) <
                                    std::fabs(microSlidingFitResult.GetGlobalMaxLayerPosition().GetX() - nearestTPCBoundaryX));

        ////////////////
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

        // Reject clusters that do not cross TPC boundary
        if (std::fabs(clusterMergeDirection.GetX()) < std::numeric_limits<float>::epsilon())
        {
            std::cout << "MERGE DIRECTION HAS NO X COMPONENT" << std::endl;
            PandoraMonitoringApi::ViewEvent(this->GetPandora());                 
            continue;
        }
            
        const CartesianVector &endpointPosition(isEndUpstream ? microSlidingFitResult.GetGlobalMinLayerPosition() : microSlidingFitResult.GetGlobalMaxLayerPosition());
        const float endpointSeparation((endpointPosition - clusterMergePoint).GetMagnitude());

        /////////////////
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &endpointPosition, "ENDPOINT", BLUE, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &clusterMergePoint, "MERGE POINT", BLUE, 2);
        std::cout << "Endpoint Separation: " << endpointSeparation << std::endl;
        /////////////////            

        // Reject if no clustering error
        if (endpointSeparation < std::numeric_limits<float>::epsilon())
        {
            std::cout << "MERGE POINT AND ENDPOINT ARE THE SAME" << std::endl;
            PandoraMonitoringApi::ViewEvent(this->GetPandora());            
            continue;
        }

        // Reject if not close enough to the TPC boundary
        const bool isMergePointInAllowanceRegion(std::fabs(nearestTPCBoundaryX - clusterMergePoint.GetX()) < m_maxDistanceFromTPC);
        const bool isEndpointInAllowanceRegion(std::fabs(nearestTPCBoundaryX - endpointPosition.GetX()) < m_maxDistanceFromTPC);

        if (!(isMergePointInAllowanceRegion || isEndpointInAllowanceRegion))
        {
            std::cout << "NOT CLOSE ENOUGH TO THE BOUNDARY" << std::endl;
            PandoraMonitoringApi::ViewEvent(this->GetPandora());            
            continue;
        }
            
        const float predictedGradient(clusterMergeDirection.GetZ() / clusterMergeDirection.GetX());
        const float predictedIntercept(clusterMergePoint.GetZ() - (predictedGradient * clusterMergePoint.GetX()));
        const CartesianVector fitEndpointPosition(endpointPosition.GetX(), 0.f, predictedIntercept + (predictedGradient * endpointPosition.GetX()));
        const float deltaZ(std::fabs(endpointPosition.GetZ() - fitEndpointPosition.GetZ()));

        ///////////////////
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &fitEndpointPosition, "FIT ENDPOINT", ORANGE, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &clusterMergePoint, "clusterMergePoint", VIOLET, 2);
        const CartesianVector start(clusterMergePoint + (clusterMergeDirection*40));
        const CartesianVector end(clusterMergePoint - (clusterMergeDirection*40));
        PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &start, &end, "CLUSTER DIRECTION", GREEN, 4, 2);            
        std::cout << "deltaZ: " << deltaZ << std::endl;
        /////////////////        

        // Reject if no significant endpoint deviation
        if (deltaZ < m_minZOffset)
        {
            std::cout << "DELTA Z IS TOO LOW" << std::endl;
            PandoraMonitoringApi::ViewEvent(this->GetPandora());            
            continue;
        }

        // Reject if cannot find a delta ray bend
        if(!this->IsDeltaRay(microSlidingFitResult, clusterMergeDirection, isEndUpstream, clusterMergePoint))
        {
            std::cout << "NOT A DELTA RAY" << std::endl;
            PandoraMonitoringApi::ViewEvent(this->GetPandora());              
            continue;
        }

        // ATTN: Temporarily set the other merge point to define extrapolate hits search region        
        const CartesianVector extrapolatedHitsEndpoint(nearestTPCBoundaryX, 0.f, predictedIntercept + (predictedGradient * nearestTPCBoundaryX));

        if (isEndUpstream ? clusterMergePoint.GetZ() < extrapolatedHitsEndpoint.GetZ() : clusterMergePoint.GetZ() > extrapolatedHitsEndpoint.GetZ())
        {
            std::cout << "EXTRAPOLATED ENDPOINT IS NOT IN FORWARD DIRECTION" << std::endl;
            PandoraMonitoringApi::ViewEvent(this->GetPandora());            
            continue;
        }

        clusterAssociation = isEndUpstream ?
            ClusterEndpointAssociation(extrapolatedHitsEndpoint, clusterMergeDirection, clusterMergePoint, clusterMergeDirection * (-1.f), pCluster, true) :
            ClusterEndpointAssociation(clusterMergePoint, clusterMergeDirection, extrapolatedHitsEndpoint, clusterMergeDirection * (-1.f), pCluster, false);          

       ////////////////////////////////
        ClusterList mainCluster({clusterAssociation.GetMainTrackCluster()});
        const CartesianVector &upstream(clusterAssociation.GetUpstreamMergePoint());
        const CartesianVector &downstream(clusterAssociation.GetDownstreamMergePoint());
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &upstream, "UPSTREAM", VIOLET, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &downstream, "DOWNSTREAM", RED, 2);
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &mainCluster, "CLUSTER", BLACK);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        ////////////////////////////////   
        
        return true;
    }

    std::cout << "DID NOT FIND AN ASSOCIATION" << std::endl;    
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ExtensionPastDeltaRayAlgorithm::IsDeltaRay(const TwoDSlidingFitResult &microSlidingFitResult, const CartesianVector &clusterMergeDirection, const bool isEndUpstream,
    CartesianVector &clusterMergePoint) const 
{
    // Perform more detailed fit of hits past the clusterMergePoint
    float rL(0.f), rT(0.f);
    microSlidingFitResult.GetLocalPosition(clusterMergePoint, rL, rT);

    CartesianPointVector hitSubset;
    const OrderedCaloHitList &orderedCaloHitList(microSlidingFitResult.GetCluster()->GetOrderedCaloHitList());
    for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
    {
        for (const CaloHit *const pCaloHit : *mapEntry.second) 
        {
            float thisL(0.f), thisT(0.f);
            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            microSlidingFitResult.GetLocalPosition(pCaloHit->GetPositionVector(), thisL, thisT);

            if ((isEndUpstream && (thisL < rL)) || (!isEndUpstream && (thisL > rL)))
                hitSubset.push_back(hitPosition);
        }
    }
    /*
    for (const CartesianVector jam : hitSubset)
    {
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &jam, "DR point", BLACK, 2);
    }
    */
    try
    {
        const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const TwoDSlidingFitResult subsetFit(&hitSubset, m_deltaRayslidingFitWindow, slidingFitPitch);

        // Check for curve at the clusterMergePoint end
        if(this->IsCurvePresent(subsetFit, isEndUpstream, true, clusterMergeDirection, clusterMergePoint))
        {
            return true;
        }
        else
        {
            // Check for curve at the endpoint end
            const CartesianVector &endDirection(isEndUpstream ? subsetFit.GetGlobalMinLayerDirection() : subsetFit.GetGlobalMaxLayerDirection());
            const float endOpeningAngle(endDirection.GetOpeningAngle(clusterMergeDirection) * 180 / 3.14);
            
            if (std::fabs(endOpeningAngle) < m_thresholdMaxAngleDeviation)
            {
                std::cout << "endpoint opnening angle: " << std::to_string(std::fabs(endOpeningAngle)) << "not large enough" << std::endl;
                return false;
            }

            if(this->IsCurvePresent(subsetFit, isEndUpstream, false, clusterMergeDirection, clusterMergePoint))
                return true;
        }
    }
    catch (const StatusCodeException &)
    {
        std::cout << "CANNOT MAKE A FIT" << std::endl;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ExtensionPastDeltaRayAlgorithm::IsCurvePresent(const TwoDSlidingFitResult &subsetFit, const bool isEndUpstream, const bool isClusterMergePointEnd,
    const CartesianVector &clusterMergeDirection, CartesianVector &clusterMergePoint) const
{
    const LayerFitResultMap &clusterMicroLayerFitResultMap(subsetFit.GetLayerFitResultMap());
    const bool isStartLayerUpstream(isEndUpstream ? (isClusterMergePointEnd ? false : true) : (isClusterMergePointEnd ? true : false)); 
    const int startLayer(isStartLayerUpstream ? subsetFit.GetMinLayer() : subsetFit.GetMaxLayer());
    const int endLayer(isStartLayerUpstream ? subsetFit.GetMaxLayer() : subsetFit.GetMinLayer());
    const int loopTerminationLayer(isStartLayerUpstream ? endLayer + 1 : endLayer - 1);
    const int step(isStartLayerUpstream ? 1 : -1);
    unsigned int anomalousLayerCount(0);
    bool reachedFirstCurve(false);
    float previousOpeningAngle;

    std::cout << "CURVE SEARCH START" << std::endl;
    std::cout << "isClusterMergePointEnd: " << isClusterMergePointEnd << std::endl;
    
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
            if (isClusterMergePointEnd)
            {
                if (std::fabs(microOpeningAngle) > m_thresholdAngleDeviation)
                    reachedFirstCurve = true;
            }
            else
            {
               if (std::fabs(microOpeningAngle) < m_thresholdMaxAngleDeviation)
                    reachedFirstCurve = true;                
            }
            
            previousOpeningAngle = microOpeningAngle;
            continue;
        }

        const float layerAngleDeviation(std::fabs(microOpeningAngle - previousOpeningAngle));
        const bool isBendConsistent((isClusterMergePointEnd && (std::fabs(microOpeningAngle) > std::fabs(previousOpeningAngle))) ||
            (!isClusterMergePointEnd && (std::fabs(microOpeningAngle) < std::fabs(previousOpeningAngle))));        

        ///////////////////////////
        CartesianVector microPosition(0.f, 0.f, 0.f);
        subsetFit.GetGlobalPosition(microIter->second.GetL(), microIter->second.GetFitT(), microPosition);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &microPosition, std::to_string(microOpeningAngle), BLACK, 2);
        std::cout << "Opening angle: " << microOpeningAngle << std::endl;
        std::cout << "Layer angle deviation: " << layerAngleDeviation << std::endl;
        ///////////////////////////
        
            
        if (!isBendConsistent || (microOpeningAngle * previousOpeningAngle < 0.f) || (layerAngleDeviation < m_thresholdAngleDeviationBetweenLayers))
        {
            ++anomalousLayerCount;

            if (anomalousLayerCount > m_maxAnomalousPoints)
            {
                std::cout << "TOO MANY WOBBLES" << std::endl;
                break;
            }
        }
        else
        {
            if (isClusterMergePointEnd)
            {
                if (std::fabs(microOpeningAngle) > m_thresholdMaxAngleDeviation)
                {
                    clusterMergePoint = isEndUpstream ? subsetFit.GetGlobalMaxLayerPosition() : subsetFit.GetGlobalMinLayerPosition();
                    std::cout << "merge point pass" << std::endl;
                    return true;
                }
            }
            else
            {
               if (std::fabs(microOpeningAngle) < m_thresholdAngleDeviation)
               {
                   clusterMergePoint = isEndUpstream ? subsetFit.GetGlobalMaxLayerPosition() : subsetFit.GetGlobalMinLayerPosition();
                   std::cout << "end point pass" << std::endl;
                   return true;
               }
            }

            anomalousLayerCount = 0;
        }

        previousOpeningAngle = microOpeningAngle;
    }
    
    std::cout << "no curve found" << std::endl;
    return false;
    
}

//------------------------------------------------------------------------------------------------------------------------------------------ 
    
StatusCode ExtensionPastDeltaRayAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxDistanceFromTPC", m_maxDistanceFromTPC));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinZOffset", m_minZOffset));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DeltaRayslidingFitWindow", m_deltaRayslidingFitWindow));    

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ThresholdAngleDeviation", m_thresholdAngleDeviation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ThresholdMaxAngleDeviation", m_thresholdMaxAngleDeviation));    

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ThresholdAngleDeviationBetweenLayers", m_thresholdAngleDeviationBetweenLayers));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxAnomalousPoints", m_maxAnomalousPoints));

    return TrackExtensionRefinementAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
