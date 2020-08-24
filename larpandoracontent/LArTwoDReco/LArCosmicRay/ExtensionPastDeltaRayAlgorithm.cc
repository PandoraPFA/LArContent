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
            
        CartesianVector clusterMergePoint(0.f, 0.f, 0.f), clusterMergeDirection(0.f, 0.f, 0.f);
        if (!GetClusterMergingCoordinates(microSlidingFitResult, macroSlidingFitResult, macroSlidingFitResult, isEndUpstream, clusterMergePoint, clusterMergeDirection))
            continue;

        // Reject clusters that do not cross TPC boundary
        if (std::fabs(clusterMergeDirection.GetX()) < std::numeric_limits<float>::epsilon())
            continue;
            
        const CartesianVector &endpointPosition(isEndUpstream ? microSlidingFitResult.GetGlobalMinLayerPosition() : microSlidingFitResult.GetGlobalMaxLayerPosition());
        const float endpointSeparation((endpointPosition - clusterMergePoint).GetMagnitude());

        // Reject if no clustering error
        if (endpointSeparation < std::numeric_limits<float>::epsilon())
            continue;

        // Reject if not close enough to the TPC boundary
        const bool isMergePointInAllowanceRegion(std::fabs(nearestTPCBoundaryX - clusterMergePoint.GetX()) < m_maxDistanceFromTPC);
        const bool isEndpointInAllowanceRegion(std::fabs(nearestTPCBoundaryX - endpointPosition.GetX()) < m_maxDistanceFromTPC);

        if (!(isMergePointInAllowanceRegion || isEndpointInAllowanceRegion))
            continue;
            
        const float predictedGradient(clusterMergeDirection.GetZ() / clusterMergeDirection.GetX());
        const float predictedIntercept(clusterMergePoint.GetZ() - (predictedGradient * clusterMergePoint.GetX()));
        const CartesianVector fitEndpointPosition(endpointPosition.GetX(), 0.f, predictedIntercept + (predictedGradient * endpointPosition.GetX()));
        const float deltaZ(std::fabs(endpointPosition.GetZ() - fitEndpointPosition.GetZ()));

        // Reject if no significant endpoint deviation
        if (deltaZ < m_minZOffset)
            continue;

        // Reject if cannot find a delta ray bend
        if(!this->IsDeltaRay(microSlidingFitResult, clusterMergeDirection, isEndUpstream, clusterMergePoint))
            continue;

        // ATTN: Temporarily set the other merge point to define extrapolate hits search region        
        const CartesianVector extrapolatedHitsEndpoint(nearestTPCBoundaryX, 0.f, predictedIntercept + (predictedGradient * nearestTPCBoundaryX));

        if (isEndUpstream ? clusterMergePoint.GetZ() < extrapolatedHitsEndpoint.GetZ() : clusterMergePoint.GetZ() > extrapolatedHitsEndpoint.GetZ())
            continue;    

        clusterAssociation = isEndUpstream ?
            ClusterEndpointAssociation(extrapolatedHitsEndpoint, clusterMergeDirection, clusterMergePoint, clusterMergeDirection * (-1.f), pCluster, true) :
            ClusterEndpointAssociation(clusterMergePoint, clusterMergeDirection, extrapolatedHitsEndpoint, clusterMergeDirection * (-1.f), pCluster, false);          
        
        return true;
    }

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
                return false;

            if(this->IsCurvePresent(subsetFit, isEndUpstream, false, clusterMergeDirection, clusterMergePoint))
                return true;
        }
    }
    catch (const StatusCodeException &) {}

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
            
        if (!isBendConsistent || (microOpeningAngle * previousOpeningAngle < 0.f) || (layerAngleDeviation < m_thresholdAngleDeviationBetweenLayers))
        {
            ++anomalousLayerCount;

            if (anomalousLayerCount > m_maxAnomalousPoints)
                break;
        }
        else
        {
            if (isClusterMergePointEnd)
            {
                if (std::fabs(microOpeningAngle) > m_thresholdMaxAngleDeviation)
                {
                    clusterMergePoint = isEndUpstream ? subsetFit.GetGlobalMaxLayerPosition() : subsetFit.GetGlobalMinLayerPosition();
                    return true;
                }
            }
            else
            {
               if (std::fabs(microOpeningAngle) < m_thresholdAngleDeviation)
               {
                   clusterMergePoint = isEndUpstream ? subsetFit.GetGlobalMaxLayerPosition() : subsetFit.GetGlobalMinLayerPosition();
                   return true;
               }
            }

            anomalousLayerCount = 0;
        }

        previousOpeningAngle = microOpeningAngle;
    }

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
