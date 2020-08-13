/**
 *  @file   larpandoracontent/LArTwoDReco/LArCosmicRay/ExtensionPastDeltaRayAlgorithm.cc 
 *
 *  @brief  Implementation of the cosmic ray endpoint correction class
 *
 *  $Log: $
 */
#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArControlFlow/MultiPandoraApi.h"

#include "larpandoracontent/LArTwoDReco/LArCosmicRay/ExtensionPastDeltaRayAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArStitchingHelper.h"

using namespace pandora;

namespace lar_content
{
    
ExtensionPastDeltaRayAlgorithm::ExtensionPastDeltaRayAlgorithm() :
    m_maxDistanceFromTPC(15.f),
    m_minScaledZOffset(2.f), //0.25 m_curveThreshold(0.3f),
    m_thresholdAngleDeviation(10.f),
    m_thresholdAngleDeviationBetweenLayers(1.f),
    m_maxAnomalousPoints(2),
    m_thresholdMaxAngleDeviation(21.f),
    m_deltaRayslidingFitWindow(5)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------    
    
bool ExtensionPastDeltaRayAlgorithm::FindBestClusterAssociation(const ClusterVector &clusterVector, const SlidingFitResultMapPair &slidingFitResultMapPair,
    ClusterEndpointAssociation &clusterAssociation, const ClusterList *const /*pClusterList*/, const bool isHigherXBoundary)
{
    const float nearestTPCBoundaryX(isHigherXBoundary ? m_tpcMaxXEdge : m_tpcMinXEdge);
    
    //ATTN: Find best cluster association - method assumes clusters ordered by furthest distance from TPC
    for (const Cluster *const pCluster : clusterVector)
    {
        const TwoDSlidingFitResult &microSlidingFitResult(slidingFitResultMapPair.first->at(pCluster));
        const TwoDSlidingFitResult &macroSlidingFitResult(slidingFitResultMapPair.second->at(pCluster));


        const bool isEndUpstream = (std::fabs(microSlidingFitResult.GetGlobalMinLayerPosition().GetX() - nearestTPCBoundaryX) <
                                    std::fabs(microSlidingFitResult.GetGlobalMaxLayerPosition().GetX() - nearestTPCBoundaryX));

        // ATTN: Match the cluster to itself to get the cluster merging points
        const bool isClusterUpstream(!isEndUpstream);

        ////////////////
        /*
        ClusterList theCluster({pCluster});
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &theCluster, "CONSIDERED CLUSTER", BLACK);
        std::cout << "isEndUpstream: " << isEndUpstream << std::endl;
        */
        ////////////////
            
        CartesianVector clusterMergePoint(0.f, 0.f, 0.f), clusterMergeDirection(0.f, 0.f, 0.f);
        if (!GetClusterMergingCoordinates(microSlidingFitResult, macroSlidingFitResult, macroSlidingFitResult, isClusterUpstream, clusterMergePoint, clusterMergeDirection))
        {
            //std::cout << "CANNOT FIND MERGE POSITION" << std::endl;
            //PandoraMonitoringApi::ViewEvent(this->GetPandora());
            continue;
        }

        // WILL NOT CROSS TPC BOUNDARY
        if (std::fabs(clusterMergeDirection.GetX()) < std::numeric_limits<float>::epsilon())
        {
            //std::cout << "MERGE DIRECTION HAS NO X COMPONENT" << std::endl;
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

        const bool isMergePointInAllowanceRegion(std::fabs(nearestTPCBoundaryX - clusterMergePoint.GetX()) < m_maxDistanceFromTPC);
        const bool isEndpointInAllowanceRegion(std::fabs(nearestTPCBoundaryX - endpointPosition.GetX()) < m_maxDistanceFromTPC);

        if (!(isMergePointInAllowanceRegion || isEndpointInAllowanceRegion))
        {
            //std::cout << "NOT CLOSE ENOUGH TO THE BOUNDARY" << std::endl;
            //PandoraMonitoringApi::ViewEvent(this->GetPandora());
            continue;
        }
            
        const float predictedGradient(clusterMergeDirection.GetZ() / clusterMergeDirection.GetX());
        const float predictedIntercept(clusterMergePoint.GetZ() - (predictedGradient * clusterMergePoint.GetX()));
        const CartesianVector fitEndpointPosition(endpointPosition.GetX(), 0.f, predictedIntercept + (predictedGradient * endpointPosition.GetX()));
        const float deltaZ(std::fabs(endpointPosition.GetZ() - fitEndpointPosition.GetZ()));

        ///////////////////
        /*
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &fitEndpointPosition, "FIT ENDPOINT", ORANGE, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &clusterMergePoint, "clusterMergePoint", VIOLET, 2);
        const CartesianVector start(clusterMergePoint + (clusterMergeDirection*40));
        const CartesianVector end(clusterMergePoint - (clusterMergeDirection*40));
        PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &start, &end, "CLUSTER DIRECTION", GREEN, 4, 2);            
        std::cout << "deltaZ: " << deltaZ << std::endl;
        */
        /////////////////

        if (deltaZ < m_minScaledZOffset)
        {
            //std::cout << "DELTA Z IS TOO LOW" << std::endl;
            //PandoraMonitoringApi::ViewEvent(this->GetPandora());
            continue;
        }
           
        // SEARCH FOR DELTA RAY BEND    
        if(!this->IsDeltaRay(pCluster, clusterMergePoint, clusterMergeDirection, isEndUpstream))
        {
            //std::cout << "NOT A DELTA RAY" << std::endl;
            //PandoraMonitoringApi::ViewEvent(this->GetPandora());            
            continue;
        }

        const CartesianVector extrapolatedEndpointPosition(nearestTPCBoundaryX, 0.f, predictedIntercept + (predictedGradient * nearestTPCBoundaryX));

        // IS THIS EVER NEEDED? - I THINK IT CAN GO WRONG IF THE FIT DOESN'T MAKE SENSE
        if (isEndUpstream ? extrapolatedEndpointPosition.GetZ() > clusterMergePoint.GetZ() : extrapolatedEndpointPosition.GetZ() < clusterMergePoint.GetZ())
        {
            std::cout << "EXTRAPOLATED ENDPOINT IS NOT IN FORWARD DIRECTION" << std::endl;
            //PandoraMonitoringApi::ViewEvent(this->GetPandora());
            continue;
        }

        clusterAssociation = isEndUpstream ?
            ClusterEndpointAssociation(extrapolatedEndpointPosition, clusterMergeDirection, clusterMergePoint, clusterMergeDirection * (-1.f), pCluster, true) :
            ClusterEndpointAssociation(clusterMergePoint, clusterMergeDirection, extrapolatedEndpointPosition, clusterMergeDirection * (-1.f), pCluster, false);        

        ////////////////////////////////
        /*
        ClusterList mainCluster({clusterAssociation.GetMainTrackCluster()});
        const CartesianVector &upstream(clusterAssociation.GetUpstreamMergePoint());
        const CartesianVector &downstream(clusterAssociation.GetDownstreamMergePoint());
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &upstream, "UPSTREAM", VIOLET, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &downstream, "DOWNSTREAM", RED, 2);
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &mainCluster, "CLUSTER", BLACK);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        */
        ////////////////////////////////        
        
        return true;
    }

    //std::cout << "DID NOT FIND AN ASSOCIATION" << std::endl;
    return false;

}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ExtensionPastDeltaRayAlgorithm::IsDeltaRay(const Cluster *const pCluster, CartesianVector &clusterMergePoint,
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
                std::cout << "MERGE POINT END TOO SMOOTH" << std::endl;
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
        std::cout << "END DEVIATION ANGLE NOT BIG ENOUGH" << std::endl;
        std::cout << "endOpeningAngle: " << endOpeningAngle << std::endl;
        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
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
                std::cout << "ENDPOINT END TOO SMOOTH" << std::endl;
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

StatusCode ExtensionPastDeltaRayAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
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

    return TrackExtensionRefinementAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
