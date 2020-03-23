/**
 *  @file   EMTrackAlgorithm.cc
 *
 *  @brief  Implementation of the em track algorithm class
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArTwoDReco/EMTrackAlgorithm.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

using namespace pandora;

namespace lar_content
{


EMTrackAlgorithm::EMTrackAlgorithm() :
    m_minCaloHits(30)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EMTrackAlgorithm::Run()
{

    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

    // Get filtered clusters
    ClusterVector clusterVector;
    SelectCleanClusters(pClusterList, clusterVector);

    PandoraMonitoringApi::Create(this->GetPandora());
    PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_DEFAULT, -1.f, 1.f, 1.f);
    
    ClusterToAssociatedClusterMap clusterToAssociatedClusterMap;
    FillClusterToAssociatedClusterMap(clusterVector, clusterToAssociatedClusterMap);
    
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EMTrackAlgorithm::SelectCleanClusters(const ClusterList *pClusterList, ClusterVector &clusterVector)
{
    for (const Cluster *const pCluster : *pClusterList)
    {
        if (pCluster->GetNCaloHits() < m_minCaloHits)
            continue;

        clusterVector.push_back(pCluster);
    }
    
    // sort hits by Z, then X and then Y
    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------


void EMTrackAlgorithm::FillClusterToAssociatedClusterMap(const ClusterVector &clusterVector, ClusterToAssociatedClusterMap &clusterToAssociatedClusterMap)
{

    TwoDSlidingFitResultMap slidingFitResultMap;
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const float slidingFitWindow(20);
    
    for (const Cluster *const pCluster : clusterVector)
    {
        try {(void) slidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pCluster, TwoDSlidingFitResult(pCluster, slidingFitWindow, slidingFitPitch)));}
        catch (StatusCodeException &) {}
    }

    
    for (ClusterVector::const_iterator iterI = clusterVector.begin(); iterI != clusterVector.end(); ++iterI)
    {
        const Cluster *const pCurrentCluster = *iterI;
        TwoDSlidingFitResultMap::const_iterator fitIterI = slidingFitResultMap.find(pCurrentCluster);

        if (slidingFitResultMap.end() == fitIterI)
            continue;

        for (ClusterVector::const_iterator iterJ = iterI; iterJ != clusterVector.end(); ++iterJ)
        {
            const Cluster *const pTestCluster = *iterJ;

            if (pCurrentCluster == pTestCluster)
                continue;

            TwoDSlidingFitResultMap::const_iterator fitIterJ = slidingFitResultMap.find(pTestCluster);

            if (slidingFitResultMap.end() == fitIterJ)
                continue;

            if (!AreClustersAssociated(fitIterI->second, fitIterJ->second))
                continue;

            clusterToAssociatedClusterMap[pCurrentCluster] = pTestCluster;
        }
    }     
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EMTrackAlgorithm::AreClustersAssociated(const TwoDSlidingFitResult &currentClusterFit, const TwoDSlidingFitResult &testClusterFit)
{

    CartesianVector currentPoint(0,0,0), testPoint(0,0,0);
    CartesianVector currentDirection(0,0,0), testDirection(0,0,0);
    GetClusterMergingCoordinates(currentClusterFit, currentPoint, currentDirection, testClusterFit, testPoint, testDirection);

    ClusterList currentCluster;
    currentCluster.push_back(currentClusterFit.GetCluster());
    currentCluster.push_back(testClusterFit.GetCluster());
    //ClusterList testCluster;
    //testCluster.push_back(testClusterFit.GetCluster());   
    
    // check that opening angle is not too large
    if (std::fabs(currentDirection.GetCosOpeningAngle(testDirection)) < 0.97)
    {
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &currentCluster, "ANGLE: " + std::to_string(std::fabs(currentDirection.GetCosOpeningAngle(testDirection))), RED);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return false;
    }
    
    // check that clusters are reasonably far away
    if (std::sqrt(currentPoint.GetDistanceSquared(testPoint) < 50))
    {
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &currentCluster, "DISTANCE: " + std::to_string(std::sqrt(currentPoint.GetDistanceSquared(testPoint))), RED);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return false;
    }
        
    // check that fit allows you to get from one merge point to other merge point
    float separationDistance(std::sqrt(currentPoint.GetDistanceSquared(testPoint)));
    CartesianVector extrapolatedCurrentPoint(currentPoint + currentDirection*separationDistance);
    CartesianVector extrapolatedTestPoint(testPoint + testDirection*separationDistance);


    float m_maxXSeparation(5);
    if (extrapolatedCurrentPoint.GetX() > testPoint.GetX() + m_maxXSeparation || extrapolatedCurrentPoint.GetX() < testPoint.GetX() - m_maxXSeparation)
    {
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &currentPoint, "CURRENT MERGE POINT", BLUE, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedCurrentPoint, "EXTRAPOLATED CURRENT POINT", BLACK, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testPoint, "TEST MERGE POINT", RED, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedTestPoint, "EXTRAPOLATED TEST POINT", DARKGREEN, 2);

        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &currentCluster, "X TOO FAR", RED);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return false;
    }

    if (extrapolatedTestPoint.GetX() > currentPoint.GetX() + m_maxXSeparation || extrapolatedTestPoint.GetX() < currentPoint.GetX() - m_maxXSeparation)
    {
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &currentPoint, "CURRENT MERGE POINT", BLUE, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedCurrentPoint, "EXTRAPOLATED CURRENT POINT", BLACK, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testPoint, "TEST MERGE POINT", RED, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedTestPoint, "EXTRAPOLATED TEST POINT", DARKGREEN, 2);

        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &currentCluster, "X TOO FAR", RED);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return false;
    }

    float m_maxZSeparation(5);
    if (extrapolatedCurrentPoint.GetZ() > testPoint.GetZ() + m_maxZSeparation || extrapolatedCurrentPoint.GetZ() < testPoint.GetZ() - m_maxZSeparation)
    {
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &currentPoint, "CURRENT MERGE POINT", BLUE, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedCurrentPoint, "EXTRAPOLATED CURRENT POINT", BLACK, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testPoint, "TEST MERGE POINT", RED, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedTestPoint, "EXTRAPOLATED TEST POINT", DARKGREEN, 2);

        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &currentCluster, "Z TOO FAR", RED);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return false;
    }

    if (extrapolatedTestPoint.GetZ() > currentPoint.GetZ() + m_maxZSeparation || extrapolatedTestPoint.GetZ() < currentPoint.GetZ() - m_maxZSeparation)
    {
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &currentPoint, "CURRENT MERGE POINT", BLUE, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedCurrentPoint, "EXTRAPOLATED CURRENT POINT", BLACK, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testPoint, "TEST MERGE POINT", RED, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &extrapolatedTestPoint, "EXTRAPOLATED TEST POINT", DARKGREEN, 2);

        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &currentCluster, "Z TOO FAR", RED);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return false;
    }

    
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &currentCluster, "PASSED!", GREEN);




    
    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

    void EMTrackAlgorithm::GetClusterMergingCoordinates(const TwoDSlidingFitResult &cluster1FitResult, CartesianVector &cluster1MergePoint, CartesianVector &cluster1Direction, const TwoDSlidingFitResult &cluster2FitResult, CartesianVector &cluster2MergePoint, CartesianVector &cluster2Direction)
{

    CartesianVector minLayerPosition1(cluster1FitResult.GetGlobalMinLayerPosition()), maxLayerPosition1(cluster1FitResult.GetGlobalMaxLayerPosition());
    CartesianVector minLayerPosition2(cluster2FitResult.GetGlobalMinLayerPosition()), maxLayerPosition2(cluster2FitResult.GetGlobalMaxLayerPosition());

    CartesianVector minLayerDirection1(cluster1FitResult.GetGlobalMinLayerDirection()), maxLayerDirection1(cluster1FitResult.GetGlobalMaxLayerDirection());
    CartesianVector minLayerDirection2(cluster2FitResult.GetGlobalMinLayerDirection()), maxLayerDirection2(cluster2FitResult.GetGlobalMaxLayerDirection());

    float separationDistance(std::numeric_limits<float>::max());

     float ii(minLayerPosition1.GetDistanceSquared(minLayerPosition2));
     if (ii < separationDistance)
     {
         cluster1MergePoint = minLayerPosition1;
         cluster2MergePoint = minLayerPosition2;
         cluster1Direction = minLayerDirection1;
         cluster2Direction = minLayerDirection2;
         separationDistance = ii;
     }

     float io(minLayerPosition1.GetDistanceSquared(maxLayerPosition2));
     if (io < separationDistance)
     {
         cluster1MergePoint = minLayerPosition1;
         cluster2MergePoint = maxLayerPosition2;
         cluster1Direction = minLayerDirection1;
         cluster2Direction = maxLayerDirection2;
         separationDistance = io;
     }

     float oi(maxLayerPosition1.GetDistanceSquared(minLayerPosition2));
     if (oi < separationDistance)
     {
         cluster1MergePoint = maxLayerPosition1;
         cluster2MergePoint = minLayerPosition2;
         cluster1Direction = maxLayerDirection1;
         cluster2Direction = minLayerDirection2;
         separationDistance = oi;
     }

     float oo(maxLayerPosition1.GetDistanceSquared(maxLayerPosition2));
     if (oo < separationDistance)
     {
         cluster1MergePoint = maxLayerPosition1;
         cluster2MergePoint = maxLayerPosition2;
         cluster1Direction = maxLayerDirection1;
         cluster2Direction = maxLayerDirection2;
         separationDistance = oo;
     }                 
     
     //make sure directions are pointing at each other
     CartesianVector cluster1To2(cluster2MergePoint - cluster1MergePoint);
     CartesianVector cluster2To1(cluster1MergePoint - cluster2MergePoint);
     if (cluster1To2.GetCosOpeningAngle(cluster1Direction*(-1)) > cluster1To2.GetCosOpeningAngle(cluster1Direction))
         cluster1Direction *= -1.0;

     if (cluster2To1.GetCosOpeningAngle(cluster2Direction*(-1)) > cluster2To1.GetCosOpeningAngle(cluster2Direction))
         cluster2Direction *= -1.0;
     
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
StatusCode EMTrackAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{


    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCaloHits", m_minCaloHits));
    
    return STATUS_CODE_SUCCESS;
}


} // namespace lar_content
