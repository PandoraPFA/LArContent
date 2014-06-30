/**
 *  @file   LArContent/src/LArThreeDReco/LArClusterSplitting/ThreeDTrackSplittingAlgorithm.cc
 * 
 *  @brief  Implementation of the cosmic ray identification algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArThreeDReco/LArClusterSplitting/ThreeDTrackSplittingAlgorithm.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArPointingClusterHelper.h"
#include "LArHelpers/LArMCParticleHelper.h"

using namespace pandora;

namespace lar
{

StatusCode ThreeDTrackSplittingAlgorithm::Run()
{
    std::cout << " --- ThreeDTrackSplittingAlgorithm::Run() --- " << std::endl;



    // Get the cluster lists for each view
    const ClusterList *pClusterListU = NULL;    
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_clusterListNameU, pClusterListU));

    const ClusterList *pClusterListV = NULL;    
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_clusterListNameV, pClusterListV));

    const ClusterList *pClusterListW = NULL;    
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_clusterListNameW, pClusterListW));


    std::cout << " ThreeDTrackSplittingAlgorithm::Run() Running..." << std::endl;


    
    // Select a list of vertex clusters
    LArPointingClusterMap pointingClusterMapU, pointingClusterMapV, pointingClusterMapW;
    LArPointingClusterVertexList pointingVertexListU, pointingVertexListV, pointingVertexListW;

    this->GetListOfCleanPointingClusters(pClusterListU, pointingClusterMapU, pointingVertexListU);
    this->GetListOfCleanPointingClusters(pClusterListV, pointingClusterMapV, pointingVertexListV);
    this->GetListOfCleanPointingClusters(pClusterListW, pointingClusterMapW, pointingVertexListW);







// --- BEGIN EVENT DISPLAY ---
ClusterList tempClusterListU, tempClusterListV, tempClusterListW;
CartesianPointList tempMarkerListU, tempMarkerListV, tempMarkerListW;

this->CollectClusters(pointingClusterMapU, tempClusterListU);
this->CollectClusters(pointingClusterMapV, tempClusterListV);
this->CollectClusters(pointingClusterMapW, tempClusterListW);

this->CollectMarkers(pointingVertexListU, tempMarkerListU);
this->CollectMarkers(pointingVertexListV, tempMarkerListV);
this->CollectMarkers(pointingVertexListW, tempMarkerListW);

PANDORA_MONITORING_API(VisualizeClusters(&tempClusterListU, "Clusters (U)", GREEN));
PANDORA_MONITORING_API(VisualizeClusters(&tempClusterListV, "Clusters (V)", BLUE));
PANDORA_MONITORING_API(VisualizeClusters(&tempClusterListW, "Clusters (W)", RED));

for(CartesianPointList::const_iterator iterU = tempMarkerListU.begin(), iterEndU = tempMarkerListU.end(); iterU != iterEndU; ++iterU){
const CartesianVector& positionU = *iterU;
PANDORA_MONITORING_API(AddMarkerToVisualization(&positionU, "vertex (U)", GREEN, 2));
}

for(CartesianPointList::const_iterator iterV = tempMarkerListV.begin(), iterEndV = tempMarkerListV.end(); iterV != iterEndV; ++iterV){
const CartesianVector& positionV = *iterV;
PANDORA_MONITORING_API(AddMarkerToVisualization(&positionV, "vertex (V)", BLUE, 2));
}

for(CartesianPointList::const_iterator iterW = tempMarkerListW.begin(), iterEndW = tempMarkerListW.end(); iterW != iterEndW; ++iterW){
const CartesianVector& positionW = *iterW;
PANDORA_MONITORING_API(AddMarkerToVisualization(&positionW, "vertex (W)", RED, 2));
}

PANDORA_MONITORING_API(ViewEvent());
// --- END EVENT DISPLAY ---


    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDTrackSplittingAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    for ( ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter ) 
    {
        Cluster *pCluster = *iter;

        if( LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLength * m_minClusterLength)
            continue;

        clusterVector.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDTrackSplittingAlgorithm::GetListOfCleanPointingClusters(const ClusterList *const pClusterList, 
    LArPointingClusterMap& pointingClusterMap, LArPointingClusterVertexList& pointingVertexList) const
{  
    // Generate list of clean clusters
    ClusterVector clusterVector;
    this->GetListOfCleanClusters(pClusterList, clusterVector);
    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);

    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
    {
        LArPointingCluster pointingCluster(*iter);
        pointingClusterMap.insert(std::pair<Cluster*,LArPointingCluster>(*iter, pointingCluster));
        pointingVertexList.push_back(pointingCluster.GetInnerVertex());
        pointingVertexList.push_back(pointingCluster.GetOuterVertex());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
 
void ThreeDTrackSplittingAlgorithm::CollectClusters(const LArPointingClusterMap& inputList, ClusterList& outputList) const
{
    for (LArPointingClusterMap::const_iterator iter = inputList.begin(), iterEnd = inputList.end(); iter != iterEnd; ++iter)
        outputList.insert((iter->second).GetCluster());
}

//------------------------------------------------------------------------------------------------------------------------------------------
 
void ThreeDTrackSplittingAlgorithm::CollectClusters(const LArPointingClusterVertexList& inputList, ClusterList& outputList) const
{
    for (LArPointingClusterVertexList::const_iterator iter = inputList.begin(), iterEnd = inputList.end(); iter != iterEnd; ++iter)
        outputList.insert((*iter).GetCluster());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDTrackSplittingAlgorithm::CollectMarkers(const LArPointingClusterVertexList& inputList, CartesianPointList& outputList) const
{
    for (LArPointingClusterVertexList::const_iterator iter = inputList.begin(), iterEnd = inputList.end(); iter != iterEnd; ++iter)
        outputList.push_back((*iter).GetPosition());
}
 
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDTrackSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_clusterListNameU = "";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "ClusterListNameU", m_clusterListNameU));

    m_clusterListNameV = "";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "ClusterListNameV", m_clusterListNameV));

    m_clusterListNameW = "";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "ClusterListNameW", m_clusterListNameW));

    m_minClusterLength = 7.5f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLength", m_minClusterLength));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
