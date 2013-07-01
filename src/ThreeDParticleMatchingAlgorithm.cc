/**
 *  @file   ThreeDParticleMatchingAlgorithm.cc
 * 
 *  @brief  Implementation of the three dimensional particle creation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArClusterHelper.h"
#include "LArGeometryHelper.h"

#include "ThreeDParticleMatchingAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode ThreeDParticleMatchingAlgorithm::Run()
{
    std::cout << " *** ThreeDParticleMatchingAlgorithm::Run() *** " << std::endl;
    


    // Obtain sorted vectors of seed clusters in U, V and W views
    const ClusterList *pClusterListU(NULL), *pClusterListV(NULL), *pClusterListW(NULL);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetClusterList(*this, m_inputClusterListNameU, pClusterListU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetClusterList(*this, m_inputClusterListNameV, pClusterListV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetClusterList(*this, m_inputClusterListNameW, pClusterListW));

    

   


    LArPointingClusterVertexList clusterVertexListU, clusterVertexListV, clusterVertexListW;
   
    ProcessSingleView( pClusterListU, clusterVertexListU );
    ProcessSingleView( pClusterListV, clusterVertexListV );
    ProcessSingleView( pClusterListW, clusterVertexListW );

    


    ClusterVector clusterVectorU, clusterVectorV, clusterVectorW;

    this->GetListOfCleanClusters( pClusterListU, clusterVectorU );
    this->GetListOfCleanClusters( pClusterListV, clusterVectorV );
    this->GetListOfCleanClusters( pClusterListW, clusterVectorW );

    std::sort(clusterVectorU.begin(), clusterVectorU.end(), LArClusterHelper::SortByNOccupiedLayers);
    std::sort(clusterVectorV.begin(), clusterVectorV.end(), LArClusterHelper::SortByNOccupiedLayers);
    std::sort(clusterVectorW.begin(), clusterVectorW.end(), LArClusterHelper::SortByNOccupiedLayers);


    ClusterList tempListU, tempListV, tempListW;

    for (ClusterVector::const_iterator iterU = clusterVectorU.begin(), iterUEnd = clusterVectorU.end(); iterU != iterUEnd; ++iterU)
        tempListU.insert(*iterU);

    for (ClusterVector::const_iterator iterV = clusterVectorV.begin(), iterVEnd = clusterVectorV.end(); iterV != iterVEnd; ++iterV)
        tempListV.insert(*iterV);

    for (ClusterVector::const_iterator iterW = clusterVectorW.begin(), iterWEnd = clusterVectorW.end(); iterW != iterWEnd; ++iterW)
        tempListW.insert(*iterW);


    

    PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
    PandoraMonitoringApi::VisualizeClusters(&tempListU, "ClusterListU", RED);
    PandoraMonitoringApi::VisualizeClusters(&tempListV, "ClusterListV", BLUE);
    PandoraMonitoringApi::VisualizeClusters(&tempListW, "ClusterListW", GREEN);
    PandoraMonitoringApi::ViewEvent();



    






    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDParticleMatchingAlgorithm::ProcessSingleView( const pandora::ClusterList* const pClusterList, LArPointingClusterVertexList& pointingClusterVertexList )
{
 
    // Select a set of clean clusters
    ClusterVector clusterVector;
    this->GetListOfCleanClusters(pClusterList, clusterVector);
    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNOccupiedLayers);


    // Generate a list of pointing clusters
    LArPointingClusterMap pointingClusterMap;
    LArPointingClusterList pointingClusterList;

    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
    {
        pointingClusterMap.insert( std::pair<Cluster*,LArPointingCluster>(*iter,LArPointingCluster(*iter)) );
    }

    for (LArPointingClusterMap::const_iterator iter = pointingClusterMap.begin(), iterEnd = pointingClusterMap.end(); iter != iterEnd; ++iter)
    {
        const LArPointingCluster& pointingCluster = iter->second;
        pointingClusterList.push_back(pointingCluster);
        pointingClusterVertexList.push_back(pointingCluster.GetInnerVertex());
        pointingClusterVertexList.push_back(pointingCluster.GetOuterVertex());
    }

  


}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDParticleMatchingAlgorithm::GetListOfCleanClusters(const pandora::ClusterList* const pClusterList, pandora::ClusterVector &clusterVector)
{
   for ( ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter ) 
    {
        Cluster *pCluster = *iter;

        if ( 1 + pCluster->GetOuterPseudoLayer() - pCluster->GetInnerPseudoLayer() < 10 ) continue;

        const CartesianVector innerVertex = pCluster->GetCentroid( pCluster->GetInnerPseudoLayer() );
        const CartesianVector outerVertex = pCluster->GetCentroid( pCluster->GetOuterPseudoLayer() );

        if ( (outerVertex-innerVertex).GetMagnitudeSquared() < 15.f ) continue;

        if ( LArClusterHelper::GetLayerOccupancy(pCluster) < 0.75 ) continue;

        clusterVector.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDParticleMatchingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameU", m_inputClusterListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameV", m_inputClusterListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameW", m_inputClusterListNameW));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
