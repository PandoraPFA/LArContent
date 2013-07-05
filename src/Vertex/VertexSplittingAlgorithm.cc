/**
 *  @file   VertexSplittingAlgorithm.cc
 * 
 *  @brief  Implementation of the vertex splitting algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArClusterHelper.h"
#include "LArGeometryHelper.h"

#include "VertexSplittingAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode VertexSplittingAlgorithm::Run()
{
    std::cout << " *** VertexSplittingAlgorithm::Run() *** " << std::endl;
    


    // Obtain sorted vectors of clusters for this view
    const ClusterList *pClusterList(NULL);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetClusterList(*this, m_inputClusterListName, pClusterList));
    


    ClusterVector clusterVector;

    this->GetListOfCleanClusters( pClusterList, clusterVector );
   
    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNOccupiedLayers);
  

    ClusterList tempList;

    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
        tempList.insert(*iter);

 


    

    PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
    PandoraMonitoringApi::VisualizeClusters(&tempList, "ClusterList", RED);
    PandoraMonitoringApi::ViewEvent();



    



    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSplittingAlgorithm::GetListOfCleanClusters(const pandora::ClusterList* const pClusterList, pandora::ClusterVector &clusterVector)
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

StatusCode VertexSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListName", m_inputClusterListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
