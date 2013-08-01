/**
 *  @file   LArContent/src/LArMonitoring/EventDisplayAlgorithm.cc
 * 
 *  @brief  Implementation of the event display algorithm
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArVertexHelper.h"

#include "LArMonitoring/EventDisplayAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode EventDisplayAlgorithm::Run()
{ 
    const ClusterList *pSeedClusterList = NULL;
    const ClusterList *pNonSeedClusterList = NULL;

    if (STATUS_CODE_SUCCESS == PandoraContentApi::GetClusterList(*this, m_seedClusterListName, pSeedClusterList))
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetClusterList(*this, m_seedClusterListName, pSeedClusterList));
    }

    if (STATUS_CODE_SUCCESS == PandoraContentApi::GetClusterList(*this, m_nonSeedClusterListName, pNonSeedClusterList))
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetClusterList(*this, m_nonSeedClusterListName, pNonSeedClusterList));
    }
    else
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentClusterList(*this, pNonSeedClusterList));
    }

    ClusterVector seedClusterVector, nonSeedClusterVector;

    // Add the seed clusters
    if( NULL != pSeedClusterList )
    {
        for ( ClusterList::const_iterator iter = pSeedClusterList->begin(), iterEnd = pSeedClusterList->end(); iter != iterEnd; ++iter ) 
            seedClusterVector.push_back(*iter);
        std::sort(seedClusterVector.begin(), seedClusterVector.end(), LArClusterHelper::SortByNOccupiedLayers);
    }

    // Add the non-seed clusters
    if( NULL != pNonSeedClusterList )
    {
        for ( ClusterList::const_iterator iter = pNonSeedClusterList->begin(), iterEnd = pNonSeedClusterList->end(); iter != iterEnd; ++iter ) 
            nonSeedClusterVector.push_back(*iter);
        std::sort(nonSeedClusterVector.begin(), nonSeedClusterVector.end(), LArClusterHelper::SortByNOccupiedLayers);
    }

    
    CartesianVector theVertex = LArVertexHelper::GetVertex(m_vertexName);

    PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);

    for( unsigned int n=0; n<seedClusterVector.size(); ++n ){
      Cluster* pCluster = seedClusterVector.at(n);
      ClusterList clusterList;
      clusterList.insert(pCluster);
      PandoraMonitoringApi::VisualizeClusters(&clusterList, "Cluster", GetColor(n) );
    } 

    for( unsigned int n=0; n<nonSeedClusterVector.size(); ++n ){
      Cluster* pCluster = nonSeedClusterVector.at(n);
      ClusterList clusterList;
      clusterList.insert(pCluster);
      PandoraMonitoringApi::VisualizeClusters(&clusterList, "Cluster", GRAY );
    } 

    PandoraMonitoringApi::AddMarkerToVisualization(&theVertex, "Vertex", BLACK, 2.0);

    PandoraMonitoringApi::ViewEvent();


    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

Color EventDisplayAlgorithm::GetColor( unsigned int icolor )
{
    unsigned thisColor = icolor%10;


    switch( thisColor ){
    case 0: return RED;
    case 1: return BLUE;
    case 2: return DARKGREEN;
    case 3: return DARKMAGENTA;
    case 4: return DARKCYAN;
    case 5: return BLACK;
    case 6: return GRAY;
    case 7: return TEAL;
    case 8: return AZURE;
    case 9: return DARKYELLOW;
    default: return YELLOW;
    }

    return YELLOW;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventDisplayAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "SeedClusterListName", m_seedClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NonSeedClusterListName", m_nonSeedClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "VertexName", m_vertexName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
