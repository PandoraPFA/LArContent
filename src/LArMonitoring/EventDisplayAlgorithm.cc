/**
 *  @file   LArContent/src/LArMonitoring/EventDisplayAlgorithm.cc
 * 
 *  @brief  Implementation of the event display algorithm
 * 
 *  $Log: $
 */

#ifdef MONITORING

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
    const PfoList     *pPfoList = NULL;

    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetClusterList(*this, m_seedClusterListName, pSeedClusterList));

    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetClusterList(*this, m_nonSeedClusterListName, pNonSeedClusterList));

    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetPfoList(*this, m_particleListName, pPfoList));

    if (NULL == pSeedClusterList && NULL == pNonSeedClusterList)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentClusterList(*this, pNonSeedClusterList));  
        std::cout << " --- Loading Current Cluster List --- " << std::endl;
    }


    // Hit type
    HitType hitType(CUSTOM); 


    // Add the seed clusters
    ClusterList seedClusterList;

    if( NULL != pSeedClusterList )
    {
        for ( ClusterList::const_iterator iter = pSeedClusterList->begin(), iterEnd = pSeedClusterList->end(); iter != iterEnd; ++iter ) 
        {
            Cluster *pCluster = *iter;

            if (hitType == CUSTOM && pCluster->ContainsHitType(VIEW_U)) 
                hitType = VIEW_U;
                
            if (hitType == CUSTOM && pCluster->ContainsHitType(VIEW_V)) 
                hitType = VIEW_V;
                
            if (hitType == CUSTOM && pCluster->ContainsHitType(VIEW_W)) 
                hitType = VIEW_W;

            if (pCluster->ContainsHitType(hitType) == false)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            if (pCluster->IsAvailable())
                seedClusterList.insert(pCluster);
        }
    }

  
    // Add the non-seed clusters
    ClusterList nonSeedClusterList; 

    if( NULL != pNonSeedClusterList )
    {
        for ( ClusterList::const_iterator iter = pNonSeedClusterList->begin(), iterEnd = pNonSeedClusterList->end(); iter != iterEnd; ++iter ) 
        {
            Cluster *pCluster = *iter;

            if (hitType == CUSTOM && pCluster->ContainsHitType(VIEW_U)) 
                hitType = VIEW_U;
                
            if (hitType == CUSTOM && pCluster->ContainsHitType(VIEW_V)) 
                hitType = VIEW_V;
                
            if (hitType == CUSTOM && pCluster->ContainsHitType(VIEW_W)) 
                hitType = VIEW_W;

            if (pCluster->ContainsHitType(hitType) == false)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            if (pCluster->IsAvailable())
                nonSeedClusterList.insert(pCluster);
        }
    }

   
    // Start Drawing Stuff
    unsigned int n(0);

    PANDORA_MONITORING_API(SetEveDisplayParameters(0, 0, -1.f, 1.f));

    if ( NULL != pPfoList )
    {
        PfoVector pfoVector(pPfoList->begin(), pPfoList->end());

        for (PfoVector::iterator pIter = pfoVector.begin(), pIterEnd = pfoVector.end(); pIter != pIterEnd; ++pIter)
        {
            ParticleFlowObject *pPfo = *pIter;

            ClusterList clusterList;

            for (ClusterList::const_iterator cIter = pPfo->GetClusterList().begin(), cIterEnd = pPfo->GetClusterList().end(); cIter != cIterEnd; ++cIter)
            {
                Cluster *pCluster = *cIter;

                if ( pCluster->ContainsHitType(hitType) )
                    clusterList.insert(pCluster);
            }

            PANDORA_MONITORING_API(VisualizeClusters(&clusterList, "Particles", GetColor(n++) ));
        }

        PANDORA_MONITORING_API(VisualizeClusters(&seedClusterList, "Seeds", GRAY ));
        PANDORA_MONITORING_API(VisualizeClusters(&nonSeedClusterList, "NonSeeds", GRAY ));
    }

    else if( NULL != pSeedClusterList )
    {
        for( ClusterList::const_iterator iter = seedClusterList.begin(), iterEnd = seedClusterList.end(); iter != iterEnd; ++iter )
        {
            Cluster* pCluster = *iter; 
            ClusterList tempList;
            tempList.insert(pCluster);
            PANDORA_MONITORING_API(VisualizeClusters(&tempList, "Seed", GetColor(n++) ));
        }

        PANDORA_MONITORING_API(VisualizeClusters(&nonSeedClusterList, "NonSeeds", GRAY ));
    }

    else
    {
        for( ClusterList::const_iterator iter = nonSeedClusterList.begin(), iterEnd = nonSeedClusterList.end(); iter != iterEnd; ++iter )
        {
            Cluster* pCluster = *iter; 
            ClusterList tempList;
            tempList.insert(pCluster);
            PANDORA_MONITORING_API(VisualizeClusters(&tempList, "NonSeed", GetColor(n++) ));
        }
    }

    if( LArVertexHelper::DoesVertexExist(m_vertexName) ){
      CartesianVector theVertex = LArVertexHelper::GetVertex(m_vertexName);
      PANDORA_MONITORING_API(AddMarkerToVisualization(&theVertex, "Vertex", BLACK, 2.0));
    }

    PANDORA_MONITORING_API(ViewEvent());


    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

#ifdef MONITORING
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
#endif

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventDisplayAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_seedClusterListName = "SeedClusterListName";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "SeedClusterListName", m_seedClusterListName));

    m_nonSeedClusterListName = "NonSeedClusterListName";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "NonSeedClusterListName", m_nonSeedClusterListName));

    m_vertexName = "VertexName";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "VertexName", m_vertexName));

    m_particleListName = "ParticleListName";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "ParticleListName", m_particleListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar

#endif // #ifdef MONITORING
