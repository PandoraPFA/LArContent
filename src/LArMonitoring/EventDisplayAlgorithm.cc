/**
 *  @file   LArContent/src/LArMonitoring/EventDisplayAlgorithm.cc
 *
 *  @brief  Implementation of the event display algorithm
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArMonitoring/EventDisplayAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode EventDisplayAlgorithm::Run()
{
    const ClusterList *pClusterList = NULL;
    const PfoList     *pPfoList = NULL;

    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_clusterListName, pClusterList));

    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_particleListName, pPfoList));

    if (NULL == pClusterList)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));
        std::cout << " --- Loading Current Cluster List --- " << std::endl;
    }


    // Hit type
    HitType hitType(HIT_CUSTOM);


    // Add the seed clusters
    ClusterList clusterList;

    if( NULL != pClusterList )
    {
        for ( ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter )
        {
            Cluster *pCluster = *iter;
            const HitType clusterHitType(LArClusterHelper::GetClusterHitType(pCluster));

            if (HIT_CUSTOM == hitType)
            {
                if ((TPC_VIEW_U != clusterHitType) && (TPC_VIEW_V != clusterHitType) && (TPC_VIEW_W != clusterHitType))
                    throw StatusCodeException(STATUS_CODE_FAILURE);

                hitType = clusterHitType;
            }
            else if (hitType != clusterHitType)
            {
                throw StatusCodeException(STATUS_CODE_FAILURE);
            }

            if (pCluster->IsAvailable())
                clusterList.insert(pCluster);
        }
    }


    // Start Drawing Stuff
#ifdef MONITORING
    unsigned int n(0);
#endif

    if ( NULL != pPfoList )
    {
        PfoVector pfoVector(pPfoList->begin(), pPfoList->end());

        for (PfoVector::iterator pIter = pfoVector.begin(), pIterEnd = pfoVector.end(); pIter != pIterEnd; ++pIter)
        {
            ParticleFlowObject *pPfo = *pIter;

            ClusterList particleClusterList;

            for (ClusterList::const_iterator cIter = pPfo->GetClusterList().begin(), cIterEnd = pPfo->GetClusterList().end(); cIter != cIterEnd; ++cIter)
            {
                Cluster *pCluster = *cIter;

                if (hitType == LArClusterHelper::GetClusterHitType(pCluster))
                    particleClusterList.insert(pCluster);
            }

            PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &particleClusterList, "Particles", GetColor(n++) ));
        }

        PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clusterList, "Clusters", GRAY ));
    }

    else if( NULL != pClusterList )
    {
        for( ClusterList::const_iterator iter = clusterList.begin(), iterEnd = clusterList.end(); iter != iterEnd; ++iter )
        {
            Cluster *pCluster = *iter;
            ClusterList tempList;
            tempList.insert(pCluster);
            PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &tempList, "Clusters", GetColor(n++) ));
        }
    }

    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

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
    m_clusterListName = "ClusterListName";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterListName", m_clusterListName));

    m_particleListName = "ParticleListName";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ParticleListName", m_particleListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
