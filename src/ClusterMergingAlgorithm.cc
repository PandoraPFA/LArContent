/**
 *  @file   ClusterMergingAlgorithm.cc
 * 
 *  @brief  Implementation of the cluster association algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "ClusterMergingAlgorithm.h"

#include "LArClusterHelper.h"

using namespace pandora;

namespace lar
{

StatusCode ClusterMergingAlgorithm::Run()
{
    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentClusterList(*this, pClusterList));

    ClusterVector clusterVector;
    this->GetListOfCleanClusters(pClusterList, clusterVector);
    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNOccupiedLayers);

    // Form associations between clusters
    m_clusterMergeMap.clear();

    for (ClusterVector::const_iterator iterI = clusterVector.begin(), iterEndI = clusterVector.end(); iterI != iterEndI; ++iterI)
    {
        Cluster *pClusterI = *iterI;

        for (ClusterVector::const_iterator iterJ = iterI, iterEndJ = clusterVector.end(); iterJ != iterEndJ; ++iterJ)
        {
            Cluster *pClusterJ = *iterJ;

            if (pClusterI == pClusterJ)
                continue;

            // Check on Association (note: symmetrical)
            if (this->AreClustersAssociated(pClusterI, pClusterJ))
            {
   
// ClusterList tempList1, tempList2;
// tempList1.insert(pClusterI);
// tempList2.insert(pClusterJ);
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&tempList1, "Cluster1",  BLUE);
// PandoraMonitoringApi::VisualizeClusters(&tempList2, "Cluster2",  GREEN);
// PandoraMonitoringApi::ViewEvent();
                m_clusterMergeMap[pClusterI].insert(pClusterJ);
                m_clusterMergeMap[pClusterJ].insert(pClusterI);
            }
        }
    }

    // Make merges 
    m_clusterVetoMap.clear();

    for (ClusterVector::iterator iter1 = clusterVector.begin(), iterEnd1 = clusterVector.end(); iter1 != iterEnd1; ++iter1)
    {
        Cluster *pSeedCluster = *iter1;

        ClusterList mergeList;
        this->CollectAssociatedClusters(pSeedCluster, pSeedCluster, mergeList);

        for( ClusterList::iterator iter2 = mergeList.begin(), iterEnd2 = mergeList.end(); iter2 != iterEnd2; ++iter2 )
        {
            Cluster *pAssociatedCluster = *iter2;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pSeedCluster, pAssociatedCluster));
            m_clusterVetoMap[pAssociatedCluster] = true;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterMergingAlgorithm::CollectAssociatedClusters(Cluster *pSeedCluster, Cluster *pCurrentCluster, ClusterList &associatedClusterList)
{
    ClusterVetoMap::const_iterator iter0 = m_clusterVetoMap.find(pCurrentCluster);

    if (iter0 != m_clusterVetoMap.end())
        return;

    ClusterMergeMap::const_iterator iter1 = m_clusterMergeMap.find(pCurrentCluster);

    if (iter1 == m_clusterMergeMap.end())
        return;

    for (ClusterList::const_iterator iter2 = iter1->second.begin(), iterEnd2 = iter1->second.end(); iter2 != iterEnd2; ++iter2)
    {
        Cluster* pAssociatedCluster = *iter2; 

        if (pAssociatedCluster == pSeedCluster)
            continue;

        if (!associatedClusterList.insert(pAssociatedCluster).second)
            continue;

        this->CollectAssociatedClusters(pSeedCluster, pAssociatedCluster, associatedClusterList);
    }

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ClusterMergingAlgorithm::AreClustersAssociated(const Cluster *const pClusterI, const Cluster *const pClusterJ)
{
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterMergingAlgorithm::GetListOfCleanClusters( const ClusterList *const pClusterList, ClusterVector &clusterVector )
{
    for ( ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter ) 
    {
        Cluster *pCluster = *iter;

        if ( 1 + pCluster->GetOuterPseudoLayer() - pCluster->GetInnerPseudoLayer() < 15 ) continue;

        const CartesianVector innerVertex = pCluster->GetCentroid( pCluster->GetInnerPseudoLayer() );
        const CartesianVector outerVertex = pCluster->GetCentroid( pCluster->GetOuterPseudoLayer() );

        if ( (outerVertex-innerVertex).GetMagnitudeSquared() < 25.0 ) continue;

        if ( LArClusterHelper::GetLayerOccupancy(pCluster) < 0.75 ) continue;

        clusterVector.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
