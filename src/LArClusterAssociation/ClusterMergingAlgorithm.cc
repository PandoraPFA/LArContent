/**
 *  @file   LArContent/src/LArClusterAssociation/ClusterMergingAlgorithm.cc
 * 
 *  @brief  Implementation of the cluster association algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArClusterAssociation/ClusterMergingAlgorithm.h"

#include "LArHelpers/LArClusterHelper.h"

using namespace pandora;

namespace lar
{

StatusCode ClusterMergingAlgorithm::Run()
{
    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentClusterList(*this, pClusterList));

    // Get a list of clean clusters
    ClusterVector clusterVector;
    this->GetListOfCleanClusters(pClusterList, clusterVector);
    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNOccupiedLayers);

// ClusterList tempList;
// tempList.insert(clusterVector.begin(), clusterVector.end());
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&tempList, "CleanClusters",  RED);
// PandoraMonitoringApi::ViewEvent();


    // Form the matrix of associations
    ClusterAssociationMatrix clusterAssociationMatrix;
    this->FillAssociationMatrix( clusterVector, clusterAssociationMatrix );


    // Generate the list of merges
    ClusterMergeMap clusterMergeMap;
     
    for (ClusterVector::const_iterator iterI = clusterVector.begin(), iterEndI = clusterVector.end(); iterI != iterEndI; ++iterI)
    {
        Cluster* pClusterI = *iterI;

        for (ClusterVector::const_iterator iterJ = iterI, iterEndJ = clusterVector.end(); iterJ != iterEndJ; ++iterJ)
        {
            Cluster* pClusterJ = *iterJ;

            if ( pClusterI == pClusterJ ) continue;
 
            if ( this->AreClustersAssociated(pClusterI, pClusterJ, clusterAssociationMatrix) )
	    {

// ClusterList tempList1, tempList2;
// tempList1.insert(pClusterI);
// tempList2.insert(pClusterJ);
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&tempList1, "Cluster1",  BLUE);
// PandoraMonitoringApi::VisualizeClusters(&tempList2, "Cluster2",  GREEN);
// PandoraMonitoringApi::ViewEvent();

	        clusterMergeMap[pClusterI].insert(pClusterJ);
                clusterMergeMap[pClusterJ].insert(pClusterI);
	    }
	}
    }

    // Make merges 
    ClusterVetoMap clusterVetoMap;

    for (ClusterVector::iterator iter1 = clusterVector.begin(), iterEnd1 = clusterVector.end(); iter1 != iterEnd1; ++iter1)
    {
        Cluster *pSeedCluster = *iter1;

        ClusterList mergeList;
        this->CollectAssociatedClusters(pSeedCluster, pSeedCluster, clusterMergeMap, clusterVetoMap, mergeList);

        for( ClusterList::iterator iter2 = mergeList.begin(), iterEnd2 = mergeList.end(); iter2 != iterEnd2; ++iter2 )
        {
            Cluster *pAssociatedCluster = *iter2;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pSeedCluster, pAssociatedCluster));
            clusterVetoMap[pAssociatedCluster] = true;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ClusterMergingAlgorithm::AreClustersAssociated( Cluster* pCluster1, Cluster* pCluster2, ClusterAssociationMatrix& clusterAssociationMatrix ) const
{
    // Search for the pCluster1 <-> pCluster2 association
    ClusterAssociationMatrix::iterator iter1A = clusterAssociationMatrix.find(pCluster1);

    if ( clusterAssociationMatrix.end() != iter1A ) 
    {
        ClusterAssociationMap& clusterAssociationMap1(iter1A->second);
        ClusterAssociationMap::iterator iter1B = clusterAssociationMap1.find(pCluster2);

        if ( clusterAssociationMap1.end() != iter1B )
	{
            ClusterAssociation& clusterAssociation12(iter1B->second);

            if ( clusterAssociation12.GetStrength() > ClusterAssociation::UNASSOCIATED )
	        return true;
	}
    }

    // Search for the pCluster2 <-> pCluster1 association
    ClusterAssociationMatrix::iterator iter2A = clusterAssociationMatrix.find(pCluster2);

    if ( clusterAssociationMatrix.end() != iter2A ) 
    {
        ClusterAssociationMap& clusterAssociationMap2(iter2A->second);
        ClusterAssociationMap::iterator iter2B = clusterAssociationMap2.find(pCluster1);

        if ( clusterAssociationMap2.end() != iter2B )
	{
            ClusterAssociation& clusterAssociation21(iter2B->second);

            if ( clusterAssociation21.GetStrength() > ClusterAssociation::UNASSOCIATED )
	        return true;
	}
    }


    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterMergingAlgorithm::CollectAssociatedClusters(Cluster *pSeedCluster, Cluster *pCurrentCluster, ClusterMergeMap& clusterMergeMap, ClusterVetoMap& clusterVetoMap, ClusterList &associatedClusterList)
{
    ClusterVetoMap::const_iterator iter0 = clusterVetoMap.find(pCurrentCluster);

    if (iter0 != clusterVetoMap.end())
        return;

    ClusterMergeMap::const_iterator iter1 = clusterMergeMap.find(pCurrentCluster);

    if (iter1 == clusterMergeMap.end())
        return;

    for (ClusterList::const_iterator iter2 = iter1->second.begin(), iterEnd2 = iter1->second.end(); iter2 != iterEnd2; ++iter2)
    {
        Cluster* pAssociatedCluster = *iter2; 

        if (pAssociatedCluster == pSeedCluster)
            continue;

        if (!associatedClusterList.insert(pAssociatedCluster).second)
            continue;

        this->CollectAssociatedClusters(pSeedCluster, pAssociatedCluster, clusterMergeMap, clusterVetoMap, associatedClusterList);
    }

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ClusterMergingAlgorithm::ClusterAssociation::ClusterAssociation() : 
    m_parent(NONE), m_daughter(NONE), m_association(NOTHING), m_strength(UNASSOCIATED), m_fom(0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

ClusterMergingAlgorithm::ClusterAssociation::ClusterAssociation( StrengthType strength ) :
    m_parent(NONE), m_daughter(NONE), m_association(NOTHING), m_strength(strength), m_fom(0.f)
{

}

//------------------------------------------------------------------------------------------------------------------------------------------

ClusterMergingAlgorithm::ClusterAssociation::ClusterAssociation( ClusterMergingAlgorithm::ClusterAssociation::VertexType parent, ClusterMergingAlgorithm::ClusterAssociation::VertexType daughter, ClusterMergingAlgorithm::ClusterAssociation::AssociationType association, ClusterMergingAlgorithm::ClusterAssociation::StrengthType strength, float fom) :
    m_parent(parent), m_daughter(daughter), m_association(association), m_strength(strength), m_fom(fom)
{
}
 
//------------------------------------------------------------------------------------------------------------------------------------------

ClusterMergingAlgorithm::ClusterAssociation::VertexType ClusterMergingAlgorithm::ClusterAssociation::GetParent()
{
    return m_parent;
}
  
//------------------------------------------------------------------------------------------------------------------------------------------
   
ClusterMergingAlgorithm::ClusterAssociation::VertexType ClusterMergingAlgorithm::ClusterAssociation::GetDaughter()
{
    return m_daughter;
}

//------------------------------------------------------------------------------------------------------------------------------------------

ClusterMergingAlgorithm::ClusterAssociation::AssociationType ClusterMergingAlgorithm::ClusterAssociation::GetAssociation()
{
    return m_association;
}
   
//------------------------------------------------------------------------------------------------------------------------------------------
      
ClusterMergingAlgorithm::ClusterAssociation::StrengthType ClusterMergingAlgorithm::ClusterAssociation::GetStrength()
{
    return m_strength;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ClusterMergingAlgorithm::ClusterAssociation::GetFigureOfMerit()
{
    return m_fom;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
