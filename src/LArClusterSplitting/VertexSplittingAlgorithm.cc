/**
 *  @file   LArContent/src/LArVertex/VertexSplittingAlgorithm.cc
 * 
 *  @brief  Implementation of the vertex splitting algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"
#include "LArHelpers/LArVertexHelper.h"

#include "LArClusterSplitting/VertexSplittingAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode VertexSplittingAlgorithm::Run()
{
    
    




    // Obtain sorted vectors of clusters for this view
    const ClusterList *pClusterList(NULL);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetClusterList(*this, m_inputClusterListName, pClusterList));
   

    std::list<Cluster*> internalClusterList(pClusterList->begin(), pClusterList->end());


    for (std::list<Cluster*>::iterator iter = internalClusterList.begin(); iter != internalClusterList.end(); ++iter)
    {  
        Cluster* pCluster = *iter;

        if (!this->IsPossibleSplit(pCluster))
            continue;

        unsigned int splitLayer(std::numeric_limits<unsigned int>::max());

        if (STATUS_CODE_SUCCESS != this->FindBestSplitLayer(pCluster,splitLayer))
            continue;

	if ((splitLayer <= pCluster->GetInnerPseudoLayer()) || (splitLayer >= pCluster->GetOuterPseudoLayer()))
            continue;

        std::list<Cluster*> daughters;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->SplitCluster(pCluster, splitLayer, daughters));	  

// ClusterList tempList(daughters.begin(),daughters.end());
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&tempList, "SplitCluster", AUTOITER);
// PandoraMonitoringApi::ViewEvent();

	*iter = NULL;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexSplittingAlgorithm::IsPossibleSplit(const Cluster *const pCluster) const
{
    if ( LArClusterHelper::GetLengthSquared(pCluster) < 4.0 * m_minSplitDisplacementSquared )
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexSplittingAlgorithm::FindBestSplitLayer(const Cluster* const pCluster, unsigned int& splitLayer )
{ 
    // Nothing to do
    if ( LArVertexHelper::DoesCurrentVertexExist() == false )
        return STATUS_CODE_NOT_FOUND;

    // Find nearest hit to vertex
    const CartesianVector &theVertex(LArVertexHelper::GetCurrentVertex());
    const CartesianVector innerCentroid(pCluster->GetCentroid(pCluster->GetInnerPseudoLayer()));
    const CartesianVector outerCentroid(pCluster->GetCentroid(pCluster->GetOuterPseudoLayer()));

    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    bool foundSplit(false);
    CartesianVector splitPosition(0.f,0.f,0.f);
    float minDisplacementSquared(m_minSplitDisplacementSquared);

    for (OrderedCaloHitList::const_iterator iterI = orderedCaloHitList.begin(), iterEndI = orderedCaloHitList.end(); iterI != iterEndI; ++iterI)
    {
        const unsigned int thisLayer    = iterI->first;
        const CaloHitList *pCaloHitList = iterI->second;

        for (CaloHitList::const_iterator iterJ = pCaloHitList->begin(), iterEndJ = pCaloHitList->end(); iterJ != iterEndJ; ++iterJ)
        {
            const CaloHit* pCaloHit = *iterJ;

            float thisDisplacementSquared((pCaloHit->GetPositionVector() - theVertex).GetMagnitudeSquared());

            if( thisDisplacementSquared > minDisplacementSquared )
	        continue;
	    
            foundSplit = false;
            minDisplacementSquared = thisDisplacementSquared;

            if( (pCaloHit->GetPositionVector()-innerCentroid).GetMagnitudeSquared() < m_minSplitDisplacementSquared
             || (pCaloHit->GetPositionVector()-outerCentroid).GetMagnitudeSquared() < m_minSplitDisplacementSquared )
	        continue;

            splitPosition = pCaloHit->GetPositionVector();
            splitLayer = thisLayer;
            foundSplit = true;	    
	}
    }

    if ( false == foundSplit ) 
        return STATUS_CODE_NOT_FOUND;


// ClusterList tempList;
// Cluster* tempCluster = (Cluster*)pCluster;
// tempList.insert(tempCluster);
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&tempList, "Cluster", GREEN);
// PandoraMonitoringApi::AddMarkerToVisualization(&theVertex, "Vertex", RED, 1.75); 
// PandoraMonitoringApi::AddMarkerToVisualization(&splitPosition, "Split", BLUE, 1.75);
// PandoraMonitoringApi::ViewEvent();


    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexSplittingAlgorithm::SplitCluster(Cluster *const pCluster, const unsigned int splitLayer, std::list<Cluster*>& daughters)
{
    // Begin cluster fragmentation operations
    ClusterList clusterList;
    clusterList.insert(pCluster);
    std::string clusterListToSaveName, clusterListToDeleteName;

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*this, clusterList, clusterListToDeleteName,
        clusterListToSaveName));

    // Create new clusters
    Cluster *pCluster1(NULL);
    Cluster *pCluster2(NULL);

    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(); iter != orderedCaloHitList.end(); ++iter)
    {
        const unsigned int thisLayer(iter->first);

        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            CaloHit *pCaloHit = *hitIter;
            Cluster *&pClusterToModify((thisLayer < splitLayer) ? pCluster1 : pCluster2);

            if (NULL == pClusterToModify)
            {
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, pCaloHit, pClusterToModify));
		daughters.push_back(pClusterToModify);
            }
            else
            {
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddCaloHitToCluster(*this, pClusterToModify, pCaloHit));
            }
        }
    }

    // End cluster fragmentation operations
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, clusterListToSaveName, clusterListToDeleteName));

    return STATUS_CODE_SUCCESS;
}



//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListName", m_inputClusterListName));

    m_minSplitDisplacement = 2.5;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinSplitDisplacement", m_minSplitDisplacement));
    m_minSplitDisplacementSquared = m_minSplitDisplacement * m_minSplitDisplacement;


    return STATUS_CODE_SUCCESS;
}

} // namespace lar
