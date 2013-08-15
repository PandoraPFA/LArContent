/**
 *  @file   LArContent/src/ClusterSplitting/ClusterSplittingAlgorithm.cc
 * 
 *  @brief  Implementation of the particle seed algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArParticleIdHelper.h"
#include "LArHelpers/LArVertexHelper.h"

#include "LArClusterSplitting/ClusterSplittingAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode ClusterSplittingAlgorithm::Run()
{
    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentClusterList(*this, pClusterList));

    ClusterSplittingList internalClusterList(pClusterList->begin(), pClusterList->end());

    for (ClusterSplittingList::iterator iter = internalClusterList.begin(); iter != internalClusterList.end(); ++iter)
    {
        Cluster* pCluster = *iter;

        if (!this->IsPossibleSplit(pCluster))
            continue;

        unsigned int splitLayer(std::numeric_limits<unsigned int>::max());

        if (STATUS_CODE_SUCCESS != this->FindBestSplitLayer(pCluster,splitLayer))
            continue;

        if ((splitLayer <= pCluster->GetInnerPseudoLayer()) || (splitLayer >= pCluster->GetOuterPseudoLayer()))
            continue;

// const CartesianVector& bestPosition = pCluster->GetCentroid(splitLayer);
// Cluster* tempCluster = (Cluster*)(pCluster);
// ClusterList tempList;
// tempList.insert(tempCluster);
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&tempList, "Cluster", GREEN);
// PandoraMonitoringApi::AddMarkerToVisualization(&bestPosition, "Split", RED, 1.75);
// PandoraMonitoringApi::ViewEvent();
        ClusterSplittingList clusterSplittingList;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->SplitCluster(pCluster, splitLayer, clusterSplittingList));

// ClusterList tempList(clusterSplittingList.begin(),clusterSplittingList.end());
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&tempList, "SplitCluster", AUTOITER);
// PandoraMonitoringApi::ViewEvent();
        internalClusterList.splice(internalClusterList.end(), clusterSplittingList);
        *iter = NULL;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterSplittingAlgorithm::SplitCluster(Cluster *const pCluster, const unsigned int splitLayer, ClusterSplittingList &clusterSplittingList) const
{
    // Begin cluster fragmentation operations
    ClusterList clusterList;
    clusterList.insert(pCluster);
    std::string clusterListToSaveName, clusterListToDeleteName;

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*this, clusterList, clusterListToDeleteName,
        clusterListToSaveName));

    // Create new clusters
    Cluster *pCluster1(NULL), *pCluster2(NULL);
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(); iter != orderedCaloHitList.end(); ++iter)
    {
        const unsigned int thisLayer(iter->first);

        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            CaloHit *pCaloHit = *hitIter;
            Cluster *pClusterToModify((thisLayer < splitLayer) ? pCluster1 : pCluster2);

            if (NULL == pClusterToModify)
            {
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, pCaloHit, pClusterToModify));
                clusterSplittingList.push_back(pClusterToModify);
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

StatusCode ClusterSplittingAlgorithm::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar
