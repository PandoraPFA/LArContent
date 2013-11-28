/**
 *  @file   LArContent/src/ClusterSplitting/BranchSplittingAlgorithm.cc
 * 
 *  @brief  Implementation of the branch splitting algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArParticleIdHelper.h"
#include "LArHelpers/LArVertexHelper.h"

#include "LArClusterSplitting/BranchSplittingAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode BranchSplittingAlgorithm::Run()
{
    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentClusterList(*this, pClusterList));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode BranchSplittingAlgorithm::ReplaceBranch(Cluster *const pBranchCluster, Cluster *const pReplacementCluster,
    const CartesianVector &branchStartPosition, const CartesianVector &replacementStartPosition) const
{
    // Begin cluster fragmentation operations
    ClusterList clusterList;
    clusterList.insert(pBranchCluster);
    clusterList.insert(pReplacementCluster);
    std::string clusterListToSaveName, clusterListToDeleteName;

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*this, clusterList, clusterListToDeleteName,
        clusterListToSaveName));

    // Create new clusters
    Cluster *pPrincipalCluster(NULL), *pResidualCluster(NULL);


    // Entire replacement cluster goes into principal cluster
    const OrderedCaloHitList &replacementClusterCaloHitList(pReplacementCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = replacementClusterCaloHitList.begin(); iter != replacementClusterCaloHitList.end(); ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            CaloHit *pCaloHit = *hitIter;

            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddCaloHitToCluster(*this, pPrincipalCluster, pCaloHit));
	}
    }

    // Split the branch cluster between the principal and residual cluster
    const CartesianVector branchProjection((replacementStartPosition - branchStartPosition).GetUnitVector());

    const OrderedCaloHitList &branchClusterCaloHitList(pBranchCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = branchClusterCaloHitList.begin(); iter != branchClusterCaloHitList.end(); ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            CaloHit *pCaloHit = *hitIter;

            if (branchProjection.GetDotProduct((pCaloHit->GetPositionVector() - branchStartPosition)) > 0)
	    {
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddCaloHitToCluster(*this, pResidualCluster, pCaloHit));
	    }
            else{
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddCaloHitToCluster(*this, pPrincipalCluster, pCaloHit));
	    }
	}
    }

    //
    // Delete pBranchCluster and pReplacementCluster
    //
    // Keep pPrincipalCluster and pResidualCluster
    //

    // End cluster fragmentation operations
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, clusterListToSaveName, clusterListToDeleteName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode BranchSplittingAlgorithm::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar
