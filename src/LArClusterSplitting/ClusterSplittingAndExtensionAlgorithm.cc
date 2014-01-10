/**
 *  @file   LArContent/src/ClusterSplitting/ClusterSplittingAndExtensionAlgorithm.cc
 * 
 *  @brief  Implementation of the branch splitting algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArClusterSplitting/ClusterSplittingAndExtensionAlgorithm.h"

#include "LArHelpers/LArPointingClusterHelper.h"

using namespace pandora;

namespace lar
{

StatusCode ClusterSplittingAndExtensionAlgorithm::Run()
{
    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentClusterList(*this, pClusterList));

    typedef std::map<pandora::Cluster*,LArClusterHelper::TwoDSlidingFitResult> TwoDSlidingFitResultMap;
    TwoDSlidingFitResultMap branchSlidingFitResultMap, replacementSlidingFitResultMap;

    bool carryOn(true);

    while (carryOn)
    {
        carryOn = false;
        
        // Get ordered list of candidate clusters
        ClusterVector clusterVector;

        for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
        {
            Cluster *pCluster = *iter;

            if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLength * m_minClusterLength)
                continue;

            clusterVector.push_back(pCluster);
	}

        std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);


        // Calculate sliding fit results 
        for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
        {
	    // Possible branch clusters (use a soft sliding fit for these)
	    if (branchSlidingFitResultMap.end() == branchSlidingFitResultMap.find(*iter))
	    {                
                LArClusterHelper::TwoDSlidingFitResult branchSlidingFitResult;
                LArClusterHelper::LArTwoDSlidingFit(*iter, m_shortHalfWindowLayers, branchSlidingFitResult);

                if (!branchSlidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(*iter, branchSlidingFitResult)).second)
                    throw StatusCodeException(STATUS_CODE_FAILURE);
	    }

            // Possible replacement clusters (use a hard linear fit for these)
	    if (replacementSlidingFitResultMap.end() == replacementSlidingFitResultMap.find(*iter))
	    {
                LArClusterHelper::TwoDSlidingFitResult replacementSlidingFitResult;
                LArClusterHelper::LArTwoDSlidingFit(*iter, m_longHalfWindowLayers, replacementSlidingFitResult);

                if (!replacementSlidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(*iter, replacementSlidingFitResult)).second)
                    throw StatusCodeException(STATUS_CODE_FAILURE);
	    }
	}


        // First loop over possible branch clusters
        for (ClusterVector::iterator iterI = clusterVector.begin(), iterEndI = clusterVector.end(); iterI != iterEndI; ++iterI)
        {
            Cluster* pBranchCluster = *iterI;

            if (NULL == pBranchCluster)
	        continue;

            TwoDSlidingFitResultMap::const_iterator iterBranch = branchSlidingFitResultMap.find(*iterI);
            if (branchSlidingFitResultMap.end() == iterBranch)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            const LArClusterHelper::TwoDSlidingFitResult &branchSlidingFit(iterBranch->second);

            // Second loop over possible replacement clusters
            for (ClusterVector::iterator iterJ = clusterVector.begin(), iterEndJ = clusterVector.end(); iterJ != iterEndJ; ++iterJ)
            {
                Cluster* pReplacementCluster = *iterJ;

                if (NULL == pReplacementCluster)
	            continue;

                if (pBranchCluster == pReplacementCluster)
	            continue;

                TwoDSlidingFitResultMap::const_iterator iterReplacement = replacementSlidingFitResultMap.find(*iterJ);
                if (replacementSlidingFitResultMap.end() == iterReplacement)
                    throw StatusCodeException(STATUS_CODE_FAILURE);

                const LArClusterHelper::TwoDSlidingFitResult &replacementSlidingFit(iterReplacement->second);

                CartesianVector branchStartPosition(0.f, 0.f, 0.f);
                CartesianVector branchStartDirection(0.f, 0.f, 0.f);

                if (STATUS_CODE_SUCCESS != this->FindBestSplitPosition(branchSlidingFit, replacementSlidingFit, 
                    branchStartPosition, branchStartDirection))
	            continue;

// --- BEGIN DISPLAY ---
// ClusterList tempList1, tempList2;
// tempList1.insert(pBranchCluster);
// tempList2.insert(pReplacementCluster);
// PANDORA_MONITORING_API(SetEveDisplayParameters(false, false, -1, 1));
// PANDORA_MONITORING_API(VisualizeClusters(&tempList1, "BranchCluster", RED));
// PANDORA_MONITORING_API(VisualizeClusters(&tempList2, "ReplacementCluster", BLUE));
// PANDORA_MONITORING_API(AddMarkerToVisualization(&branchStartPosition, "BranchStartPosition", BLACK, 2.75));
// PANDORA_MONITORING_API(ViewEvent());
// --- END DISPLAY ---

                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReplaceBranch(pBranchCluster, pReplacementCluster,
	            branchStartPosition, branchStartDirection));

                branchSlidingFitResultMap.erase(*iterI);
                branchSlidingFitResultMap.erase(*iterJ);

                replacementSlidingFitResultMap.erase(*iterI);
                replacementSlidingFitResultMap.erase(*iterJ);

                *iterI = NULL;
                *iterJ = NULL;
                carryOn = true;

                break;
	    }
	}
    }
    
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterSplittingAndExtensionAlgorithm::ReplaceBranch(Cluster *const pBranchCluster, Cluster *const pReplacementCluster,
    const CartesianVector &branchStartPosition, const CartesianVector &branchStartDirection) const
{
    ClusterList clusterList;
    clusterList.insert(pBranchCluster);
    clusterList.insert(pReplacementCluster);

    std::string clusterListToSaveName, clusterListToDeleteName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*this, clusterList, clusterListToDeleteName,
        clusterListToSaveName));

    // Entire replacement cluster goes into new principal cluster
    CaloHitList principalCaloHitList;
    pReplacementCluster->GetOrderedCaloHitList().GetCaloHitList(principalCaloHitList);

    // Distribute hits in branch cluster between new principal and residual clusters
    CaloHitList caloHitsToDistribute;
    pBranchCluster->GetOrderedCaloHitList().GetCaloHitList(caloHitsToDistribute);

    CaloHitList residualCaloHitList;

    for (CaloHitList::const_iterator iter = caloHitsToDistribute.begin(), iterEnd = caloHitsToDistribute.end(); iter != iterEnd; ++iter)
    {
        CaloHit *pCaloHit = *iter;

        if (branchStartDirection.GetDotProduct((pCaloHit->GetPositionVector() - branchStartPosition)) > 0)
        {
            residualCaloHitList.insert(pCaloHit);
        }
        else
        {
            principalCaloHitList.insert(pCaloHit);
        }
    }

    Cluster *pPrincipalCluster(NULL), *pResidualCluster(NULL);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, &principalCaloHitList, pPrincipalCluster));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, &residualCaloHitList, pResidualCluster));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, clusterListToSaveName, clusterListToDeleteName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterSplittingAndExtensionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{ 
    m_shortHalfWindowLayers = 10;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShortHalfWindow", m_shortHalfWindowLayers));

    m_longHalfWindowLayers = 20;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "LongHalfWindow", m_longHalfWindowLayers));

    m_minClusterLength = 7.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLength", m_minClusterLength));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
