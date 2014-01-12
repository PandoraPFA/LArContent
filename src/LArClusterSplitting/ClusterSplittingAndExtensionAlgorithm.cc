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


        // Loop over each possible pair of clusters
        for (ClusterVector::iterator iterI = clusterVector.begin(), iterEndI = clusterVector.end(); iterI != iterEndI; ++iterI)
        {
            Cluster* pClusterI = *iterI;

            if (NULL == pClusterI)
	        continue;

            for (ClusterVector::iterator iterJ = iterI, iterEndJ = clusterVector.end(); iterJ != iterEndJ; ++iterJ)
            {
                Cluster* pClusterJ = *iterJ;

                if (NULL == pClusterJ)
	            continue;

                if (pClusterI == pClusterJ)
		    continue;

                // Get the branch and replacement sliding fits for this pair of clusters
                TwoDSlidingFitResultMap::const_iterator iterBranchI = branchSlidingFitResultMap.find(*iterI);
                TwoDSlidingFitResultMap::const_iterator iterBranchJ = branchSlidingFitResultMap.find(*iterJ);

                TwoDSlidingFitResultMap::const_iterator iterReplacementI = replacementSlidingFitResultMap.find(*iterI);
                TwoDSlidingFitResultMap::const_iterator iterReplacementJ = replacementSlidingFitResultMap.find(*iterJ);

                if (branchSlidingFitResultMap.end() == iterBranchI || branchSlidingFitResultMap.end() == iterBranchJ ||
                    replacementSlidingFitResultMap.end() == iterReplacementI || replacementSlidingFitResultMap.end() == iterReplacementJ)
                    throw StatusCodeException(STATUS_CODE_FAILURE);

                const LArClusterHelper::TwoDSlidingFitResult &branchSlidingFitI(iterBranchI->second);
                const LArClusterHelper::TwoDSlidingFitResult &branchSlidingFitJ(iterBranchJ->second);

                const LArClusterHelper::TwoDSlidingFitResult &replacementSlidingFitI(iterReplacementI->second);
                const LArClusterHelper::TwoDSlidingFitResult &replacementSlidingFitJ(iterReplacementJ->second);

                // Search for a split in clusterI
                float branchChisqI(0.f);
                CartesianVector branchSplitPositionI(0.f, 0.f, 0.f);
                CartesianVector branchSplitDirectionI(0.f, 0.f, 0.f);

                try{
                    this->FindBestSplitPosition(branchSlidingFitI, replacementSlidingFitJ, branchSplitPositionI, branchSplitDirectionI);
                    branchChisqI = this->CalculateBranchChi2(pClusterI, branchSplitPositionI, branchSplitDirectionI);
		}
                catch (StatusCodeException &)
		{
		}

                // Search for a split in clusterJ
                float branchChisqJ(0.f);
                CartesianVector branchSplitPositionJ(0.f, 0.f, 0.f);
                CartesianVector branchSplitDirectionJ(0.f, 0.f, 0.f);

                try{
                    this->FindBestSplitPosition(branchSlidingFitJ, replacementSlidingFitI, branchSplitPositionJ, branchSplitDirectionJ);
                    branchChisqJ = this->CalculateBranchChi2(pClusterJ, branchSplitPositionJ, branchSplitDirectionJ); 
		}
                catch (StatusCodeException &)
		{
		}

                // Re-calculate chi2 values if both clusters have a split
                if (branchChisqI > 0.f && branchChisqJ > 0.f)  
		{
		    const CartesianVector relativeDirection((branchSplitPositionJ - branchSplitPositionI).GetUnitVector());

                    if (branchSplitDirectionI.GetDotProduct(relativeDirection) > 0.f && 
                        branchSplitDirectionJ.GetDotProduct(relativeDirection) < 0.f )
		    {
                        try{
                            const float newBranchChisqI(this->CalculateBranchChi2(pClusterI, branchSplitPositionI, relativeDirection));
                            const float newBranchChisqJ(this->CalculateBranchChi2(pClusterJ, branchSplitPositionJ, relativeDirection * -1.f));
                            branchChisqI = newBranchChisqI;
                            branchChisqJ = newBranchChisqJ;
			}
                        catch (StatusCodeException &)
		        {
		        }
		    }
	        }

              
                // Select the overall best split position
                Cluster* pBranchCluster = NULL;
                Cluster* pReplacementCluster = NULL;
                CartesianVector branchSplitPosition(0.f, 0.f, 0.f);
                CartesianVector branchSplitDirection(0.f, 0.f, 0.f);

                if (branchChisqI > branchChisqJ)
		{
		    pBranchCluster = pClusterI;
                    pReplacementCluster = pClusterJ;
                    branchSplitPosition = branchSplitPositionI;
                    branchSplitDirection = branchSplitDirectionI;
		}

                else if (branchChisqJ > branchChisqI)
		{
                    pBranchCluster = pClusterJ;
                    pReplacementCluster = pClusterI;
                    branchSplitPosition = branchSplitPositionJ;
                    branchSplitDirection = branchSplitDirectionJ;
		}
                
                else
		    continue;

// --- BEGIN DISPLAY ---
// ClusterList tempList1, tempList2;
// tempList1.insert(pBranchCluster);
// tempList2.insert(pReplacementCluster);
// PANDORA_MONITORING_API(SetEveDisplayParameters(false, false, -1, 1));
// PANDORA_MONITORING_API(VisualizeClusters(&tempList1, "BranchCluster", BLUE));
// PANDORA_MONITORING_API(VisualizeClusters(&tempList2, "ReplacementCluster", GREEN));
// if(branchChisqI > 0.f) 
// PANDORA_MONITORING_API(AddMarkerToVisualization(&branchSplitPositionI, "BranchSplitPositionI", RED, 2.75));
// if(branchChisqJ > 0.f) 
// PANDORA_MONITORING_API(AddMarkerToVisualization(&branchSplitPositionJ, "BranchSplitPositionJ", RED, 2.75));
// PANDORA_MONITORING_API(AddMarkerToVisualization(&branchSplitPosition, "BranchSplitPosition", BLACK, 3.5));
// PANDORA_MONITORING_API(ViewEvent());
// --- END DISPLAY ---

                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReplaceBranch(pBranchCluster, pReplacementCluster,
	            branchSplitPosition, branchSplitDirection));

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

float ClusterSplittingAndExtensionAlgorithm::CalculateBranchChi2(const Cluster* const pCluster, const CartesianVector &splitPosition, const CartesianVector &splitDirection) const
{
    CaloHitList principalCaloHitList, branchCaloHitList; 

    this->SplitBranchCluster(pCluster, splitPosition, splitDirection, principalCaloHitList, branchCaloHitList);

    float totalChi2(0.f);
    float totalHits(0.f);

    for (CaloHitList::const_iterator iter = branchCaloHitList.begin(), iterEnd = branchCaloHitList.end(); iter != iterEnd; ++iter)
    {
        CaloHit *pCaloHit = *iter;

        const CartesianVector hitPosition(pCaloHit->GetPositionVector());
        const CartesianVector projectedPosition(splitPosition + splitDirection * splitDirection.GetDotProduct(hitPosition - splitPosition));

        totalChi2 += (hitPosition - projectedPosition).GetMagnitudeSquared();
        totalHits += 1.f;
    }

    if (totalHits > 0.f)
        return std::sqrt(totalChi2/totalHits);

    throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterSplittingAndExtensionAlgorithm::SplitBranchCluster(const Cluster* const pCluster, const CartesianVector &splitPosition, const CartesianVector &splitDirection, CaloHitList &principalCaloHitList, CaloHitList &branchCaloHitList) const
{
    // Distribute hits in branch cluster between new principal and residual clusters
    CaloHitList caloHitsToDistribute;
    pCluster->GetOrderedCaloHitList().GetCaloHitList(caloHitsToDistribute);

    for (CaloHitList::const_iterator iter = caloHitsToDistribute.begin(), iterEnd = caloHitsToDistribute.end(); iter != iterEnd; ++iter)
    {
        CaloHit *pCaloHit = *iter;

        if (splitDirection.GetDotProduct((pCaloHit->GetPositionVector() - splitPosition)) > 0.f)
        {
            branchCaloHitList.insert(pCaloHit);
        }
        else
        {
            principalCaloHitList.insert(pCaloHit);
        }
    }

    if (branchCaloHitList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterSplittingAndExtensionAlgorithm::ReplaceBranch(Cluster *const pBranchCluster, Cluster *const pReplacementCluster,
    const CartesianVector &branchSplitPosition, const CartesianVector &branchSplitDirection) const
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
    CaloHitList residualCaloHitList;
    this->SplitBranchCluster(pBranchCluster, branchSplitPosition, branchSplitDirection, principalCaloHitList, residualCaloHitList);

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
