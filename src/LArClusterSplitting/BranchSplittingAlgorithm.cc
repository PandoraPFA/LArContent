/**
 *  @file   LArContent/src/ClusterSplitting/BranchSplittingAlgorithm.cc
 * 
 *  @brief  Implementation of the branch splitting algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArClusterSplitting/BranchSplittingAlgorithm.h"

#include "LArHelpers/LArVertexHelper.h"
#include "LArHelpers/LArGeometryHelper.h"

using namespace pandora;

namespace lar
{

StatusCode BranchSplittingAlgorithm::Run()
{
    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentClusterList(*this, pClusterList));


    typedef std::map<pandora::Cluster*,LArClusterHelper::TwoDSlidingFitResult> TwoDSlidingFitResultMap;
    TwoDSlidingFitResultMap branchSlidingFitResultMap, replacementSlidingFitResultMap;


    bool carryOn(true);

    while (carryOn)
    {
        carryOn = false;
        
        ClusterVector clusterVector;
        this->GetListOfCleanClusters(pClusterList, clusterVector);
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

// ClusterList tempList1, tempList2;
// tempList1.insert(pBranchCluster);
// tempList2.insert(pReplacementCluster);
// PANDORA_MONITORING_API(SetEveDisplayParameters(false, false, -1, 1));
// PANDORA_MONITORING_API(VisualizeClusters(&tempList1, "BranchCluster", RED));
// PANDORA_MONITORING_API(VisualizeClusters(&tempList2, "ReplacementCluster", BLUE));
// PANDORA_MONITORING_API(AddMarkerToVisualization(&branchStartPosition, "BranchStartPosition", BLACK, 2.75));
// PANDORA_MONITORING_API(ViewEvent());

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

void BranchSplittingAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLength * m_minClusterLength)
            continue;

        clusterVector.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode BranchSplittingAlgorithm::FindBestSplitPosition(const LArClusterHelper::TwoDSlidingFitResult &branchSlidingFit, const LArClusterHelper::TwoDSlidingFitResult &principalSlidingFit, CartesianVector &branchStartPosition, CartesianVector &branchStartDirection) const
{
    bool foundSplit(false);

    for (unsigned int principalForward = 0; principalForward < 2; ++principalForward)
    {
        const CartesianVector principalVertex(1==principalForward ? principalSlidingFit.GetGlobalMinLayerPosition() : principalSlidingFit.GetGlobalMaxLayerPosition());
        const CartesianVector principalEnd(1==principalForward ? principalSlidingFit.GetGlobalMaxLayerPosition() : principalSlidingFit.GetGlobalMinLayerPosition());
        const CartesianVector principalDirection(1==principalForward ? principalSlidingFit.GetGlobalMinLayerDirection() : principalSlidingFit.GetGlobalMaxLayerDirection() * -1.f);

	if (LArClusterHelper::GetClosestDistance(principalVertex, branchSlidingFit.GetCluster()) > m_maxLongitudinalDisplacement)
            continue;


        for (unsigned int branchForward = 0; branchForward < 2; ++branchForward)
        {
            const CartesianVector branchVertex(1==branchForward ? branchSlidingFit.GetGlobalMinLayerPosition() : branchSlidingFit.GetGlobalMaxLayerPosition());
            const CartesianVector branchEnd(1==branchForward ? branchSlidingFit.GetGlobalMaxLayerPosition() : branchSlidingFit.GetGlobalMinLayerPosition()); 
            const CartesianVector branchDirection(1==branchForward ? branchSlidingFit.GetGlobalMinLayerDirection() : branchSlidingFit.GetGlobalMaxLayerDirection() * -1.f); 

            const float vertex_to_vertex((principalVertex - branchVertex).GetMagnitudeSquared());
            const float vertex_to_end((principalVertex - branchEnd).GetMagnitudeSquared());
            const float end_to_vertex((principalEnd - branchVertex).GetMagnitudeSquared());
            const float end_to_end((principalEnd - branchEnd).GetMagnitudeSquared());

            // sign convention for vertexProjection: positive means that clusters overlap
            const float vertexProjection(+branchDirection.GetDotProduct(principalVertex - branchVertex));
            const float cosRelativeAngle(-branchDirection.GetDotProduct(principalDirection));

            if (vertex_to_vertex > std::min(end_to_end, std::min(vertex_to_end, end_to_vertex)))
	        continue;

            if (end_to_end < std::max(vertex_to_vertex, std::max(vertex_to_end, end_to_vertex)))
	        continue;

            if (vertexProjection < 0.f && cosRelativeAngle > m_minCosRelativeAngle)
	        continue;

            if (cosRelativeAngle < 0.f)
	        continue;

// ClusterList tempList1, tempList2;
// Cluster* pBranchCluster = (Cluster*)(branchSlidingFit.GetCluster());
// Cluster* pReplacementCluster = (Cluster*)(principalSlidingFit.GetCluster());
// tempList1.insert(pBranchCluster);
// tempList2.insert(pReplacementCluster);
// PANDORA_MONITORING_API(SetEveDisplayParameters(false, false, -1, 1));
// PANDORA_MONITORING_API(VisualizeClusters(&tempList1, "BranchCluster", GREEN));
// PANDORA_MONITORING_API(VisualizeClusters(&tempList2, "PrincipalCluster", BLUE));
// PANDORA_MONITORING_API(AddMarkerToVisualization(&principalVertex, "PrincipalVertex", BLACK, 2.75));
// PANDORA_MONITORING_API(AddMarkerToVisualization(&branchVertex, "BranchVertex", BLACK, 2.75));
// PANDORA_MONITORING_API(ViewEvent());

            const float stepSize(branchSlidingFit.GetL(m_stepSizeLayers));
            const float deltaL(1==branchForward ? +branchSlidingFit.GetL(m_shortHalfWindowLayers) : -branchSlidingFit.GetL(m_shortHalfWindowLayers));

	    float branchDistance(std::max(0.f,vertexProjection) + 0.5f * stepSize);
          
	    while (!foundSplit)
	    {
	        branchDistance += stepSize;

	        const CartesianVector linearProjection(branchVertex + branchDirection * branchDistance);
 
                if (principalDirection.GetDotProduct(linearProjection - principalVertex) < -m_maxLongitudinalDisplacement)
		    break;
		
                if ((linearProjection - branchVertex).GetMagnitudeSquared() > (linearProjection - branchEnd).GetMagnitudeSquared())
		    break;
		
                try{
                    float localL(0.f), localT(0.f);
                    CartesianVector truncatedPosition(0.f,0.f,0.f);
                    CartesianVector forwardDirection(0.f,0.f,0.f);
                    branchSlidingFit.GetLocalPosition(linearProjection, localL, localT);
                    branchSlidingFit.GetGlobalFitPosition(localL, truncatedPosition);
                    branchSlidingFit.GetGlobalFitDirection(localL + deltaL, forwardDirection);

                    CartesianVector truncatedDirection(1==branchForward ? forwardDirection : forwardDirection * -1.f);
                    const float cosTheta(-truncatedDirection.GetDotProduct(principalDirection));

                    float rT1(0.f), rL1(0.f), rT2(0.f), rL2(0.f);
                    LArVertexHelper::GetImpactParameters(truncatedPosition, truncatedDirection, principalVertex, rL1, rT1);
                    LArVertexHelper::GetImpactParameters(principalVertex, principalDirection, truncatedPosition, rL2, rT2);

                    if ((cosTheta > m_minCosRelativeAngle) && (rT1 < m_maxTransverseDisplacement) && (rT2 < m_maxTransverseDisplacement))
		    {
                        foundSplit = true;
                        branchStartPosition = truncatedPosition;
                        branchStartDirection = truncatedDirection * -1.f;
		    }
		}
                catch (StatusCodeException &)
		{
		}
	    }

            if (foundSplit)
                return STATUS_CODE_SUCCESS;
        }
    }

    return STATUS_CODE_NOT_FOUND;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode BranchSplittingAlgorithm::ReplaceBranch(Cluster *const pBranchCluster, Cluster *const pReplacementCluster,
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

StatusCode BranchSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{ 
    m_shortHalfWindowLayers = 10;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShortHalfWindow", m_shortHalfWindowLayers));

    m_longHalfWindowLayers = 20;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "LongHalfWindow", m_longHalfWindowLayers));

    m_stepSizeLayers = 3;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "StepSize", m_stepSizeLayers));

    m_minClusterLength = 7.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLength", m_minClusterLength));

    m_maxTransverseDisplacement = 1.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTransverseDisplacement", m_maxTransverseDisplacement));

    m_maxLongitudinalDisplacement = 10.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxLongitudinalDisplacement", m_maxLongitudinalDisplacement));

    m_minCosRelativeAngle = 0.985;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCosRelativeAngle", m_minCosRelativeAngle));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
