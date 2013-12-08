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



    ClusterVector clusterVector;
    this->GetListOfCleanClusters(pClusterList, clusterVector);
    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNOccupiedLayers);


    typedef std::map<pandora::Cluster*,LArClusterHelper::TwoDSlidingFitResult> TwoDSlidingFitResultMap;

    

    TwoDSlidingFitResultMap slidingFitResultMap;

    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
    {
        LArClusterHelper::TwoDSlidingFitResult slidingFitResult;
        LArClusterHelper::LArTwoDSlidingFit(*iter, m_halfWindowLayers, slidingFitResult);

        if (!slidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(*iter, slidingFitResult)).second)
            throw StatusCodeException(STATUS_CODE_FAILURE);
    }


    unsigned int ctrI(0);

    for (ClusterVector::iterator iterI = clusterVector.begin(), iterEndI = clusterVector.end(); iterI != iterEndI; ++iterI)
    {
        Cluster* pBranchCluster = *iterI;   ++ctrI;

        if (NULL == pBranchCluster)
	    continue;

        TwoDSlidingFitResultMap::const_iterator iterFitI = slidingFitResultMap.find(*iterI);
        if (slidingFitResultMap.end() == iterFitI)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        const LArClusterHelper::TwoDSlidingFitResult &branchSlidingFit(iterFitI->second);

        unsigned int ctrJ(0);

        for (ClusterVector::iterator iterJ = clusterVector.begin(), iterEndJ = clusterVector.end(); iterJ != iterEndJ; ++iterJ)
        {
            Cluster* pReplacementCluster = *iterJ;   ++ctrJ;

            if (NULL == pReplacementCluster)
	        continue;

            std::cout << " [" << ctrI << "][" << ctrJ << "] " << std::endl;

            if (pBranchCluster == pReplacementCluster)
	        continue;

            TwoDSlidingFitResultMap::const_iterator iterFitJ = slidingFitResultMap.find(*iterJ);
            if (slidingFitResultMap.end() == iterFitJ)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            const LArClusterHelper::TwoDSlidingFitResult &replacementSlidingFit(iterFitJ->second);

            CartesianVector replacementStartPosition(0.f, 0.f, 0.f);
            CartesianVector branchStartPosition(0.f, 0.f, 0.f);

            if (STATUS_CODE_SUCCESS != this->FindBestSplitPosition(branchSlidingFit, replacementSlidingFit, 
                branchStartPosition, replacementStartPosition))
	        continue;

	    std::cout << "   foundSplit: [" << ctrI << "][" << ctrJ << "] " << std::endl;

            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReplaceBranch(pBranchCluster, pReplacementCluster,
	        branchStartPosition, replacementStartPosition));

            *iterI = NULL;
            *iterJ = NULL;

            break;
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

        if (1 + pCluster->GetOuterPseudoLayer() - pCluster->GetInnerPseudoLayer() < 5)
            continue;

        if (LArClusterHelper::GetLengthSquared(pCluster) < 25.f)
            continue;

        clusterVector.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode BranchSplittingAlgorithm::FindBestSplitPosition(const LArClusterHelper::TwoDSlidingFitResult &branchSlidingFit, const LArClusterHelper::TwoDSlidingFitResult &principalSlidingFit, CartesianVector& branchStartPosition, CartesianVector& principalStartPosition) const
{
    bool foundSplit(false);

    const float principalLengthSquared((principalSlidingFit.GetGlobalMaxLayerPosition() - principalSlidingFit.GetGlobalMinLayerPosition()).GetMagnitudeSquared());
    const float branchLengthSquared((branchSlidingFit.GetGlobalMaxLayerPosition() - branchSlidingFit.GetGlobalMinLayerPosition()).GetMagnitudeSquared());

    if (principalLengthSquared < m_minClusterLength * m_minClusterLength || branchLengthSquared < m_minClusterLength * m_minClusterLength)
        return STATUS_CODE_NOT_FOUND;;

    for (unsigned int principalForward = 0; principalForward < 2; ++principalForward)
    {
        const CartesianVector principalVertex(0==principalForward ? principalSlidingFit.GetGlobalMinLayerPosition() : principalSlidingFit.GetGlobalMaxLayerPosition());
        const CartesianVector principalEnd(0==principalForward ? principalSlidingFit.GetGlobalMaxLayerPosition() : principalSlidingFit.GetGlobalMinLayerPosition());
        const CartesianVector principalDirection(0==principalForward ? principalSlidingFit.GetGlobalMinLayerDirection() : principalSlidingFit.GetGlobalMaxLayerDirection() * -1.0);

	if (LArClusterHelper::GetClosestDistance(principalVertex, branchSlidingFit.GetCluster()) > 25.f)
            continue;

        for (unsigned int branchForward = 0; branchForward < 2; ++branchForward)
        {
            const CartesianVector branchVertex(0==branchForward ? branchSlidingFit.GetGlobalMinLayerPosition() : branchSlidingFit.GetGlobalMaxLayerPosition());
            const CartesianVector branchEnd(0==branchForward ? branchSlidingFit.GetGlobalMaxLayerPosition() : branchSlidingFit.GetGlobalMinLayerPosition()); 
            const CartesianVector branchDirection(0==branchForward ? branchSlidingFit.GetGlobalMinLayerDirection() : branchSlidingFit.GetGlobalMaxLayerDirection() * -1.0); 

            const float vertex_to_vertex((principalVertex - branchVertex).GetMagnitudeSquared());
            const float vertex_to_end((principalVertex - branchEnd).GetMagnitudeSquared());
            const float end_to_vertex((principalEnd - branchVertex).GetMagnitudeSquared());
            const float end_to_end((principalEnd - branchEnd).GetMagnitudeSquared());

            const float cosRelativeAngle(-branchDirection.GetDotProduct(principalDirection));

            if (vertex_to_vertex > std::min(end_to_end, std::min(vertex_to_end, end_to_vertex)))
	        continue;

            if (end_to_end < std::max(vertex_to_vertex, std::max(vertex_to_end, end_to_vertex)))
	        continue;

            if ((branchVertex - principalVertex).GetMagnitudeSquared() < m_maxTransverseDisplacement * m_maxTransverseDisplacement)
	        continue;

	    if (cosRelativeAngle > m_minCosRelativeAngle || cosRelativeAngle < 0.f)
	        continue;

// ClusterList tempList1, tempList2;
// Cluster* pBranchCluster = (Cluster*)(branchSlidingFit.GetCluster());
// Cluster* pRelacementCluster = (Cluster*)(principalSlidingFit.GetCluster());
// tempList1.insert(pBranchCluster);
// tempList2.insert(pRelacementCluster);
// PANDORA_MONITORING_API(SetEveDisplayParameters(false, false, -1, 1));
// PANDORA_MONITORING_API(VisualizeClusters(&tempList1, "BranchCluster", GREEN));
// PANDORA_MONITORING_API(VisualizeClusters(&tempList2, "PrincipalCluster", BLUE));
// PANDORA_MONITORING_API(AddMarkerToVisualization(&principalVertex, "PrincipalVertex", BLACK, 2.75));
// PANDORA_MONITORING_API(AddMarkerToVisualization(&branchVertex, "BranchVertex", BLACK, 2.75));
// PANDORA_MONITORING_API(ViewEvent());

            const float stepSize(branchSlidingFit.GetL(m_stepSizeLayers));
            const float deltaL(0==branchForward ? +branchSlidingFit.GetL(m_halfWindowLayers) : -branchSlidingFit.GetL(m_halfWindowLayers));

	    float branchDistance(0.5 * stepSize);
          
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
                    CartesianVector truncatedDirection(0.f,0.f,0.f);
                    branchSlidingFit.GetLocalPosition(linearProjection, localL, localT);
                    branchSlidingFit.GetGlobalFitPosition(localL, truncatedPosition);
                    branchSlidingFit.GetGlobalFitDirection(localL + deltaL, truncatedDirection);

// ClusterList tempList1, tempList2;
// Cluster* pBranchCluster = (Cluster*)(branchSlidingFit.GetCluster());
// Cluster* pRelacementCluster = (Cluster*)(principalSlidingFit.GetCluster());
// tempList1.insert(pBranchCluster);
// tempList2.insert(pRelacementCluster);
// PANDORA_MONITORING_API(SetEveDisplayParameters(false, false, -1, 1));
// PANDORA_MONITORING_API(VisualizeClusters(&tempList1, "BranchCluster", RED));
// PANDORA_MONITORING_API(VisualizeClusters(&tempList2, "PrincipalCluster", BLUE));
// PANDORA_MONITORING_API(AddMarkerToVisualization(&truncatedPosition, "TestPosition", BLACK, 4.0));
// for(int nCheck =-10; nCheck<=10; ++nCheck)
// {
// CartesianVector testPosition(truncatedPosition + truncatedDirection * 0.5 * (float)(nCheck+0.5));
// PANDORA_MONITORING_API(AddMarkerToVisualization(&testPosition, "TestPosition", BLACK, 2.0));
// }
// PANDORA_MONITORING_API(ViewEvent());

                    const float cosTheta(0==branchForward ? -truncatedDirection.GetDotProduct(principalDirection) : +truncatedDirection.GetDotProduct(principalDirection));

                    float rT1(0.f), rL1(0.f), rT2(0.f), rL2(0.f);
                    LArVertexHelper::GetImpactParameters(truncatedPosition, truncatedDirection, principalVertex, rL1, rT1);
                    LArVertexHelper::GetImpactParameters(principalVertex, principalDirection, truncatedPosition, rL2, rT2);

// std::cout << "   cosTheta=" << cosTheta << " rT1=" << rT1 << " rT2=" << rT2 << std::endl; 

                    if ((cosTheta > m_minCosRelativeAngle) && (rT1 < m_maxTransverseDisplacement) && (rT2 < m_maxTransverseDisplacement))
		    {
                        foundSplit = true;
                        branchStartPosition = truncatedPosition;
                        principalStartPosition = principalVertex;
		    }
		}
                catch (StatusCodeException &)
		{
		}
	    }
	}
    }
		    
if(foundSplit)
{
ClusterList tempList1, tempList2;
Cluster* pBranchCluster = (Cluster*)(branchSlidingFit.GetCluster());
Cluster* pRelacementCluster = (Cluster*)(principalSlidingFit.GetCluster());
tempList1.insert(pBranchCluster);
tempList2.insert(pRelacementCluster);
PANDORA_MONITORING_API(SetEveDisplayParameters(false, false, -1, 1));
PANDORA_MONITORING_API(VisualizeClusters(&tempList1, "BranchCluster", RED));
PANDORA_MONITORING_API(VisualizeClusters(&tempList2, "ReplacementCluster", BLUE));
PANDORA_MONITORING_API(AddMarkerToVisualization(&principalStartPosition, "ReplacementStartPosition", BLACK, 2.75));
PANDORA_MONITORING_API(AddMarkerToVisualization(&branchStartPosition, "BranchStartPosition", BLACK, 2.75));
PANDORA_MONITORING_API(ViewEvent());
}

    if (!foundSplit)
        return STATUS_CODE_NOT_FOUND;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode BranchSplittingAlgorithm::ReplaceBranch(Cluster *const pBranchCluster, Cluster *const pReplacementCluster,
    const CartesianVector &branchStartPosition, const CartesianVector &replacementStartPosition) const
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
    const CartesianVector branchProjection((replacementStartPosition - branchStartPosition).GetUnitVector());

    for (CaloHitList::const_iterator iter = caloHitsToDistribute.begin(), iterEnd = caloHitsToDistribute.end(); iter != iterEnd; ++iter)
    {
        CaloHit *pCaloHit = *iter;

        if (branchProjection.GetDotProduct((pCaloHit->GetPositionVector() - branchStartPosition)) > 0)
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
    m_halfWindowLayers = 10;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "HalfWindow", m_halfWindowLayers));

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
