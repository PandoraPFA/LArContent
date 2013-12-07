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



    LArClusterHelper::TwoDSlidingFitResultMap slidingFitResultMap;

    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
    {
        LArClusterHelper::TwoDSlidingFitResult slidingFitResult;
        LArClusterHelper::LArTwoDSlidingFit(*iter, 20, slidingFitResult);

        if (!slidingFitResultMap.insert(LArClusterHelper::TwoDSlidingFitResultMap::value_type(*iter, slidingFitResult)).second)
            throw StatusCodeException(STATUS_CODE_FAILURE);
    }


    for (ClusterVector::const_iterator iterI = clusterVector.begin(), iterEndI = clusterVector.end(); iterI != iterEndI; ++iterI)
    {
        Cluster* pBranchCluster = *iterI; 

        LArClusterHelper::TwoDSlidingFitResultMap::const_iterator iterFitI = slidingFitResultMap.find(*iterI);
        if (slidingFitResultMap.end() == iterFitI)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        const LArClusterHelper::TwoDSlidingFitResult &branchSlidingFit(iterFitI->second);

        for (ClusterVector::const_iterator iterJ = clusterVector.begin(), iterEndJ = clusterVector.end(); iterJ != iterEndJ; ++iterJ)
        {
            Cluster* pReplacementCluster = *iterJ;

            if (pBranchCluster == pReplacementCluster)
	        continue;

            LArClusterHelper::TwoDSlidingFitResultMap::const_iterator iterFitJ = slidingFitResultMap.find(*iterJ);
            if (slidingFitResultMap.end() == iterFitJ)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            const LArClusterHelper::TwoDSlidingFitResult &replacementSlidingFit(iterFitJ->second);

            CartesianVector replacementStartPosition(0.f, 0.f, 0.f);
            CartesianVector branchStartPosition(0.f, 0.f, 0.f);

            this->FindBestBranchSplitPosition(branchSlidingFit, replacementSlidingFit, branchStartPosition, replacementStartPosition);
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

void BranchSplittingAlgorithm::FindBestBranchSplitPosition(const LArClusterHelper::TwoDSlidingFitResult &branchSlidingFit, const LArClusterHelper::TwoDSlidingFitResult &principalSlidingFit, CartesianVector& branchStartPosition, CartesianVector& principalStartPosition) const
{

    const float principalLengthSquared((principalSlidingFit.GetGlobalMaxLayerPosition() - principalSlidingFit.GetGlobalMinLayerPosition()).GetMagnitudeSquared());
    const float branchLengthSquared((branchSlidingFit.GetGlobalMaxLayerPosition() - branchSlidingFit.GetGlobalMinLayerPosition()).GetMagnitudeSquared());

    if ( principalLengthSquared < 75.f || branchLengthSquared < 75.f )
        return;

    for( unsigned int principalForward = 0; principalForward < 2; ++principalForward )
    {
        const CartesianVector principalVertex(0==principalForward ? principalSlidingFit.GetGlobalMinLayerPosition() : principalSlidingFit.GetGlobalMaxLayerPosition());
        const CartesianVector principalEnd(0==principalForward ? principalSlidingFit.GetGlobalMaxLayerPosition() : principalSlidingFit.GetGlobalMinLayerPosition());
        const CartesianVector principalDirection(0==principalForward ? principalSlidingFit.GetGlobalMinLayerDirection() : principalSlidingFit.GetGlobalMaxLayerDirection() * -1.0);

	if (LArClusterHelper::GetClosestDistance(principalVertex, branchSlidingFit.GetCluster()) > 25.f)
            continue;

        for( unsigned int branchForward = 0; branchForward < 2; ++branchForward )
        {
            const CartesianVector branchVertex(0==branchForward ? branchSlidingFit.GetGlobalMinLayerPosition() : branchSlidingFit.GetGlobalMaxLayerPosition());
            const CartesianVector branchEnd(0==branchForward ? branchSlidingFit.GetGlobalMaxLayerPosition() : branchSlidingFit.GetGlobalMinLayerPosition()); 
            const CartesianVector branchDirection(0==branchForward ? branchSlidingFit.GetGlobalMinLayerDirection() : principalSlidingFit.GetGlobalMaxLayerDirection() * -1.0); 
        
            const float vertex_to_vertex((principalVertex - branchVertex).GetMagnitudeSquared());
            const float vertex_to_end((principalVertex - branchEnd).GetMagnitudeSquared());
            const float end_to_vertex((principalEnd - branchVertex).GetMagnitudeSquared());
            const float end_to_end((principalEnd - branchEnd).GetMagnitudeSquared());

            if (end_to_end < std::max(vertex_to_vertex, std::max(vertex_to_end, end_to_vertex)))
	        continue;


            const CartesianVector linearPosition(branchVertex + branchDirection * 7.5f);

            CartesianVector truncatedPosition(0.f,0.f,0.f);
            CartesianVector truncatedDirection(0.f,0.f,0.f);

            float localL(0.f), localT(0.f);
            CartesianVector projectedPosition(0.f,0.f,0.f);
            branchSlidingFit.GetGlobalFitProjection(linearPosition, projectedPosition);
            branchSlidingFit.GetLocalPosition(projectedPosition, localL, localT);
            branchSlidingFit.GetGlobalFitPosition(localL,truncatedPosition);
            branchSlidingFit.GetGlobalFitDirection(localL,truncatedDirection);
            

            const float cosTheta(0==branchForward ? -truncatedDirection.GetDotProduct(principalDirection) : +truncatedDirection.GetDotProduct(principalDirection));

            float rT1(0.f), rL1(0.f), rT2(0.f), rL2(0.f);
            LArVertexHelper::GetImpactParameters(truncatedPosition, truncatedDirection, principalVertex, rL1, rT1);
            LArVertexHelper::GetImpactParameters(principalVertex, principalDirection, truncatedPosition, rL2, rT2);

	    std::cout << " check: cosTheta=" << cosTheta << " rL1" << rL1 << " rT1=" << rT1 << " rL2=" << rL2 << " rT2=" << rT2 << std::endl; 

    
	    
// OLD SELECTION CRITERIA
// if ((cosTheta > 0.985f) && (truncatedVertexI.GetRms() < 0.5f) && (std::fabs(rL1) < 10.f) && (std::fabs(rL2) < 10.f) &&
//     (rT1 < m_spatialResolution) && (rT2 < m_spatialResolution))

ClusterList tempList1, tempList2;
Cluster* pBranchCluster = (Cluster*)(branchSlidingFit.GetCluster());
Cluster* pRelacementCluster = (Cluster*)(principalSlidingFit.GetCluster());
tempList1.insert(pBranchCluster);
tempList2.insert(pRelacementCluster);
PANDORA_MONITORING_API(SetEveDisplayParameters(false, false, -1, 1));
PANDORA_MONITORING_API(VisualizeClusters(&tempList1, "BranchCluster", RED));
PANDORA_MONITORING_API(VisualizeClusters(&tempList2, "PrincipalCluster", BLUE));
PANDORA_MONITORING_API(AddMarkerToVisualization(&principalVertex, "ReplacementStartPosition", BLACK, 2.75));
PANDORA_MONITORING_API(AddMarkerToVisualization(&truncatedPosition, "BranchStartPosition", BLACK, 2.75));
PANDORA_MONITORING_API(ViewEvent());


        }
    }


    // HANDLE DELTA RAYS
    // bool m_handleDeltaRays(true);

    // if (m_handleDeltaRays && (clusterLengthI > 10.f) && (clusterLengthJ > 10.f))
    // {
    //     // Remove 10 layers and re-calculate impact parameters
    //     if (strengthType < LongitudinalAssociation::STRONG)
    //     {
    //         const LArPointingCluster::Vertex truncatedVertexI(targetVertexI.GetCluster(), targetVertexI.IsInnerVertex(), 10);
    //         const float cosTheta(-truncatedVertexI.GetDirection().GetDotProduct(vertexDirectionJ));

    //         float rT1(0.f), rL1(0.f), rT2(0.f), rL2(0.f);
    //         LArVertexHelper::GetImpactParameters(truncatedVertexI.GetPosition(), truncatedVertexI.GetDirection(), vertexPositionJ, rL1, rT1);
    //         LArVertexHelper::GetImpactParameters(vertexPositionJ, vertexDirectionJ, truncatedVertexI.GetPosition(), rL2, rT2);

    //         if ((cosTheta > 0.985f) && (truncatedVertexI.GetRms() < 0.5f) && (std::fabs(rL1) < 10.f) && (std::fabs(rL2) < 10.f) &&
    //             (rT1 < m_spatialResolution) && (rT2 < m_spatialResolution))
    //         {
    //             associationType = LongitudinalAssociation::EMISSION;
    //             strengthType = LongitudinalAssociation::STRONG;
    //         }
    //     }

    //     if (strengthType < LongitudinalAssociation::STRONG)
    //     {
    //         const LArPointingCluster::Vertex truncatedVertexJ(targetVertexJ.GetCluster(), targetVertexJ.IsInnerVertex(), 10);
    //         const float cosTheta(-truncatedVertexJ.GetDirection().GetDotProduct(vertexDirectionI));

    //         float rT1(0.f), rL1(0.f), rT2(0.f), rL2(0.f);
    //         LArVertexHelper::GetImpactParameters(truncatedVertexJ.GetPosition(), truncatedVertexJ.GetDirection(), vertexPositionI, rL1, rT1);
    //         LArVertexHelper::GetImpactParameters(vertexPositionI, vertexDirectionI, truncatedVertexJ.GetPosition(), rL2, rT2);

    //         if ((cosTheta > 0.985f) && (truncatedVertexJ.GetRms() < 0.5f) && (std::fabs(rL1) < 10.f) && (std::fabs(rL2) < 10.f) &&
    //             (rT1 < m_spatialResolution) && (rT2 < m_spatialResolution))
    //         {
    //             associationType = LongitudinalAssociation::EMISSION;
    //             strengthType = LongitudinalAssociation::STRONG;
    //         }
    //     }
    // }


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

StatusCode BranchSplittingAlgorithm::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar
