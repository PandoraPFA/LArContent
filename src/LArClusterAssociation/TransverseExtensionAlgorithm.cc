/**
 *  @file   LArContent/src/LArClusterAssociation/TransverseExtensionAlgorithm.cc
 * 
 *  @brief  Implementation of the transverse extension algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArClusterAssociation/TransverseExtensionAlgorithm.h"

#include "LArHelpers/LArVertexHelper.h"
#include "LArHelpers/LArClusterHelper.h"

using namespace pandora;

namespace lar
{

  
void TransverseExtensionAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
   for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
        clusterVector.push_back(*iter);

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);
}
 
//------------------------------------------------------------------------------------------------------------------------------------------
   
void TransverseExtensionAlgorithm::FillClusterAssociationMatrix(const ClusterVector &clusterVector, ClusterAssociationMatrix &clusterAssociationMatrix) const
{  
    // Step 1: Separate input clusters into candidate parent and daughter clusters.
    //         Convert parent clusters into pointing clusters.
    ClusterVector daughterClusterList;
    LArPointingClusterList parentClusterList;

    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
    {
        Cluster* const pCluster(*iter);

        if (LArClusterHelper::GetLengthSquared(pCluster) > m_maxClusterLength * m_maxClusterLength )
            parentClusterList.push_back(LArPointingCluster(pCluster));
        else
	    daughterClusterList.push_back(pCluster);	
    }

  
    // Step 2: Form associations between parent and daughter clusters
    for (LArPointingClusterList::const_iterator iter1 = parentClusterList.begin(), iterEnd1 = parentClusterList.end(); iter1 != iterEnd1; ++iter1)
    {
        const LArPointingCluster &parentPointingCluster = *iter1;

        // First, associate pointing clusters with each other
        for (LArPointingClusterList::const_iterator iter2 = parentClusterList.begin(), iterEnd2 = parentClusterList.end(); iter2 != iterEnd2; ++iter2)
        {
            const LArPointingCluster &daughterPointingCluster = *iter2;

            if (parentPointingCluster.GetCluster() == daughterPointingCluster.GetCluster())
                continue;

            this->FillAssociationMatrix(parentPointingCluster, daughterPointingCluster, clusterAssociationMatrix);
        }

        // Next, associate pointing clusters with candidate daughter clusters
        for (ClusterVector::const_iterator iter2 = daughterClusterList.begin(), iterEnd2 = daughterClusterList.end(); iter2 != iterEnd2; ++iter2)
	{
	    const Cluster* pDaughterCluster = *iter2;

            if (parentPointingCluster.GetCluster() == pDaughterCluster)
	        continue;

            this->FillAssociationMatrix(parentPointingCluster, pDaughterCluster, clusterAssociationMatrix);
	}
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseExtensionAlgorithm::FillAssociationMatrix(const LArPointingCluster &clusterI, const LArPointingCluster &clusterJ, ClusterAssociationMatrix &clusterAssociationMatrix) const
{
    const Cluster *const pClusterI(clusterI.GetCluster());
    const Cluster *const pClusterJ(clusterJ.GetCluster());

    if (pClusterI == pClusterJ)
        return;

    for (unsigned int useInnerI=0; useInnerI<2; ++useInnerI)
    {
        for (unsigned int useInnerJ=0; useInnerJ<2; ++useInnerJ)
	{
            
            const LArPointingCluster::Vertex &targetVertexI(useInnerI==1 ? clusterI.GetInnerVertex() : clusterI.GetOuterVertex());
            const LArPointingCluster::Vertex &oppositeVertexI(useInnerI==1 ? clusterI.GetOuterVertex() : clusterI.GetInnerVertex());

            const LArPointingCluster::Vertex &targetVertexJ(useInnerJ==1 ? clusterJ.GetInnerVertex() : clusterJ.GetOuterVertex());
            const LArPointingCluster::Vertex &oppositeVertexJ(useInnerJ==1 ? clusterJ.GetOuterVertex() : clusterJ.GetInnerVertex());

            const float distSquared_targetI_to_targetJ((targetVertexI.GetPosition() - targetVertexJ.GetPosition()).GetMagnitudeSquared());    
            const float distSquared_targetI_to_oppositeJ((targetVertexI.GetPosition() - oppositeVertexJ.GetPosition()).GetMagnitudeSquared());

            const float distSquared_oppositeI_to_targetJ((oppositeVertexI.GetPosition() - targetVertexJ.GetPosition()).GetMagnitudeSquared());
            const float distSquared_oppositeI_to_oppositeJ((oppositeVertexI.GetPosition() - oppositeVertexJ.GetPosition()).GetMagnitudeSquared());

            const float distSquared_targetI_to_oppositeI((targetVertexI.GetPosition() - oppositeVertexI.GetPosition()).GetMagnitudeSquared());
            const float distSquared_targetJ_to_oppositeJ((targetVertexJ.GetPosition() - oppositeVertexJ.GetPosition()).GetMagnitudeSquared());

            if (distSquared_oppositeI_to_targetJ < distSquared_targetI_to_targetJ ||
                distSquared_oppositeI_to_targetJ < distSquared_targetI_to_oppositeI ||
                distSquared_targetI_to_oppositeJ < distSquared_targetI_to_targetJ ||
                distSquared_targetI_to_oppositeJ < distSquared_targetJ_to_oppositeJ)    
	        continue;  

            const CartesianVector &vertexPositionI(targetVertexI.GetPosition());
            const CartesianVector &vertexPositionJ(targetVertexJ.GetPosition());
            const CartesianVector &vertexDirectionI(targetVertexI.GetDirection());
 
            float rT(0.f), rL(0.f);
            LArVertexHelper::GetImpactParameters(vertexPositionI, vertexDirectionI, vertexPositionJ, rL, rT);
            
            if (rL > -2.5f && rL < m_maxLongitudinalDisplacement && rT < 3.f * m_maxTransverseDisplacement)
	    {

	    }
                  
	}
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseExtensionAlgorithm::FillAssociationMatrix(const LArPointingCluster &parentCluster, const Cluster* const pDaughterCluster, ClusterAssociationMatrix &clusterAssociationMatrix) const
{

}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseExtensionAlgorithm::FillClusterMergeMap(const ClusterAssociationMatrix &clusterAssociationMatrix, ClusterMergeMap &clusterMergeMap) const
{
  std::cout << " *** TransverseExtensionAlgorithm::FillClusterMergeMap(...) *** " << std::endl;

}
//------------------------------------------------------------------------------------------------------------------------------------------

// bool TransverseExtensionAlgorithm::IsEndAssociated(const LArClusterHelper::TwoDSlidingFitResult &slidingFitResult, const Cluster *const pCluster) const
// {

//     CartesianVector firstCoordinate(0.f, 0.f, 0.f);
//     CartesianVector secondCoordinate(0.f, 0.f, 0.f);
//     LArClusterHelper::GetExtremalCoordinatesXZ(pCluster, firstCoordinate, secondCoordinate);

//     for (unsigned int iForward = 0; iForward < 2; ++iForward)
//     {
//         const CartesianVector vertex(1==iForward ? slidingFitResult.GetGlobalMinLayerPosition() : slidingFitResult.GetGlobalMaxLayerPosition());
//         const CartesianVector direction(1==iForward ? slidingFitResult.GetGlobalMinLayerDirection() : slidingFitResult.GetGlobalMaxLayerDirection() * -1.f);

//         float firstL(0.f), firstT(0.f), secondT(0.f), secondL(0.f);
//         LArVertexHelper::GetImpactParameters(vertex, direction, firstCoordinate, firstL, firstT);
//         LArVertexHelper::GetImpactParameters(vertex, direction, secondCoordinate, secondL, secondT);

//         const float innerL(firstL < secondL ? firstL : secondL);
//         const float innerT(firstL < secondL ? firstT : secondT);
//         const float outerT(firstL > secondL ? firstT : secondT);

//         if ( innerL > 0.f && innerL < m_maxLongitudinalDisplacement && 
//              innerT < m_maxTransverseDisplacement && outerT < 1.5f * m_maxTransverseDisplacement )
// 	    return true;
//     }

//     return false;
// }



//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TransverseExtensionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_maxClusterLength = 5.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxClusterLength", m_maxClusterLength));

    m_maxLongitudinalDisplacement = 10.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxLongitudinalDisplacement", m_maxLongitudinalDisplacement));

    m_maxTransverseDisplacement = 1.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTransverseDisplacement", m_maxTransverseDisplacement));

    return ClusterExtensionAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
