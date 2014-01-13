/**
 *  @file   LArContent/src/LArClusterAssociation/EndToEndMergingAlgorithm.cc
 *
 *  @brief  Implementation of the longitudinal extension algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArClusterAssociation/EndToEndMergingAlgorithm.h"

#include "LArHelpers/LArClusterHelper.h"

using namespace pandora;

namespace lar
{

void EndToEndMergingAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster* pCluster = *iter;

        if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLength * m_minClusterLength)
	    continue;

        clusterVector.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EndToEndMergingAlgorithm::FillClusterAssociationMatrix(const ClusterVector &clusterVector, ClusterAssociationMatrix &clusterAssociationMatrix) const
{
    for (ClusterVector::const_iterator iterI = clusterVector.begin(), iterEndI = clusterVector.end(); iterI != iterEndI; ++iterI)
    {
	const Cluster *pParentCluster = *iterI;

	for (ClusterVector::const_iterator iterJ = clusterVector.begin(), iterEndJ = clusterVector.end(); iterJ != iterEndJ; ++iterJ)
	{
	    const Cluster *pDaughterCluster = *iterJ;

	    if (pParentCluster == pDaughterCluster)
		continue;

	    this->FillClusterAssociationMatrix(pParentCluster, pDaughterCluster, clusterAssociationMatrix);
	}
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EndToEndMergingAlgorithm::FillClusterAssociationMatrix(const Cluster* const pParentCluster, const Cluster* const pDaughterCluster, ClusterAssociationMatrix &clusterAssociationMatrix) const
{
    CartesianVector innerCoordinateP(0.f,0.f,0.f), outerCoordinateP(0.f,0.f,0.f);
    CartesianVector innerCoordinateD(0.f,0.f,0.f), outerCoordinateD(0.f,0.f,0.f);

    LArClusterHelper::GetExtremalCoordinatesXZ(pParentCluster, innerCoordinateP, outerCoordinateP);
    LArClusterHelper::GetExtremalCoordinatesXZ(pDaughterCluster, innerCoordinateD, outerCoordinateD);

    for (unsigned int useInnerD = 0; useInnerD < 2; ++useInnerD)
    {
	const CartesianVector daughterVertex(useInnerD == 1 ? innerCoordinateD : outerCoordinateD);
	const CartesianVector daughterEnd(useInnerD == 1 ? outerCoordinateD : innerCoordinateD);
	const CartesianVector projectedVertex(LArClusterHelper::GetClosestPosition(daughterVertex, pParentCluster));

	const float daughterVertexDistanceSquared((projectedVertex - daughterVertex).GetMagnitudeSquared());
	const float daughterEndDistanceSquared((projectedVertex - daughterEnd).GetMagnitudeSquared());
        const float daughterLengthSquared((daughterEnd - daughterVertex).GetMagnitudeSquared());

	if (daughterVertexDistanceSquared > std::min(m_maxLongitudinalDisplacement * m_maxLongitudinalDisplacement,
	    std::min(daughterEndDistanceSquared, daughterLengthSquared)))
	    continue;

        const ClusterAssociation::VertexType daughterVertexType(useInnerD == 1 ? ClusterAssociation::INNER : ClusterAssociation::OUTER);
	const float figureOfMerit(daughterVertexDistanceSquared);

        ClusterAssociation::AssociationType associationType(ClusterAssociation::WEAK);
        ClusterAssociation::VertexType parentVertexType(ClusterAssociation::UNDEFINED);

	for (unsigned int useInnerP = 0; useInnerP < 2; ++useInnerP)
	{
	    const CartesianVector parentVertex(useInnerP == 1 ? innerCoordinateP : outerCoordinateP);
	    const CartesianVector parentEnd(useInnerP == 1 ? outerCoordinateP : innerCoordinateP);

	    const float parentVertexDistanceSquared((projectedVertex - parentVertex).GetMagnitudeSquared());
	    const float parentEndDistanceSquared((projectedVertex - parentEnd).GetMagnitudeSquared());

            if (parentVertexDistanceSquared < parentEndDistanceSquared)
	        parentVertexType = (useInnerP == 1 ? ClusterAssociation::INNER : ClusterAssociation::OUTER);

	    if (parentVertexDistanceSquared > std::min(m_maxTransverseDisplacement * m_maxTransverseDisplacement,
	        std::min(parentEndDistanceSquared, daughterEndDistanceSquared)))
		continue;

	    const CartesianVector daughterDirection((daughterEnd - daughterVertex).GetUnitVector());
	    const CartesianVector parentDirection((parentEnd - projectedVertex).GetUnitVector());

            const float forwardDistance(daughterDirection.GetDotProduct((daughterVertex - projectedVertex)));
	    const float sidewaysDistance(daughterDirection.GetCrossProduct((daughterVertex - projectedVertex)).GetMagnitude());

            if (forwardDistance < 0.f || forwardDistance > m_maxLongitudinalDisplacement || sidewaysDistance > m_maxTransverseDisplacement)
	        continue;

	    if (-parentDirection.GetDotProduct(daughterDirection) < 0.25)
	        continue;

            // NEED SOMETHING ELSE HERE TO AVOID MAKING LARGE-ANGLE JOINS BETWEEN LONG TRACKS...

            associationType = ClusterAssociation::STRONG;
            break;
	}

        if (parentVertexType > ClusterAssociation::UNDEFINED)
	{
            (void) clusterAssociationMatrix[pParentCluster].insert(ClusterAssociationMap::value_type(pDaughterCluster, 
	        ClusterAssociation(parentVertexType, daughterVertexType, associationType, figureOfMerit)));
	}
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EndToEndMergingAlgorithm::FillClusterMergeMap(const ClusterAssociationMatrix &clusterAssociationMatrix, ClusterMergeMap &clusterMergeMap) const
{

}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EndToEndMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_minClusterLength = 1.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
	"MinClusterLength", m_minClusterLength));

    m_maxLongitudinalDisplacement = 2.5f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
	"MaxLongitudinalDisplacement", m_maxLongitudinalDisplacement));

    m_maxTransverseDisplacement = 1.5f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
	"MaxTransverseDisplacement", m_maxTransverseDisplacement));

    return ClusterExtensionAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
