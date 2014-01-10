/**
 *  @file   LArContent/src/ClusterSplitting/DeltaRaySplittingAlgorithm.cc
 * 
 *  @brief  Implementation of the branch splitting algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArClusterSplitting/DeltaRaySplittingAlgorithm.h"

#include "LArHelpers/LArPointingClusterHelper.h"

using namespace pandora;

namespace lar
{

StatusCode DeltaRaySplittingAlgorithm::FindBestSplitPosition(const LArClusterHelper::TwoDSlidingFitResult &branchSlidingFit, const LArClusterHelper::TwoDSlidingFitResult &principalSlidingFit, CartesianVector &branchStartPosition, CartesianVector &branchStartDirection) const
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

            const unsigned int halfWindowLayers(branchSlidingFit.GetLayerFitHalfWindow());
            const float stepSize(branchSlidingFit.GetL(m_stepSizeLayers));
            const float deltaL(1==branchForward ? +branchSlidingFit.GetL(halfWindowLayers) : -branchSlidingFit.GetL(halfWindowLayers));

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
                    LArPointingClusterHelper::GetImpactParameters(truncatedPosition, truncatedDirection, principalVertex, rL1, rT1);
                    LArPointingClusterHelper::GetImpactParameters(principalVertex, principalDirection, truncatedPosition, rL2, rT2);

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

StatusCode DeltaRaySplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{ 
    m_stepSizeLayers = 3;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "StepSize", m_stepSizeLayers));

    m_maxTransverseDisplacement = 1.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTransverseDisplacement", m_maxTransverseDisplacement));

    m_maxLongitudinalDisplacement = 10.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxLongitudinalDisplacement", m_maxLongitudinalDisplacement));

    m_minCosRelativeAngle = 0.985;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCosRelativeAngle", m_minCosRelativeAngle));

    return BranchSplittingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
