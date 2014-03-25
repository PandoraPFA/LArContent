/**
 *  @file   LArContent/src/LArTwoDReco/ClusterSplitting/DeltaRaySplittingAlgorithm.cc
 *
 *  @brief  Implementation of the branch splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArPointingClusterHelper.h"

#include "LArTwoDReco/LArClusterSplitting/DeltaRaySplittingAlgorithm.h"

using namespace pandora;

namespace lar
{

void DeltaRaySplittingAlgorithm::FindBestSplitPosition(const TwoDSlidingFitResult &branchSlidingFit, const TwoDSlidingFitResult &principalSlidingFit, 
    CartesianVector &principalStartPosition, CartesianVector &branchSplitPosition, CartesianVector &branchSplitDirection) const
{
    // Conventions:
    // (1) Delta ray is split from the branch cluster
    // (2) Delta ray occurs where the vertex of the principal cluster meets the vertex of the branch cluster
    //     Method loops over the inner and outer positions of the principal and branch clusters, trying all
    //     possible assignments of vertex and end position until a split is found

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

            // Require vertices to be closest two ends
            const float vertex_to_vertex((principalVertex - branchVertex).GetMagnitudeSquared());
            const float vertex_to_end((principalVertex - branchEnd).GetMagnitudeSquared());
            const float end_to_vertex((principalEnd - branchVertex).GetMagnitudeSquared());
            const float end_to_end((principalEnd - branchEnd).GetMagnitudeSquared());

            // (sign convention for vertexProjection: positive means that clusters overlap)
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

            // Serach for a split by winding back the branch cluster sliding fit
            bool foundSplit(false);

            const float halfWindowLength(branchSlidingFit.GetLayerFitHalfWindowLength());
            const float deltaL(1==branchForward ? +halfWindowLength : -halfWindowLength);

            float branchDistance(std::max(0.f,vertexProjection) + 0.5f * m_stepSize);

            while (!foundSplit)
            {
                branchDistance += m_stepSize;

                const CartesianVector linearProjection(branchVertex + branchDirection * branchDistance);

                if (principalDirection.GetDotProduct(linearProjection - principalVertex) < -m_maxLongitudinalDisplacement)
                    break;

                if ((linearProjection - branchVertex).GetMagnitudeSquared() > (linearProjection - branchEnd).GetMagnitudeSquared())
                    break;

                try
                {
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
                        principalStartPosition = principalVertex;
                        branchSplitPosition = truncatedPosition;
                        branchSplitDirection = truncatedDirection * -1.f;
                    }
                }
                catch (StatusCodeException &)
                {
                }
            }

            if (foundSplit)
                return;
        }
    }

    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DeltaRaySplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_stepSize = 1.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "StepSize", m_stepSize));

    m_maxTransverseDisplacement = 1.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTransverseDisplacement", m_maxTransverseDisplacement));

    m_maxLongitudinalDisplacement = 10.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxLongitudinalDisplacement", m_maxLongitudinalDisplacement));

    m_minCosRelativeAngle = 0.985;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCosRelativeAngle", m_minCosRelativeAngle));

    return TwoDSlidingFitSplittingAndSplicingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
