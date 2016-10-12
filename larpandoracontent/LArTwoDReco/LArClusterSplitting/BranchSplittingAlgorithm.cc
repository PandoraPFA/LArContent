/**
 *  @file   larpandoracontent/LArTwoDReco/ClusterSplitting/BranchSplittingAlgorithm.cc
 *
 *  @brief  Implementation of the branch splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/BranchSplittingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

BranchSplittingAlgorithm::BranchSplittingAlgorithm() :
    m_maxTransverseDisplacement(1.5f),
    m_maxLongitudinalDisplacement(10.f),
    m_minLongitudinalExtension(3.f),
    m_minCosRelativeAngle(0.966f),
    m_projectionAngularAllowance(20.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BranchSplittingAlgorithm::FindBestSplitPosition(const TwoDSlidingFitResult &branchSlidingFit, const TwoDSlidingFitResult &principalSlidingFit, 
    CartesianVector &principalStartPosition, CartesianVector &branchSplitPosition, CartesianVector &branchSplitDirection) const
{
    // Conventions:
    // (1) Delta ray is split from the branch cluster
    // (2) Delta ray occurs where the vertex of the principal cluster meets the vertex of the branch cluster
    //     Method loops over the inner and outer positions of the principal and branch clusters, trying all
    //     possible assignments of vertex and end position until a split is found
    for (unsigned int principalForward = 0; principalForward < 2; ++principalForward)
    {
        const CartesianVector principalVertexPosition(1==principalForward ? principalSlidingFit.GetGlobalMinLayerPosition() : principalSlidingFit.GetGlobalMaxLayerPosition());
        const CartesianVector principalEndPosition(1!=principalForward ? principalSlidingFit.GetGlobalMinLayerPosition() : principalSlidingFit.GetGlobalMaxLayerPosition());
        const CartesianVector principalVertexDirection(1==principalForward ? principalSlidingFit.GetGlobalMinLayerDirection() : principalSlidingFit.GetGlobalMaxLayerDirection() * -1.f);

        CartesianVector projectedBranchPosition(0.f,0.f,0.f);
        bool projectedPositionFound(false), projectedPositionFail(false);

        for (unsigned int branchForward = 0; branchForward < 2; ++branchForward)
        {
            const CartesianVector branchVertexPosition(1==branchForward ? branchSlidingFit.GetGlobalMinLayerPosition() : branchSlidingFit.GetGlobalMaxLayerPosition());
            const CartesianVector branchEndPosition(1!=branchForward ? branchSlidingFit.GetGlobalMinLayerPosition() : branchSlidingFit.GetGlobalMaxLayerPosition());
            const CartesianVector branchEndDirection(1!=branchForward ? branchSlidingFit.GetGlobalMinLayerDirection() : branchSlidingFit.GetGlobalMaxLayerDirection() * -1.f);

            if (principalVertexDirection.GetDotProduct(branchEndDirection) < 0.5f)
                continue;

            if ((principalEndPosition - branchEndPosition).GetMagnitudeSquared() < (principalVertexPosition - branchVertexPosition).GetMagnitudeSquared())
                continue;

            // Project the principal vertex onto the branch cluster
            try
            {
                if (!projectedPositionFound && !projectedPositionFail)
                {
                    projectedBranchPosition = LArPointingClusterHelper::GetProjectedPosition(principalVertexPosition, principalVertexDirection, branchSlidingFit.GetCluster(), m_projectionAngularAllowance);
                    projectedPositionFound = true;
                }
            }
            catch (StatusCodeException &)
            {
                projectedPositionFail = true;
            }

            if (!projectedPositionFound || projectedPositionFail)
                continue;

            const float projectedDistanceSquared((projectedBranchPosition - principalVertexPosition).GetMagnitudeSquared());

            if (projectedDistanceSquared > m_maxLongitudinalDisplacement * m_maxLongitudinalDisplacement)
                continue;

            const float commonDistanceSquared((projectedBranchPosition - branchEndPosition).GetMagnitudeSquared());

            if (projectedDistanceSquared > commonDistanceSquared)
                continue;

            const float replacementDistanceSquared((projectedBranchPosition - principalEndPosition).GetMagnitudeSquared());

            if (replacementDistanceSquared < m_minLongitudinalExtension * m_minLongitudinalExtension)
                continue;

            const float branchDistanceSquared((projectedBranchPosition - branchVertexPosition).GetMagnitudeSquared());

            if (branchDistanceSquared > 4.f * replacementDistanceSquared)
                continue;

            // Require that principal vertex and branch projection have good (and improved) pointing
            bool foundSplit(false);

            const float halfWindowLength(branchSlidingFit.GetLayerFitHalfWindowLength());
            const float deltaL(1==branchForward ? +halfWindowLength : -halfWindowLength);

            float localL(0.f), localT(0.f);
            CartesianVector forwardDirection(0.f,0.f,0.f);
            branchSlidingFit.GetLocalPosition(projectedBranchPosition, localL, localT);

            if (STATUS_CODE_SUCCESS != branchSlidingFit.GetGlobalFitDirection(localL + deltaL, forwardDirection))
                continue;

            CartesianVector projectedBranchDirection(1==branchForward ? forwardDirection : forwardDirection * -1.f);
            const float cosTheta(-projectedBranchDirection.GetDotProduct(principalVertexDirection));

            try
            {
                const float currentCosTheta(branchSlidingFit.GetCosScatteringAngle(localL)); 

                if (cosTheta < currentCosTheta)
                    continue;
            }
            catch (StatusCodeException &)
            {
            }

            float rT1(0.f), rL1(0.f), rT2(0.f), rL2(0.f);
            LArPointingClusterHelper::GetImpactParameters(projectedBranchPosition, projectedBranchDirection, principalVertexPosition, rL1, rT1);
            LArPointingClusterHelper::GetImpactParameters(principalVertexPosition, principalVertexDirection, projectedBranchPosition, rL2, rT2);

            if ((cosTheta > m_minCosRelativeAngle) && (rT1 < m_maxTransverseDisplacement) && (rT2 < m_maxTransverseDisplacement))
            {
                foundSplit = true;
                principalStartPosition = principalVertexPosition;
                branchSplitPosition = projectedBranchPosition;
                branchSplitDirection = projectedBranchDirection * -1.f;
            }

            if (foundSplit)
                return;
        }
    }

    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode BranchSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTransverseDisplacement", m_maxTransverseDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxLongitudinalDisplacement", m_maxLongitudinalDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinLongitudinalExtension", m_minLongitudinalExtension));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCosRelativeAngle", m_minCosRelativeAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ProjectionAngularAllowance", m_projectionAngularAllowance));

    return TwoDSlidingFitSplittingAndSplicingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
