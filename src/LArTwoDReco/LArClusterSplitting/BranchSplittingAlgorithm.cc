/**
 *  @file   LArContent/src/LArTwoDReco/ClusterSplitting/BranchSplittingAlgorithm.cc
 *
 *  @brief  Implementation of the branch splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArPointingClusterHelper.h"

#include "LArTwoDReco/LArClusterSplitting/BranchSplittingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

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
            CartesianVector projectedBranchPosition(0.f,0.f,0.f);
            float projectedDistanceSquared(std::numeric_limits<float>::max());
            float branchDistanceSquared(0.f), replacementDistanceSquared(0.f), commonDistanceSquared(0.f);

            try
            {
                projectedBranchPosition = LArPointingClusterHelper::GetProjectedPosition(principalVertexPosition, principalVertexDirection,
                    branchSlidingFit.GetCluster(), m_projectionAngularAllowance);
                projectedDistanceSquared   = (projectedBranchPosition - principalVertexPosition).GetMagnitudeSquared();
                replacementDistanceSquared = (projectedBranchPosition - principalEndPosition).GetMagnitudeSquared();
                branchDistanceSquared      = (projectedBranchPosition - branchVertexPosition).GetMagnitudeSquared(); 
                commonDistanceSquared      = (projectedBranchPosition - branchEndPosition).GetMagnitudeSquared(); 
            }
            catch (StatusCodeException &)
            {
            }

            if (projectedDistanceSquared > m_maxLongitudinalDisplacement * m_maxLongitudinalDisplacement)
                continue;

            if (projectedDistanceSquared > commonDistanceSquared)
                continue;

            if (replacementDistanceSquared < m_minLongitudinalExtension * m_minLongitudinalExtension)
                continue;

            if (branchDistanceSquared > 4.f * replacementDistanceSquared)
                continue;

            // Require that principal vertex and branch projection have good (and improved) pointing
            bool foundSplit(false);

            const float halfWindowLength(branchSlidingFit.GetLayerFitHalfWindowLength());
            const float deltaL(1==branchForward ? +halfWindowLength : -halfWindowLength);

            try
            {
                float localL(0.f), localT(0.f);
                CartesianVector forwardDirection(0.f,0.f,0.f);
                branchSlidingFit.GetLocalPosition(projectedBranchPosition, localL, localT);
                branchSlidingFit.GetGlobalFitDirection(localL + deltaL, forwardDirection);

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
            }
            catch (StatusCodeException &)
            {
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
    m_maxTransverseDisplacement = 1.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTransverseDisplacement", m_maxTransverseDisplacement));

    m_maxLongitudinalDisplacement = 10.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxLongitudinalDisplacement", m_maxLongitudinalDisplacement));

    m_minLongitudinalExtension = 3.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinLongitudinalExtension", m_minLongitudinalExtension));

    m_minCosRelativeAngle = 0.966;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCosRelativeAngle", m_minCosRelativeAngle));

    m_projectionAngularAllowance = 20.f; // degrees
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ProjectionAngularAllowance", m_projectionAngularAllowance));

    return TwoDSlidingFitSplittingAndSplicingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
