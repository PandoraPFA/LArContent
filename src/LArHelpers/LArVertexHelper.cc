/**
 *  @file   LArContent/src/LArHelpers/LArVertexHelper.cc
 *
 *  @brief  Implementation of the vertex helper class.
 *
 *  $Log: $
 */

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"
#include "LArHelpers/LArPointingClusterHelper.h"
#include "LArHelpers/LArVertexHelper.h"

#include <limits>

using namespace pandora;

namespace lar_content
{

LArVertexHelper::ClusterDirection LArVertexHelper::GetClusterDirectionInZ(const Pandora &pandora, const Vertex *const pVertex,
    const Cluster *const pCluster, const float tanAngle, const float apexShift)
{
    if ((VERTEX_3D != pVertex->GetVertexType()) || (tanAngle < std::numeric_limits<float>::epsilon()))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
    const CartesianVector theVertex2D(LArGeometryHelper::ProjectPosition(pandora, pVertex->GetPosition(), hitType));

    try
    {
        const LArPointingCluster pointingCluster(pCluster);
        const float length((pointingCluster.GetInnerVertex().GetPosition() - pointingCluster.GetOuterVertex().GetPosition()).GetMagnitude());
        const bool innerIsAtLowerZ(pointingCluster.GetInnerVertex().GetPosition().GetZ() < pointingCluster.GetOuterVertex().GetPosition().GetZ());

        float rLInner(std::numeric_limits<float>::max()), rTInner(std::numeric_limits<float>::max());
        float rLOuter(std::numeric_limits<float>::max()), rTOuter(std::numeric_limits<float>::max());
        LArPointingClusterHelper::GetImpactParameters(pointingCluster.GetInnerVertex(), theVertex2D, rLInner, rTInner);
        LArPointingClusterHelper::GetImpactParameters(pointingCluster.GetOuterVertex(), theVertex2D, rLOuter, rTOuter);

        const bool innerIsVertexAssociated(rLInner > (rTInner / tanAngle) - (length * apexShift));
        const bool outerIsVertexAssociated(rLOuter > (rTInner / tanAngle) - (length * apexShift));

        if (innerIsVertexAssociated == outerIsVertexAssociated)
            return DIRECTION_UNKNOWN;

        if ((innerIsVertexAssociated && innerIsAtLowerZ) || (outerIsVertexAssociated && !innerIsAtLowerZ))
            return DIRECTION_FORWARD_IN_Z;

        if ((innerIsVertexAssociated && !innerIsAtLowerZ) || (outerIsVertexAssociated && innerIsAtLowerZ))
            return DIRECTION_BACKWARD_IN_Z;
    }
    catch (StatusCodeException &)
    {
        return DIRECTION_UNKNOWN;
    }

    throw StatusCodeException(STATUS_CODE_FAILURE);
}

} // namespace lar_content
