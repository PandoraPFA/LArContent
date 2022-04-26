/**
 *  @file   larpandoracontent/LArHelpers/LArVertexHelper.cc
 *
 *  @brief  Implementation of the vertex helper class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArHelpers/LArVertexHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"

#include <algorithm>
#include <limits>

using namespace pandora;

namespace lar_content
{

LArVertexHelper::ClusterDirection LArVertexHelper::GetClusterDirectionInZ(
    const Pandora &pandora, const Vertex *const pVertex, const Cluster *const pCluster, const float tanAngle, const float apexShift)
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

//-----------------------------------------------------------------------------------------------------------------------------------------

bool LArVertexHelper::IsInFiducialVolume(const CartesianVector &vertex, const std::string &detector)
{
    if (detector == "dune_fd_hd")
    {
        const float x{vertex.GetX()};
        const float y{vertex.GetY()};
        const float z{vertex.GetZ()};
        return -310 < x && x < 310 && -550 < y && y < 550 && 50 < z && z < 1244;
    }
    else
        return false;
}

} // namespace lar_content
