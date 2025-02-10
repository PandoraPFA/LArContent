/**
 *  @file   larpandoracontent/LArHelpers/LArVertexHelper.cc
 *
 *  @brief  Implementation of the vertex helper class.
 *
 *  $Log: $
 */

#include "Geometry/LArTPC.h"
#include "Managers/GeometryManager.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArVertexHelper.h"

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

bool LArVertexHelper::IsInFiducialVolume(const Pandora &pandora, const CartesianVector &vertex, const std::string &detector)
{
    const LArTPCMap &larTPCMap(pandora.GetGeometry()->GetLArTPCMap());

    if (larTPCMap.empty())
    {
        std::cout << "LArVertexHelper::IsInFiducialVolume - LArTPC description not registered with Pandora as required " << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);
    }

    float tpcMinX{std::numeric_limits<float>::max()}, tpcMaxX{-std::numeric_limits<float>::max()};
    float tpcMinY{std::numeric_limits<float>::max()}, tpcMaxY{-std::numeric_limits<float>::max()};
    float tpcMinZ{std::numeric_limits<float>::max()}, tpcMaxZ{-std::numeric_limits<float>::max()};

    for (const auto &[volumeId, pLArTPC] : larTPCMap)
    {
        (void)volumeId;
        const float centreX{pLArTPC->GetCenterX()}, halfWidthX{0.5f * pLArTPC->GetWidthX()};
        const float centreY{pLArTPC->GetCenterY()}, halfWidthY{0.5f * pLArTPC->GetWidthY()};
        const float centreZ{pLArTPC->GetCenterZ()}, halfWidthZ{0.5f * pLArTPC->GetWidthZ()};
        tpcMinX = std::min(tpcMinX, centreX - halfWidthX);
        tpcMaxX = std::max(tpcMaxX, centreX + halfWidthX);
        tpcMinY = std::min(tpcMinY, centreY - halfWidthY);
        tpcMaxY = std::max(tpcMaxY, centreY + halfWidthY);
        tpcMinZ = std::min(tpcMinZ, centreZ - halfWidthZ);
        tpcMaxZ = std::max(tpcMaxZ, centreZ + halfWidthZ);
    }

    if (detector == "dune_fd_hd")
    {
        const float x{vertex.GetX()};
        const float y{vertex.GetY()};
        const float z{vertex.GetZ()};
        return (tpcMinX + 50.f) < x && x < (tpcMaxX - 50.f) && (tpcMinY + 50.f) < y && y < (tpcMaxY - 50.f) && (tpcMinZ + 50.f) < z &&
            z < (tpcMaxZ - 150.f);
    }
    else if (detector == "dune_nd")
    {
        const float x{vertex.GetX()};
        const float y{vertex.GetY()};
        const float z{vertex.GetZ()};
        return tpcMinX < x && x < tpcMaxX && tpcMinY < y && y < tpcMaxY && tpcMinZ < z && z < tpcMaxZ;
    }
    else
    {
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void LArVertexHelper::GetTrueVertexPosition(const CartesianVector &trueVertex, float &x, float &y, float &z)
{
    x = trueVertex.GetX();
    y = trueVertex.GetY();
    z = trueVertex.GetZ();
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void LArVertexHelper::GetTrueVertexPosition(const CartesianVector &trueVertex, const LArTransformationPlugin *const pTransform, float &x, float &u, float &v, float &w)
{
    x = trueVertex.GetX();
    u = static_cast<float>(pTransform->YZtoU(trueVertex.GetY(), trueVertex.GetZ()));
    v = static_cast<float>(pTransform->YZtoV(trueVertex.GetY(), trueVertex.GetZ()));
    w = static_cast<float>(pTransform->YZtoW(trueVertex.GetY(), trueVertex.GetZ()));
}

} // namespace lar_content
