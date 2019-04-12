/**
 *  @file   larpandoracontent/LArUtility/KDTreeLinkerToolsT.cc
 *
 *  @brief  Implementation of the kd tree linker tools template class
 *
 *  $Log: $
 */

#include "larpandoracontent/LArUtility/KDTreeLinkerToolsT.h"

namespace lar_content
{

std::pair<float, float> minmax(const float a, const float b)
{
    return ((b < a) ? std::pair<float, float>(b, a) : std::pair<float, float>(a, b));
}

//------------------------------------------------------------------------------------------------------------------------------------------

KDTreeBox build_2d_kd_search_region(const pandora::CaloHit *const point, const float x_span, const float z_span)
{
    return build_2d_kd_search_region(point->GetPositionVector(), x_span, z_span);
}

//------------------------------------------------------------------------------------------------------------------------------------------

KDTreeBox build_2d_kd_search_region(const pandora::CartesianVector &pos, const float x_span, const float z_span)
{
    const auto x_side = minmax(pos.GetX() + x_span, pos.GetX() - x_span);
    const auto z_side = minmax(pos.GetZ() + z_span, pos.GetZ() - z_span);

    return KDTreeBox(x_side.first, x_side.second, z_side.first, z_side.second);
}

//------------------------------------------------------------------------------------------------------------------------------------------

KDTreeCube build_3d_kd_search_region(const pandora::CaloHit *const point, const float x_span, const float y_span, const float z_span)
{
    return build_3d_kd_search_region(point->GetPositionVector(), x_span, y_span, z_span);
}

//------------------------------------------------------------------------------------------------------------------------------------------

KDTreeCube build_3d_kd_search_region(const pandora::CartesianVector &pos, const float x_span, const float y_span, const float z_span)
{
    const auto x_side = minmax(pos.GetX() + x_span, pos.GetX() - x_span);
    const auto y_side = minmax(pos.GetY() + y_span, pos.GetY() - y_span);
    const auto z_side = minmax(pos.GetZ() + z_span, pos.GetZ() - z_span);

    return KDTreeCube(x_side.first, x_side.second, y_side.first, y_side.second, z_side.first, z_side.second);
}

} // namespace lar_content
