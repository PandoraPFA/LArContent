/**
 *  @file   LArContent/src/LArUtility/KDTreeLinkerToolsT.cc
 * 
 *  @brief  Implementation of the kd tree linker tools template class
 * 
 *  $Log: $
 */

#include "larpandoracontent/LArContent/LArUtility/KDTreeLinkerToolsT.h"

namespace lar_content
{

std::pair<float, float> minmax(const float a, const float b)
{
    return ((b < a) ? std::pair<float, float>(b, a) : std::pair<float, float>(a, b));
}

//------------------------------------------------------------------------------------------------------------------------------------------

KDTreeTesseract fill_and_bound_4d_kd_tree(const pandora::Algorithm *const caller, const pandora::CaloHitList &points,
    std::vector<KDTreeNodeInfoT<const pandora::CaloHit*, 4> > &nodes, bool passthru)
{
    std::array<float, 4> minpos{ {0.f, 0.f, 0.f, 0.f} }, maxpos{ {0.f, 0.f, 0.f, 0.f} };

    unsigned i = 0;

    for (const pandora::CaloHit *const point : points)
    {
        if (!passthru && !PandoraContentApi::IsAvailable(*caller, point))
            continue;

        const pandora::CartesianVector &pos = kdtree_type_adaptor<const pandora::CaloHit>::position(point);
        const float layer = static_cast<float>(point->GetPseudoLayer());
        nodes.emplace_back(point, pos.GetX(), pos.GetY(), pos.GetZ(), layer);

        if (0 == i)
        {
            minpos[0] = pos.GetX(); minpos[1] = pos.GetY(); minpos[2] = pos.GetZ(); minpos[3] = layer;
            maxpos[0] = pos.GetX(); maxpos[1] = pos.GetY(); maxpos[2] = pos.GetZ(); maxpos[3] = layer;
        }
        else
        {
            minpos[0] = std::min(pos.GetX(), minpos[0]);
            minpos[1] = std::min(pos.GetY(), minpos[1]);
            minpos[2] = std::min(pos.GetZ(), minpos[2]);
            minpos[3] = std::min(layer, minpos[3]);
            maxpos[0] = std::max(pos.GetX(), maxpos[0]);
            maxpos[1] = std::max(pos.GetY(), maxpos[1]);
            maxpos[2] = std::max(pos.GetZ(), maxpos[2]);
            maxpos[3] = std::max(layer, maxpos[3]);
        }

        ++i;
    }

    return KDTreeTesseract(minpos[0], maxpos[0], minpos[1], maxpos[1], minpos[2], maxpos[2], minpos[3], maxpos[3]);
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

//------------------------------------------------------------------------------------------------------------------------------------------

KDTreeTesseract build_4d_kd_search_region(const pandora::CaloHit *const point, const float x_span, const float y_span, const float z_span,
    const float search_layer)
{
    return build_4d_kd_search_region(point->GetPositionVector(), x_span, y_span, z_span, search_layer);
}

//------------------------------------------------------------------------------------------------------------------------------------------

KDTreeTesseract build_4d_kd_search_region(const pandora::CartesianVector &pos, const float x_span, const float y_span, const float z_span,
    const float search_layer)
{
    const auto x_side = minmax(pos.GetX() + x_span, pos.GetX() - x_span);
    const auto y_side = minmax(pos.GetY() + y_span, pos.GetY() - y_span);
    const auto z_side = minmax(pos.GetZ() + z_span, pos.GetZ() - z_span);
    const auto layer_side = minmax(search_layer + 0.5f, search_layer - 0.5f);

    return KDTreeTesseract(x_side.first, x_side.second, y_side.first, y_side.second, z_side.first, z_side.second, layer_side.first, layer_side.second);
}

} // namespace lar_content
