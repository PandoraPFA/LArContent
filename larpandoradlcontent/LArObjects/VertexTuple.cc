/**
 *  @file   larpandoradlcontent/LArObjects/VerteTuple.cc
 *
 *  @brief  Implementation of the vertex tuple class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoradlcontent/LArObjects/VertexTuple.h"

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

VertexTuple::VertexTuple(const Pandora &pandora, const CartesianVector &vertexU, const CartesianVector &vertexV, const CartesianVector &vertexW) :
    m_pos{0.f, 0.f, 0.f},
    m_chi2{0.f}
{
    LArGeometryHelper::MergeThreePositions3D(pandora, TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, vertexU, vertexV, vertexW, m_pos, m_chi2);
    if (m_chi2 <= 1.f)
    {
        m_components.emplace_back(vertexU);
        m_components.emplace_back(vertexV);
        m_components.emplace_back(vertexW);
    }
    else
    {
        CartesianVector vertexUV(0.f, 0.f, 0.f);
        float chi2UV{0.f};
        LArGeometryHelper::MergeTwoPositions3D(pandora, TPC_VIEW_U, TPC_VIEW_V, vertexU, vertexV, vertexUV, chi2UV);

        CartesianVector vertexUW(0.f, 0.f, 0.f);
        float chi2UW{0.f};
        LArGeometryHelper::MergeTwoPositions3D(pandora, TPC_VIEW_U, TPC_VIEW_W, vertexU, vertexW, vertexUW, chi2UW);

        CartesianVector vertexVW(0.f, 0.f, 0.f);
        float chi2VW{0.f};
        LArGeometryHelper::MergeTwoPositions3D(pandora, TPC_VIEW_V, TPC_VIEW_W, vertexV, vertexW, vertexVW, chi2VW);

        if (chi2UV < m_chi2)
        {
            m_pos = vertexUV;
            m_chi2 = chi2UV;
            m_components.emplace_back(vertexU);
            m_components.emplace_back(vertexV);
        }
        if (chi2UW < m_chi2)
        {
            m_pos = vertexUW;
            m_chi2 = chi2UW;
            m_components.clear();
            m_components.emplace_back(vertexU);
            m_components.emplace_back(vertexW);
        }
        if (chi2VW < m_chi2)
        {
            m_pos = vertexVW;
            m_chi2 = chi2VW;
            m_components.clear();
            m_components.emplace_back(vertexV);
            m_components.emplace_back(vertexW);
        }
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

VertexTuple::VertexTuple(const Pandora &pandora, const CartesianVector &vertex1, const CartesianVector &vertex2, const HitType view1, const HitType view2) :
    m_pos{0.f, 0.f, 0.f},
    m_chi2{0.f}
{
    m_components.emplace_back(vertex1);
    m_components.emplace_back(vertex2);

    LArGeometryHelper::MergeTwoPositions3D(pandora, view1, view2, vertex1, vertex2, m_pos, m_chi2);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

const CartesianVector &VertexTuple::GetPosition() const
{
    return m_pos;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

const CartesianPointVector &VertexTuple::GetComponents() const
{
    return m_components;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

float VertexTuple::GetChi2() const
{
    return m_chi2;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::string VertexTuple::ToString() const
{
    const float x{m_pos.GetX()}, y{m_pos.GetY()}, z{m_pos.GetZ()};
    return "3D pos: (" + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) + ")   X2 = " + std::to_string(m_chi2);
}

} // namespace lar_dl_content
