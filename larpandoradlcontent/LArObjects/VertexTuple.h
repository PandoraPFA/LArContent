/**
 *  @file   larpandoradlcontent/LArObjects/VertexTuple.h
 *
 *  @brief  Header file for the vertex tuple object.
 *
 *  $Log: $
 */
#ifndef LAR_DL_VERTEX_TUPLE_H
#define LAR_DL_VERTEX_TUPLE_H 1

using namespace lar_content;

namespace lar_dl_content
{
/**
 *  @brief  VertexTuple class
 */
class VertexTuple
{
public:
    VertexTuple(const pandora::Pandora &pandora, const pandora::CartesianVector &vertexU, const pandora::CartesianVector &vertexV,
        const pandora::CartesianVector &vertexW);

    VertexTuple(const pandora::Pandora &pandora, const pandora::CartesianVector &vertex1, const pandora::CartesianVector &vertex2,
        const pandora::HitType view1, const pandora::HitType view2);

    const pandora::CartesianVector &GetPosition() const;
    const pandora::CartesianPointVector &GetComponents() const;
    float GetChi2() const;
    std::string ToString() const;

private:
    pandora::CartesianVector m_pos;             ///< Calculated 3D position
    float m_chi2;                               ///< Chi squared of calculated position
    pandora::CartesianPointVector m_components; ///< The 2D vertices that contributed to the 3D vertex
};

} // namespace lar_dl_content

#endif // LAR_DL_VERTEX_TUPLE_H
