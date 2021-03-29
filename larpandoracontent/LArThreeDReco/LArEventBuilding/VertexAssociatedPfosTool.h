/**
 *  @file   larpandoracontent/LArThreeDReco/LArEventBuilding/VertexAssociatedPfosTool.h
 *
 *  @brief  Header file for the vertex associated pfos tool class.
 *
 *  $Log: $
 */
#ifndef LAR_VERTEX_ASSOCIATED_PFOS_TOOL_H
#define LAR_VERTEX_ASSOCIATED_PFOS_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArEventBuilding/NeutrinoHierarchyAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  VertexAssociatedPfosTool class
 */
class VertexAssociatedPfosTool : public PfoRelationTool
{
public:
    /**
     *  @brief  Default constructor
     */
    VertexAssociatedPfosTool();

    void Run(const NeutrinoHierarchyAlgorithm *const pAlgorithm, const pandora::Vertex *const pNeutrinoVertex,
        NeutrinoHierarchyAlgorithm::PfoInfoMap &pfoInfoMap);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_minVertexLongitudinalDistance; ///< Vertex association check: min longitudinal distance cut
    float m_maxVertexLongitudinalDistance; ///< Vertex association check: max longitudinal distance cut
    float m_maxVertexTransverseDistance;   ///< Vertex association check: max transverse distance cut
    float m_vertexAngularAllowance;        ///< Vertex association check: pointing angular allowance in degrees
};

} // namespace lar_content

#endif // #ifndef LAR_VERTEX_ASSOCIATED_PFOS_TOOL_H
