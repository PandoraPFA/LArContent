/**
 *  @file   larpandoracontent/LArThreeDReco/LArEventBuilding/EndAssociatedPfosTool.h
 *
 *  @brief  Header file for the end associated pfos tool class.
 *
 *  $Log: $
 */
#ifndef LAR_END_ASSOCIATED_PFOS_TOOL_H
#define LAR_END_ASSOCIATED_PFOS_TOOL_H 1

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArThreeDReco/LArEventBuilding/NeutrinoHierarchyAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  EndAssociatedPfosTool class
 */
class EndAssociatedPfosTool : public PfoRelationTool
{
public:
    /**
     *  @brief  Default constructor
     */
    EndAssociatedPfosTool();

    void Run(const NeutrinoHierarchyAlgorithm *const pAlgorithm, const pandora::Vertex *const pNeutrinoVertex,
        NeutrinoHierarchyAlgorithm::PfoInfoMap &pfoInfoMap);

private:
    /**
     *  @brief  Whether a daughter 3D cluster is in close proximity to the endpoint of a parent 3D cluster
     *
     *  @param  parentEndpoint the parent endpoint position
     *  @param  pParentCluster3D the address of the parent 3D cluster
     *  @param  pDaughterCluster3D the address of the daughter 3D cluster
     *
     *  @return boolean
     */
    bool IsCloseToParentEndpoint(const pandora::CartesianVector &parentEndpoint, const pandora::Cluster *const pParentCluster3D,
        const pandora::Cluster *const pDaughterCluster3D) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_minNeutrinoVertexDistance;     ///< Min distance between candidate parent endpoint and neutrino vertex
    float m_minVertexLongitudinalDistance; ///< Vertex association check: min longitudinal distance cut
    float m_maxVertexLongitudinalDistance; ///< Vertex association check: max longitudinal distance cut
    float m_maxVertexTransverseDistance;   ///< Vertex association check: max transverse distance cut
    float m_vertexAngularAllowance;        ///< Vertex association check: pointing angular allowance in degrees
    float m_maxParentEndpointDistance;     ///< Max distance between candidate parent endpoint and candidate daughter
};

} // namespace lar_content

#endif // #ifndef LAR_END_ASSOCIATED_PFOS_TOOL_H
