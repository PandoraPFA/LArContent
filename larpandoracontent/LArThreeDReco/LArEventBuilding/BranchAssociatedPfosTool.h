/**
 *  @file   larpandoracontent/LArThreeDReco/LArEventBuilding/BranchAssociatedPfosTool.h
 *
 *  @brief  Header file for the branch associated pfos tool class.
 *
 *  $Log: $
 */
#ifndef LAR_BRANCH_ASSOCIATED_PFOS_TOOL_H
#define LAR_BRANCH_ASSOCIATED_PFOS_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArEventBuilding/NeutrinoHierarchyAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  BranchAssociatedPfosTool class
 */
class BranchAssociatedPfosTool : public PfoRelationTool
{
public:
    /**
     *  @brief  Default constructor
     */
    BranchAssociatedPfosTool();

    void Run(const NeutrinoHierarchyAlgorithm *const pAlgorithm, const pandora::Vertex *const pNeutrinoVertex,
        NeutrinoHierarchyAlgorithm::PfoInfoMap &pfoInfoMap);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_minNeutrinoVertexDistance;   ///< Branch association: min distance from branch vertex to neutrino vertex
    float m_trackBranchAdditionFraction; ///< Branch association: min fraction of length along parent track before association allowed
    float m_maxParentClusterDistance;    ///< Branch association: max distance from branch vertex to a hit in parent 3D cluster
};

} // namespace lar_content

#endif // #ifndef LAR_BRANCH_ASSOCIATED_PFOS_TOOL_H
