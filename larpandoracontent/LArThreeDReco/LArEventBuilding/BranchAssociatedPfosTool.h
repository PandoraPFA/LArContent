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
     *  @brief  Factory class for instantiating algorithm tool
     */
    class Factory : public pandora::AlgorithmToolFactory
    {
    public:
        pandora::AlgorithmTool *CreateAlgorithmTool() const;
    };

    /**
     *  @brief  Default constructor
     */
    BranchAssociatedPfosTool();

    void Run(NeutrinoHierarchyAlgorithm *const pAlgorithm, const pandora::Vertex *const pNeutrinoVertex, NeutrinoHierarchyAlgorithm::PfoInfoMap &pfoInfoMap);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float       m_minNeutrinoVertexDistance;        ///< Branch association: min distance from branch vertex to neutrino vertex
    float       m_maxParentClusterDistance;         ///< Branch association: max distance from branch vertex to a hit in parent 3D cluster
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *BranchAssociatedPfosTool::Factory::CreateAlgorithmTool() const
{
    return new BranchAssociatedPfosTool();
}

} // namespace lar_content

#endif // #ifndef LAR_BRANCH_ASSOCIATED_PFOS_TOOL_H
