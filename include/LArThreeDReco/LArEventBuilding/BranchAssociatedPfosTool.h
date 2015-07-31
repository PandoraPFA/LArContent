/**
 *  @file   LArContent/include/LArThreeDReco/LArEventBuilding/BranchAssociatedPfosTool.h
 * 
 *  @brief  Header file for the branch associated pfos tool class.
 * 
 *  $Log: $
 */
#ifndef LAR_BRANCH_ASSOCIATED_PFOS_TOOL_H
#define LAR_BRANCH_ASSOCIATED_PFOS_TOOL_H 1

#include "LArThreeDReco/LArEventBuilding/PfoHierarchyAlgorithm.h"

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

    void Run(PfoHierarchyAlgorithm *const pAlgorithm, const pandora::Vertex *const pNeutrinoVertex, PfoHierarchyAlgorithm::PfoInfoMap &pfoInfoMap);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *BranchAssociatedPfosTool::Factory::CreateAlgorithmTool() const
{
    return new BranchAssociatedPfosTool();
}

} // namespace lar_content

#endif // #ifndef LAR_BRANCH_ASSOCIATED_PFOS_TOOL_H
