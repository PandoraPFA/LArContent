/**
 *  @file   LArContent/include/LArThreeDReco/LArEventBuilding/VertexAssociatedPfosTool.h
 * 
 *  @brief  Header file for the vertex associated pfos tool class.
 * 
 *  $Log: $
 */
#ifndef LAR_VERTEX_ASSOCIATED_PFOS_TOOL_H
#define LAR_VERTEX_ASSOCIATED_PFOS_TOOL_H 1

#include "LArThreeDReco/LArEventBuilding/PfoHierarchyAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  VertexAssociatedPfosTool class
 */
class VertexAssociatedPfosTool : public PfoRelationTool
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

inline pandora::AlgorithmTool *VertexAssociatedPfosTool::Factory::CreateAlgorithmTool() const
{
    return new VertexAssociatedPfosTool();
}

} // namespace lar_content

#endif // #ifndef LAR_VERTEX_ASSOCIATED_PFOS_TOOL_H
