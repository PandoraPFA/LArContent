/**
 *  @file   LArContent/include/LArThreeDReco/LArEventBuilding/EndAssociatedPfosTool.h
 * 
 *  @brief  Header file for the end associated pfos tool class.
 * 
 *  $Log: $
 */
#ifndef LAR_END_ASSOCIATED_PFOS_TOOL_H
#define LAR_END_ASSOCIATED_PFOS_TOOL_H 1

#include "LArThreeDReco/LArEventBuilding/PfoHierarchyAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  EndAssociatedPfosTool class
 */
class EndAssociatedPfosTool : public PfoRelationTool
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

inline pandora::AlgorithmTool *EndAssociatedPfosTool::Factory::CreateAlgorithmTool() const
{
    return new EndAssociatedPfosTool();
}

} // namespace lar_content

#endif // #ifndef LAR_END_ASSOCIATED_PFOS_TOOL_H
