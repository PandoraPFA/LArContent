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

    /**
     *  @brief  Default constructor
     */
    EndAssociatedPfosTool();

    void Run(PfoHierarchyAlgorithm *const pAlgorithm, const pandora::Vertex *const pNeutrinoVertex, PfoHierarchyAlgorithm::PfoInfoMap &pfoInfoMap);

private:
    /**
     *  @brief  Query the pfo info map and separate/extract pfos currently either acting as parents or associated with the neutrino vertex
     * 
     *  @param  pfoInfoMap the pfo info map
     *  @param  assignedPfos to receive the list of assigned pfos
     *  @param  unassignedPfos to receive the list of unassigned pfos
     */
    void SeparatePfos(const PfoHierarchyAlgorithm::PfoInfoMap &pfoInfoMap, pandora::PfoList &assignedPfos, pandora::PfoList &unassignedPfos) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float       m_minVertexLongitudinalDistance;        ///< Vertex association check: min longitudinal distance cut
    float       m_maxVertexLongitudinalDistance;        ///< Vertex association check: max longitudinal distance cut
    float       m_maxVertexTransverseDistance;          ///< Vertex association check: max transverse distance cut
    float       m_vertexAngularAllowance;               ///< Vertex association check: pointing angular allowance in degrees
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *EndAssociatedPfosTool::Factory::CreateAlgorithmTool() const
{
    return new EndAssociatedPfosTool();
}

} // namespace lar_content

#endif // #ifndef LAR_END_ASSOCIATED_PFOS_TOOL_H
