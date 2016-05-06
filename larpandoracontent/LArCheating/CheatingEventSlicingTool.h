/**
 *  @file   LArContent/LArCheating/CheatingEventSlicingTool.h
 * 
 *  @brief  Header file for the cheating event slicing tool class.
 * 
 *  $Log: $
 */
#ifndef LAR_CHEATING_EVENT_SLICING_TOOL_H
#define LAR_CHEATING_EVENT_SLICING_TOOL_H 1

#include "larpandoracontent/LArUtility/NeutrinoParentAlgorithm.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  CheatingEventSlicingTool class
 */
class CheatingEventSlicingTool : public SlicingTool
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

    void Slice(const NeutrinoParentAlgorithm *const pAlgorithm, const NeutrinoParentAlgorithm::HitTypeToNameMap &caloHitListNames,
        const NeutrinoParentAlgorithm::HitTypeToNameMap &clusterListNames, NeutrinoParentAlgorithm::SliceList &sliceList);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::unordered_map<const pandora::MCParticle*, NeutrinoParentAlgorithm::Slice> MCParticleToSliceMap;

    /**
     *  @brief  Fill slices using hits from a specified view
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  hitType the hit type (i.e. view)
     *  @param  caloHitListNames the hit type to calo hit list name map
     *  @param  mcParticleToSliceMap the parent mc particle to slice map
     */
    void FillSlices(const NeutrinoParentAlgorithm *const pAlgorithm, const pandora::HitType hitType,
        const NeutrinoParentAlgorithm::HitTypeToNameMap &caloHitListNames, MCParticleToSliceMap &mcParticleToSliceMap) const;

    std::string     m_mcParticleListName;           ///< The name of the three d mc particle list name
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *CheatingEventSlicingTool::Factory::CreateAlgorithmTool() const
{
    return new CheatingEventSlicingTool();
}

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_EVENT_SLICING_TOOL_H
