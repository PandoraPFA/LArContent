/**
 *  @file   larpandoracontent/LArCheating/CheatingEventSlicingTool.h
 *
 *  @brief  Header file for the cheating event slicing tool class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_EVENT_SLICING_TOOL_H
#define LAR_CHEATING_EVENT_SLICING_TOOL_H 1

#include "larpandoracontent/LArUtility/ParentSlicingBaseAlgorithm.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  CheatingEventSlicingTool class
 */
class CheatingEventSlicingTool : public SlicingTool
{
public:
    void Slice(const ParentSlicingBaseAlgorithm *const pAlgorithm, const ParentSlicingBaseAlgorithm::HitTypeToNameMap &caloHitListNames,
        const ParentSlicingBaseAlgorithm::HitTypeToNameMap &clusterListNames, ParentSlicingBaseAlgorithm::SliceList &sliceList);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::unordered_map<const pandora::MCParticle*, ParentSlicingBaseAlgorithm::Slice> MCParticleToSliceMap;

    /**
     *  @brief  Fill slices using hits from a specified view
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  hitType the hit type (i.e. view)
     *  @param  caloHitListNames the hit type to calo hit list name map
     *  @param  mcParticleToSliceMap the parent mc particle to slice map
     */
    void FillSlices(const ParentSlicingBaseAlgorithm *const pAlgorithm, const pandora::HitType hitType,
        const ParentSlicingBaseAlgorithm::HitTypeToNameMap &caloHitListNames, MCParticleToSliceMap &mcParticleToSliceMap) const;

    std::string     m_mcParticleListName;           ///< The name of the three d mc particle list name
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_EVENT_SLICING_TOOL_H
