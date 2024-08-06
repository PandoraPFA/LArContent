/**
 *  @file   larpandoracontent/LArControlFlow/StitchingBaseTool.h
 *
 *  @brief  Header file for the stitching tool base class.
 *
 *  $Log: $
 */
#ifndef LAR_STITCHING_BASE_TOOL_H
#define LAR_STITCHING_BASE_TOOL_H 1

#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"

namespace lar_content
{

class MasterAlgorithm;

typedef std::unordered_map<const pandora::ParticleFlowObject *, const pandora::LArTPC *> PfoToLArTPCMap;
typedef std::unordered_map<const pandora::ParticleFlowObject *, float> PfoToFloatMap;
/**
 *  @brief  StitchingBaseTool class
 */
class StitchingBaseTool : public pandora::AlgorithmTool
{
public:
    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  pMultiPfoList the list of pfos in multiple lar tpcs
     *  @param  pfoToLArTPCMap the pfo to lar tpc map
     *  @param  stitchedPfosToX0Map a map of cosmic-ray pfos that have been stitched between lar tpcs to the X0 shift
     */
    virtual void Run(const MasterAlgorithm *const pAlgorithm, const pandora::PfoList *const pMultiPfoList, PfoToLArTPCMap &pfoToLArTPCMap,
        PfoToFloatMap &stitchedPfosToX0Map) = 0;
};

} // namespace lar_content

#endif // #ifndef LAR_STITCHING_BASE_TOOL_H
