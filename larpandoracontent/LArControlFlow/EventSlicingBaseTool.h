/**
 *  @file   larpandoracontent/LArControlFlow/EventSlicingBaseTool.h
 *
 *  @brief  Header file for the event slicing tool base class.
 *
 *  $Log: $
 */
#ifndef LAR_EVENT_SLICING_BASE_TOOL_H
#define LAR_EVENT_SLICING_BASE_TOOL_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArObjects/LArSlice.h"

namespace lar_content
{

typedef std::map<pandora::HitType, std::string> HitTypeToNameMap;
/**
 *  @brief  EventSlicingBaseTool class
 */
class EventSlicingBaseTool : public pandora::AlgorithmTool
{
public:
    /**
     *  @brief  Run the slicing tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  caloHitListNames the hit type to calo hit list name map
     *  @param  clusterListNames the hit type to cluster list name map
     *  @param  sliceList to receive the populated slice list
     */
    virtual void RunSlicing(const pandora::Algorithm *const pAlgorithm, const HitTypeToNameMap &caloHitListNames,
        const HitTypeToNameMap &clusterListNames, SliceList &sliceList) = 0;
};

}

#endif // #ifndef LAR_EVENT_SLICING_BASE_TOOL_H
