/**
 *  @file   LArContent/include/LArThreeDReco/LArEventBuilding/EventSlicingTool.h
 * 
 *  @brief  Header file for the event slicing tool class.
 * 
 *  $Log: $
 */
#ifndef LAR_EVENT_SLICING_TOOL_H
#define LAR_EVENT_SLICING_TOOL_H 1

#include "LArUtility/NeutrinoParentAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  EventSlicingTool class
 */
class EventSlicingTool : public SlicingTool
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

    void Slice(NeutrinoParentAlgorithm *const pAlgorithm, const NeutrinoParentAlgorithm::HitTypeToNameMap &caloHitListNames,
        const NeutrinoParentAlgorithm::HitTypeToNameMap &clusterListNames, NeutrinoParentAlgorithm::SliceList &sliceList);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *EventSlicingTool::Factory::CreateAlgorithmTool() const
{
    return new EventSlicingTool();
}

} // namespace lar_content

#endif // #ifndef LAR_EVENT_SLICING_TOOL_H
