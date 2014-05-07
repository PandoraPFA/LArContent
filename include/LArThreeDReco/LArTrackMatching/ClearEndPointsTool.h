/**
 *  @file   LArContent/include/LArThreeDReco/LArTrackMatching/ClearEndPointsTool.h
 *
 *  @brief  Header file for the clear tracks tool class.
 *
 *  $Log: $
 */
#ifndef CLEAR_END_POINTS_TOOL_H
#define CLEAR_END_POINTS_TOOL_H 1

#include "LArThreeDReco/LArTrackMatching/ThreeDLongitudinalTracksAlgorithm.h"

namespace lar
{

/**
 *  @brief  ClearEndPointsTool class
 */
class ClearEndPointsTool : public LongitudinalTensorTool
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

private:
    bool Run(ThreeDLongitudinalTracksAlgorithm *pAlgorithm, TensorType &overlapTensor);
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Create three dimensional particles for a given tensor element list
     *
     *  @param  pAlgorithm address of the calling algorithm (ultimately responsible for the particles)
     *  @param  elementList the tensor element list
     *  @param  particlesMade receive boolean indicating whether particles have been made
     */
    void CreateThreeDParticles(ThreeDLongitudinalTracksAlgorithm *pAlgorithm, const TensorType::ElementList &elementList, bool &particlesMade) const;


};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *ClearEndPointsTool::Factory::CreateAlgorithmTool() const
{
    return new ClearEndPointsTool();
}

} // namespace lar

#endif // #ifndef CLEAR_END_POINTS_TOOL_H
