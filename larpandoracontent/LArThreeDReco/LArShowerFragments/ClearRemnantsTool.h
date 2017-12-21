/**
 *  @file   larpandoracontent/LArThreeDReco/LArShowerFragments/ClearRemnantsTool.h
 *
 *  @brief  Header file for the clear remnants tool class.
 *
 *  $Log: $
 */
#ifndef CLEAR_REMNANTS_TOOL_H
#define CLEAR_REMNANTS_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArShowerFragments/ThreeDRemnantsAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  ClearRemnantsTool class
 */
class ClearRemnantsTool : public RemnantTensorTool
{
public:
    bool Run(ThreeDRemnantsAlgorithm *const pAlgorithm, TensorType &overlapTensor);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Create three dimensional particles for a given tensor element list
     *
     *  @param  pAlgorithm address of the calling algorithm (ultimately responsible for the particles)
     *  @param  elementList the tensor element list
     *  @param  particlesMade receive boolean indicating whether particles have been made
     */
    void CreateThreeDParticles(ThreeDRemnantsAlgorithm *const pAlgorithm, const TensorType::ElementList &elementList, bool &particlesMade) const;
};

} // namespace lar_content

#endif // #ifndef CLEAR_REMNANTS_TOOL_H
