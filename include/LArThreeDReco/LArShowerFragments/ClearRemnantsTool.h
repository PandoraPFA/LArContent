/**
 *  @file   LArContent/include/LArThreeDReco/LArShowerFragments/ClearRemnantsTool.h
 *
 *  @brief  Header file for the clear remnants tool class.
 *
 *  $Log: $
 */
#ifndef CLEAR_REMNANTS_TOOL_H
#define CLEAR_REMNANTS_TOOL_H 1

#include "LArThreeDReco/LArShowerFragments/ThreeDRemnantsAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  ClearRemnantsTool class
 */
class ClearRemnantsTool : public RemnantTensorTool
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

    bool Run(ThreeDRemnantsAlgorithm *pAlgorithm, TensorType &overlapTensor);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Create three dimensional particles for a given tensor element list
     *
     *  @param  pAlgorithm address of the calling algorithm (ultimately responsible for the particles)
     *  @param  elementList the tensor element list
     *  @param  particlesMade receive boolean indicating whether particles have been made
     */
    void CreateThreeDParticles(ThreeDRemnantsAlgorithm *pAlgorithm, const TensorType::ElementList &elementList, bool &particlesMade) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *ClearRemnantsTool::Factory::CreateAlgorithmTool() const
{
    return new ClearRemnantsTool();
}

} // namespace lar_content

#endif // #ifndef CLEAR_REMNANTS_TOOL_H
