/**
 *  @file   LArContent/src/LArThreeDReco/LArTrackMatching/OvershootTracksTool.cc
 * 
 *  @brief  Implementation of the overshoot tracks tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "LArThreeDReco/LArTrackMatching/OvershootTracksTool.h"

using namespace pandora;

namespace lar
{

StatusCode OvershootTracksTool::Run(const SlidingFitResultMap &slidingFitResultMap, TensorType &overlapTensor, ProtoParticleVector &protoParticleVector)
{
    std::cout << "OvershootTracksTool::Run() " << std::endl;
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode OvershootTracksTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar
