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

bool OvershootTracksTool::Run(const SlidingFitResultMap &slidingFitResultMap, TrackOverlapTensor &overlapTensor, ProtoParticleVector &protoParticleVector)
{
    std::cout << "OvershootTracksTool::Run() " << std::endl;
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode OvershootTracksTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar
