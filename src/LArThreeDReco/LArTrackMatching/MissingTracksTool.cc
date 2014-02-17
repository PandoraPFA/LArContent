/**
 *  @file   LArContent/src/LArThreeDReco/LArTrackMatching/MissingTracksTool.cc
 * 
 *  @brief  Implementation of the missing tracks tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "LArThreeDReco/LArTrackMatching/MissingTracksTool.h"

using namespace pandora;

namespace lar
{

bool MissingTracksTool::Run(const SlidingFitResultMap &slidingFitResultMap, TrackOverlapTensor &overlapTensor, ProtoParticleVector &protoParticleVector)
{
    std::cout << "MissingTracksTool::Run() " << std::endl;
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MissingTracksTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar
