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

StatusCode MissingTracksTool::Run(OverlapTensor<TrackOverlapResult> &overlapTensor)
{
    std::cout << "MissingTracksTool::Run() " << std::endl;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MissingTracksTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar
