/**
 *  @file   LArContent/src/LArThreeDReco/LArTrackMatching/UndershootTracksTool.cc
 * 
 *  @brief  Implementation of the undershoot tracks tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "LArThreeDReco/LArTrackMatching/UndershootTracksTool.h"

using namespace pandora;

namespace lar
{

StatusCode UndershootTracksTool::Run(OverlapTensor<TrackOverlapResult> &overlapTensor)
{
    std::cout << "UndershootTracksTool::Run() " << std::endl;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode UndershootTracksTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar
