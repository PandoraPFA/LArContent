/**
 *  @file   LArContent/src/LArThreeDReco/LArTrackMatching/MissingTrackTool.cc
 * 
 *  @brief  Implementation of the missing track tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "LArThreeDReco/LArTrackMatching/MissingTrackTool.h"

using namespace pandora;

namespace lar
{

bool MissingTrackTool::Run(ThreeDTransverseTracksAlgorithm *pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraSettings::ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << m_algorithmToolType << std::endl;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MissingTrackTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar
