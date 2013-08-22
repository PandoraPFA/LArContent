/**
 *  @file   LArContent/src/LArThreeDSeed/ThreeDPairedTracksAlgorithm.cc
 * 
 *  @brief  Implementation of the 3D seed finding algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArThreeDSeed/ThreeDPairedTracksAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode ThreeDPairedTracksAlgorithm::Run()
{

    std::cout << "  Hello John, I am a ThreeDPairedTracksAlgorithm... " << std::endl;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDPairedTracksAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameU",
        m_inputClusterListNameU));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameV",
        m_inputClusterListNameV));
  
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameW",
        m_inputClusterListNameW));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
