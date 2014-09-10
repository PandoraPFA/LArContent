/**
 *  @file   LArContent/src/LArThreeDReco/LArCosmicRay/CosmicRayIdentificationAlgorithm.cc
 * 
 *  @brief  Implementation of the cosmic ray identification algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArThreeDReco/LArCosmicRay/CosmicRayIdentificationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode CosmicRayIdentificationAlgorithm::Run()
{
    return STATUS_CODE_NOT_INITIALIZED;  
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayIdentificationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputPfoListName", m_inputPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputPfoListName", m_outputPfoListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
