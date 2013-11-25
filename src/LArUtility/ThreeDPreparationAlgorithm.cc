/**
 *  @file   LArContent/src/LArUtility/ThreeDPreparationAlgorithm.cc
 * 
 *  @brief  Implementation of the three dimensional preparation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArVertexHelper.h"

#include "LArUtility/ThreeDPreparationAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode ThreeDPreparationAlgorithm::Run()
{
    LArVertexHelper::SetCurrentVertex(m_vertexName);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDPreparationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "VertexName", m_vertexName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
