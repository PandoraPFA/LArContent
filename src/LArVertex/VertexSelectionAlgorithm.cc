/**
 *  @file   LArContent/src/LArVertex/VertexSelectionAlgorithm.cc
 * 
 *  @brief  Implementation of the vertex selection algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArVertex/VertexSelectionAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode VertexSelectionAlgorithm::Run()
{
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexSelectionAlgorithm::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
