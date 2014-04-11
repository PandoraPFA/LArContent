/**
 *  @file   LArContent/src/LArCheating/CheatingPfoCreationAlgorithm.cc
 * 
 *  @brief  Implementation of the cheating cluster creation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArCheating/CheatingPfoCreationAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode CheatingPfoCreationAlgorithm::Run()
{

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingPfoCreationAlgorithm::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar
