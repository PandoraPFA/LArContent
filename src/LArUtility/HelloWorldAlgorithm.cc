/**
 *  @file   LArContent/src/LArUtility/HelloWorldAlgorithm.cc
 * 
 *  @brief  Implementation of the list merging algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArUtility/HelloWorldAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode HelloWorldAlgorithm::Run()
{
    std::cout << "  Hello World !!  " << std::endl;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HelloWorldAlgorithm::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar
