/**
 *  @file   LArContent/src/LArClusterAssociation/TransverseExtensionAlgorithm.cc
 * 
 *  @brief  Implementation of the transverse extension algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArClusterAssociation/TransverseExtensionAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode TransverseExtensionAlgorithm::Run()
{
    std::cout << "  Hello World !!  " << std::endl;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TransverseExtensionAlgorithm::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar
