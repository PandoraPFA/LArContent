/**
 *  @file   larpandoradlcontent/LArControlFlow/DLMasterAlgorithm.cc
 *
 *  @brief  Implementation of the master algorithm class.
 *
 *  $Log: $
 */

#include "Api/PandoraApi.h"

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoradlcontent/LArControlFlow/DLMasterAlgorithm.h"
#include "larpandoradlcontent/LArDLContent.h"

using namespace pandora;

namespace lar_dl_content
{

StatusCode DLMasterAlgorithm::Run()
{
    return MasterAlgorithm::Run();
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLMasterAlgorithm::RegisterCustomContent(const Pandora *const pPandora) const
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArDLContent::RegisterAlgorithms(*pPandora));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLMasterAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    return MasterAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_dl_content
