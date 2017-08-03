/**
 *  @file   larpandoracontent/LArCheating/CheatingNeutrinoIdTool.cc
 *
 *  @brief  Implementation of the cheating neutrino id tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArCheating/CheatingNeutrinoIdTool.h"

using namespace pandora;

namespace lar_content
{

CheatingNeutrinoIdTool::CheatingNeutrinoIdTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CheatingNeutrinoIdTool::GetNeutrinoSliceIndex(const ParentAlgorithm::SliceIndexToPropertiesMap &/*sliceIndexToPropertiesMap*/, unsigned int &neutrinoSliceIndex) const
{
    // TODO
    neutrinoSliceIndex = 0;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingNeutrinoIdTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
