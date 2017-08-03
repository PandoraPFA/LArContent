/**
 *  @file   larpandoracontent/LArThreeDReco/LArEventBuilding/NeutrinoIdTool.cc
 *
 *  @brief  Implementation of the neutrino id tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArThreeDReco/LArEventBuilding/NeutrinoIdTool.h"

using namespace pandora;

namespace lar_content
{

NeutrinoIdTool::NeutrinoIdTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoIdTool::FillNeutrinoProperties(const PfoList *const /*pPfoList*/, ParentAlgorithm::SliceProperties &/*sliceProperties*/) const
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoIdTool::FillCosmicRayProperties(const PfoList *const /*pPfoList*/, ParentAlgorithm::SliceProperties &/*sliceProperties*/) const
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool NeutrinoIdTool::GetNeutrinoSliceIndex(const ParentAlgorithm::SliceIndexToPropertiesMap &/*sliceIndexToPropertiesMap*/, unsigned int &neutrinoSliceIndex) const
{
    // TODO
    neutrinoSliceIndex = 0;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode NeutrinoIdTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
