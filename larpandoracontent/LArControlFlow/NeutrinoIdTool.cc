/**
 *  @file   larpandoracontent/LArControlFlow/NeutrinoIdTool.cc
 *
 *  @brief  Implementation of the neutrino id tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArControlFlow/NeutrinoIdTool.h"

using namespace pandora;

namespace lar_content
{

void NeutrinoIdTool::SelectOutputPfos(const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, PfoList &selectedPfos)
{
    if (nuSliceHypotheses.size() != crSliceHypotheses.size())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    for (unsigned int sliceIndex = 0, nSlices = nuSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
    {
        // TODO
        const PfoList &sliceOutput((0 == sliceIndex) ? nuSliceHypotheses.at(sliceIndex) : crSliceHypotheses.at(sliceIndex));
        selectedPfos.insert(selectedPfos.end(), sliceOutput.begin(), sliceOutput.end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode NeutrinoIdTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
