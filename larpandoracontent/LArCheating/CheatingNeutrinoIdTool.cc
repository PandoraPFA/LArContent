/**
 *  @file   larpandoracontent/LArCheating/CheatingNeutrinoIdTool.cc
 *
 *  @brief  Implementation of the cheating neutrino id tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArCheating/CheatingNeutrinoIdTool.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

using namespace pandora;

namespace lar_content
{

void CheatingNeutrinoIdTool::SelectOutputPfos(const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, PfoList &selectedPfos)
{
    if (nuSliceHypotheses.size() != crSliceHypotheses.size())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    float bestNeutrinoWeight(0.f);
    unsigned int bestSliceIndex(std::numeric_limits<unsigned int>::max());

    for (unsigned int sliceIndex = 0, nSlices = nuSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
    {
        float neutrinoWeight(0.f);
        const PfoList &neutrinoPfoList(nuSliceHypotheses.at(sliceIndex));

        for (const Pfo *const pNeutrinoPfo : neutrinoPfoList)
        {
            if (!LArPfoHelper::IsNeutrino(pNeutrinoPfo))
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

            PfoList downstreamPfos;
            LArPfoHelper::GetAllDownstreamPfos(pNeutrinoPfo, downstreamPfos);

            float thisNeutrinoWeight(0.f), thisTotalWeight(0.f);
            LArMCParticleHelper::GetNeutrinoWeight(&downstreamPfos, thisNeutrinoWeight, thisTotalWeight);
            neutrinoWeight += thisNeutrinoWeight;
        }

        if (neutrinoWeight > bestNeutrinoWeight)
        {
            neutrinoWeight = bestNeutrinoWeight;
            bestSliceIndex = sliceIndex;
        }
    }

    for (unsigned int sliceIndex = 0, nSlices = nuSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
    {
        const PfoList &sliceOutput((bestSliceIndex == sliceIndex) ? nuSliceHypotheses.at(sliceIndex) : crSliceHypotheses.at(sliceIndex));
        selectedPfos.insert(selectedPfos.end(), sliceOutput.begin(), sliceOutput.end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingNeutrinoIdTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
