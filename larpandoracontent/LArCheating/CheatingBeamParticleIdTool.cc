/**
 *  @file   larpandoracontent/LArCheating/CheatingBeamParticleIdTool.cc
 *
 *  @brief  Implementation of the cheating beam particle id tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArCheating/CheatingBeamParticleIdTool.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

using namespace pandora;

namespace lar_content
{

CheatingBeamParticleIdTool::CheatingBeamParticleIdTool() :
    m_minWeightFraction(0.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingBeamParticleIdTool::SelectOutputPfos(const pandora::Algorithm *const /*pAlgorithm*/,
    const SliceHypotheses &testBeamSliceHypotheses, const SliceHypotheses &crSliceHypotheses, PfoList &selectedPfos)
{
    if (testBeamSliceHypotheses.size() != crSliceHypotheses.size())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    for (unsigned int sliceIndex = 0, nSlices = testBeamSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
    {
        float beamParticleWeight(0.f), totalWeight(0.f);
        const PfoList &testBeamPfoList(testBeamSliceHypotheses.at(sliceIndex));

        for (const Pfo *const pTestBeamPfo : testBeamPfoList)
        {
            if (!LArPfoHelper::IsTestBeam(pTestBeamPfo))
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

            PfoList downstreamPfos;
            LArPfoHelper::GetAllDownstreamPfos(pTestBeamPfo, downstreamPfos);

            float thisBeamParticleWeight(0.f), thisTotalWeight(0.f);
            CheatingSliceIdBaseTool::GetTargetParticleWeight(&downstreamPfos, thisBeamParticleWeight, thisTotalWeight, LArMCParticleHelper::IsBeamParticle);

            beamParticleWeight += thisBeamParticleWeight;
            totalWeight += thisTotalWeight;
        }

        const float beamWeightFraction(totalWeight < std::numeric_limits<float>::epsilon() ? 0.f : beamParticleWeight / totalWeight);

        if (beamWeightFraction > m_minWeightFraction)
        {
            const PfoList &sliceOutput(testBeamSliceHypotheses.at(sliceIndex));
            selectedPfos.insert(selectedPfos.end(), sliceOutput.begin(), sliceOutput.end());
        }
        else
        {
            const PfoList &sliceOutput(crSliceHypotheses.at(sliceIndex));
            selectedPfos.insert(selectedPfos.end(), sliceOutput.begin(), sliceOutput.end());
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingBeamParticleIdTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MinimumWeightFraction", m_minWeightFraction));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
