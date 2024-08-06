/**
 *  @file   larpandoracontent/LArControlFlow/SimpleNeutrinoIdTool.cc
 *
 *  @brief  Implementation of the neutrino id tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArControlFlow/SimpleNeutrinoIdTool.h"

using namespace pandora;

namespace lar_content
{

SimpleNeutrinoIdTool::SimpleNeutrinoIdTool() :
    m_selectAllNeutrinos(true),
    m_selectOnlyFirstSliceNeutrinos(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SimpleNeutrinoIdTool::SelectOutputPfos(const Algorithm *const pAlgorithm, const SliceHypotheses &nuSliceHypotheses,
    const SliceHypotheses &crSliceHypotheses, PfoList &selectedPfos)
{
    if (nuSliceHypotheses.size() != crSliceHypotheses.size())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    for (unsigned int sliceIndex = 0, nSlices = nuSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
    {
        const PfoList &sliceOutput((m_selectAllNeutrinos || (m_selectOnlyFirstSliceNeutrinos && (0 == sliceIndex)))
                ? nuSliceHypotheses.at(sliceIndex)
                : crSliceHypotheses.at(sliceIndex));

        const float score(m_selectAllNeutrinos || (m_selectOnlyFirstSliceNeutrinos && (0 == sliceIndex)) ? 1.f : -1.f);

        for (const ParticleFlowObject *const pPfo : crSliceHypotheses.at(sliceIndex))
        {
            object_creation::ParticleFlowObject::Metadata metadata;
            metadata.m_propertiesToAdd["NuScore"] = score;
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*pAlgorithm, pPfo, metadata));
        }

        for (const ParticleFlowObject *const pPfo : nuSliceHypotheses.at(sliceIndex))
        {
            object_creation::ParticleFlowObject::Metadata metadata;
            metadata.m_propertiesToAdd["NuScore"] = score;
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*pAlgorithm, pPfo, metadata));
        }

        selectedPfos.insert(selectedPfos.end(), sliceOutput.begin(), sliceOutput.end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SimpleNeutrinoIdTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SelectAllNeutrinos", m_selectAllNeutrinos));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "SelectOnlyFirstSliceNeutrinos", m_selectOnlyFirstSliceNeutrinos));

    if (m_selectAllNeutrinos == m_selectOnlyFirstSliceNeutrinos)
    {
        std::cout << "SimpleNeutrinoIdTool::ReadSettings - exactly one of SelectAllNeutrinos and SelectOnlyFirstSliceNeutrinos must be true"
                  << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
