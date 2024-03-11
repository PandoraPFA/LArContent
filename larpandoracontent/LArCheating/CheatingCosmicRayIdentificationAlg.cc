/**
 *  @file   larpandoracontent/LArCheating/CheatingCosmicRayIdentificationAlg.cc
 *
 *  @brief  Implementation of the cheater for the cosmic ray identification algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArCheating/CheatingCosmicRayIdentificationAlg.h"
#include "larpandoracontent/LArCheating/CheatingSliceIdBaseTool.h"

using namespace pandora;

namespace lar_content
{

CheatingCosmicRayIdentificationAlg::CheatingCosmicRayIdentificationAlg() :
    m_maxNeutrinoFraction(0.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingCosmicRayIdentificationAlg::Run()
{
    const PfoList *pPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_inputPfoListName, pPfoList));

    if (!pPfoList)
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "CheatingCosmicRayIdentificationAlg: pfo list " << m_inputPfoListName << " unavailable." << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    PfoList outputPfoList, outputDaughterPfoList;

    for (const ParticleFlowObject *const pPfo : *pPfoList)
    {
        if (!pPfo->GetParentPfoList().empty())
            continue;

        PfoList downstreamPfos;
        LArPfoHelper::GetAllDownstreamPfos(pPfo, downstreamPfos);

        float thisNeutrinoWeight(0.f), thisTotalWeight(0.f);
        CheatingSliceIdBaseTool::GetTargetParticleWeight(&downstreamPfos, thisNeutrinoWeight, thisTotalWeight, LArMCParticleHelper::IsNeutrino);

        if ((thisTotalWeight < std::numeric_limits<float>::epsilon()) || ((thisNeutrinoWeight / thisTotalWeight) < m_maxNeutrinoFraction))
            outputPfoList.push_back(pPfo);
    }

    if (!outputPfoList.empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_inputPfoListName, m_outputPfoListName, outputPfoList));

        if (!outputDaughterPfoList.empty())
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                PandoraContentApi::SaveList(*this, m_inputDaughterPfoListName, m_outputDaughterPfoListName, outputDaughterPfoList));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingCosmicRayIdentificationAlg::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputPfoListName", m_inputPfoListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputPfoListName", m_outputPfoListName));

    m_inputDaughterPfoListName = m_inputPfoListName;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "InputDaughterPfoListName", m_inputDaughterPfoListName));

    m_outputDaughterPfoListName = m_outputPfoListName;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "OutputDaughterPfoListName", m_outputDaughterPfoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxNeutrinoFraction", m_maxNeutrinoFraction));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
