/**
 *  @file   larpandoracontent/LArCheating/CheatingCosmicRayRemovalAlgorithm.cc
 *
 *  @brief  Implementation of the cheating cosmic ray removal algorithm class.
 *
 *  $Log: $
 */

#include "Helpers/MCParticleHelper.h"
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArCheating/CheatingCosmicRayRemovalAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode CheatingCosmicRayRemovalAlgorithm::Run()
{
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputCaloHitListName, pCaloHitList));

    CaloHitList outputCaloHitList;

    for (const CaloHit *pCaloHit : *pCaloHitList)
    {
        try
        {
            const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
            if (!LArMCParticleHelper::IsCosmicRay(LArMCParticleHelper::GetParentMCParticle(pMCParticle)))
                outputCaloHitList.push_back(pCaloHit);
        }
        catch (const StatusCodeException &)
        {
            std::cout << "CheatingCosmicRayRemovalAlgorithm::Run - Unable to determine MCParticle origin for an input CaloHit, which will be skipped."
                      << std::endl;
            continue;
        }
    }

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, outputCaloHitList, m_outputCaloHitListName));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<CaloHit>(*this, m_outputCaloHitListName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingCosmicRayRemovalAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListName", m_inputCaloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputCaloHitListName", m_outputCaloHitListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
