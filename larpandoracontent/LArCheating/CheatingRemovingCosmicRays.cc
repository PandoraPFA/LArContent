/**
 *  @file   larpandoracontent/LArCheating/CheatingRemovingCosmicRays.cc
 * 
 *  @brief  Implementation of the cheating removing cosmic rays algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Helpers/MCParticleHelper.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArCheating/CheatingRemovingCosmicRays.h"

using namespace pandora;

namespace lar_content
{

StatusCode CheatingRemovingCosmicRays::Run()
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
            if (!LArMCParticleHelper::IsCosmicRay(MCParticleHelper::GetMainMCParticle(pCaloHit)))
                outputCaloHitList.push_back(pCaloHit);
        }
        catch (const StatusCodeException &)
        {
            std::cout << "CheatingRemovingCosmicRays::Run - Unable to determine MCParticle origin for an input CaloHit, which will be skipped." << std::endl;
            continue;
        }
    }

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, outputCaloHitList, m_outputCaloHitListName));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<CaloHit>(*this, m_outputCaloHitListName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingRemovingCosmicRays::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListName", m_inputCaloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputCaloHitListName", m_outputCaloHitListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
