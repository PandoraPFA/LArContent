/**
 *  @file   LArContent/src/LArUtility/ListPreparationAlgorithm.cc
 * 
 *  @brief  Implementation of the event preparation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArMCParticleHelper.h"

#include "LArUtility/ListPreparationAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode ListPreparationAlgorithm::Run()
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, this->ProcessCaloHits());
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, this->ProcessMCParticles());

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ListPreparationAlgorithm::ProcessCaloHits()
{
    // Split input calo hit list into different views
    const CaloHitList *pCaloHitList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputCaloHitListName, pCaloHitList));

    if (pCaloHitList->empty())
        return STATUS_CODE_NOT_INITIALIZED;

    CaloHitList filteredCaloHitListU, filteredCaloHitListV, filteredCaloHitListW;

    for (CaloHitList::const_iterator hitIter = pCaloHitList->begin(), hitIterEnd = pCaloHitList->end(); hitIter != hitIterEnd; ++hitIter)
    {
        CaloHit *pCaloHit = *hitIter;

        if (m_onlyAvailableCaloHits && !PandoraContentApi::IsAvailable(*this, pCaloHit))
            continue;

        if (pCaloHit->GetMipEquivalentEnergy() < m_mipEquivalentCut)
            continue;

        if (VIEW_U == (*hitIter)->GetHitType())
        {
            if (!filteredCaloHitListU.insert(*hitIter).second)
                return STATUS_CODE_ALREADY_PRESENT;
        }
        else if (VIEW_V == (*hitIter)->GetHitType())
        {
            if (!filteredCaloHitListV.insert(*hitIter).second)
                return STATUS_CODE_ALREADY_PRESENT;
        }
        else if (VIEW_W == (*hitIter)->GetHitType())
        {
            if (!filteredCaloHitListW.insert(*hitIter).second)
                return STATUS_CODE_ALREADY_PRESENT;
        }
    }

    // Group together views into overall lists
    CaloHitList filteredInputList;
    filteredInputList.insert(filteredCaloHitListU.begin(), filteredCaloHitListU.end());
    filteredInputList.insert(filteredCaloHitListV.begin(), filteredCaloHitListV.end());
    filteredInputList.insert(filteredCaloHitListW.begin(), filteredCaloHitListW.end());

    // Save the lists
    if (!filteredInputList.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, filteredInputList, m_filteredCaloHitListName));

    if (!filteredCaloHitListU.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, filteredCaloHitListU, m_outputCaloHitListNameU));

    if (!filteredCaloHitListV.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, filteredCaloHitListV, m_outputCaloHitListNameV));

    if (!filteredCaloHitListW.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, filteredCaloHitListW, m_outputCaloHitListNameW));

    if (!m_currentCaloHitListReplacement.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<CaloHit>(*this, m_currentCaloHitListReplacement));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ListPreparationAlgorithm::ProcessMCParticles()
{
    // Split input MC particles into different views
    const MCParticleList *pMCParticleList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputMCParticleListName, pMCParticleList));

    if (pMCParticleList->empty())
        return STATUS_CODE_NOT_INITIALIZED;

    MCParticleList mcParticleListU, mcParticleListV, mcParticleListW, mcParticleList3D;

    for (MCParticleList::const_iterator mcIter = pMCParticleList->begin(), mcIterEnd = pMCParticleList->end(); mcIter != mcIterEnd; ++mcIter)
    {
        if (m_mcNeutrinoSelection && !LArMCParticleHelper::GetPrimaryNeutrino(*mcIter))
            continue;

        if (MC_VIEW_U == (*mcIter)->GetMCParticleType())
        {
            if (!mcParticleListU.insert(*mcIter).second)
                return STATUS_CODE_ALREADY_PRESENT;
        }
        else if (MC_VIEW_V == (*mcIter)->GetMCParticleType())
        {
            if (!mcParticleListV.insert(*mcIter).second)
                return STATUS_CODE_ALREADY_PRESENT;
        }
        else if (MC_VIEW_W == (*mcIter)->GetMCParticleType())
        {
            if (!mcParticleListW.insert(*mcIter).second)
                return STATUS_CODE_ALREADY_PRESENT;
        }
        else if (MC_STANDARD == (*mcIter)->GetMCParticleType())
        {
            if (!mcParticleList3D.insert(*mcIter).second)
                return STATUS_CODE_ALREADY_PRESENT;
        }
    }

    // Save the lists
    if (!mcParticleListU.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, mcParticleListU, m_outputMCParticleListNameU));

    if (!mcParticleListV.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, mcParticleListV, m_outputMCParticleListNameV));

    if (!mcParticleListW.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, mcParticleListW, m_outputMCParticleListNameW));

    if (!mcParticleList3D.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, mcParticleList3D, m_outputMCParticleListName3D));

    if (!m_currentMCParticleListReplacement.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<MCParticle>(*this, m_currentMCParticleListReplacement));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ListPreparationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_mipEquivalentCut = 0.25f; // mips
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MipEquivalentCut", m_mipEquivalentCut));

    m_onlyAvailableCaloHits = true;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OnlyAvailableCaloHits", m_onlyAvailableCaloHits));

    m_inputCaloHitListName = "Input";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "InputCaloHitListName", m_inputCaloHitListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputCaloHitListNameU", m_outputCaloHitListNameU));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputCaloHitListNameV", m_outputCaloHitListNameV));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputCaloHitListNameW", m_outputCaloHitListNameW));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "FilteredCaloHitListName", m_filteredCaloHitListName));

    m_currentCaloHitListReplacement = "";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CurrentCaloHitListReplacement", m_currentCaloHitListReplacement));

    m_mcNeutrinoSelection = false;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MCNeutrinoSelection", m_mcNeutrinoSelection));

    m_inputMCParticleListName = "Input";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "InputMCParticleListName", m_inputMCParticleListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputMCParticleListNameU", m_outputMCParticleListNameU));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputMCParticleListNameV", m_outputMCParticleListNameV));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputMCParticleListNameW", m_outputMCParticleListNameW));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputMCParticleListName3D", m_outputMCParticleListName3D));

    m_currentMCParticleListReplacement = "";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CurrentMCParticleListReplacement", m_currentMCParticleListReplacement));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
