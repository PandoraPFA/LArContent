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
    try
    {
        this->ProcessCaloHits();
    }
    catch (StatusCodeException &)
    {
    }

    if (!m_currentCaloHitListReplacement.empty())
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::ReplaceCurrentList<CaloHit>(*this, m_currentCaloHitListReplacement))
            std::cout << "ListPreparationAlgorithm: Could not replace current calo hit list with list named: " << m_currentCaloHitListReplacement << std::endl;
    }

    try
    {
        this->ProcessMCParticles();
    }
    catch (StatusCodeException &)
    {
    }

    if (!m_currentMCParticleListReplacement.empty())
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::ReplaceCurrentList<MCParticle>(*this, m_currentMCParticleListReplacement))
            std::cout << "ListPreparationAlgorithm: Could not replace current MC particle list with list named: " << m_currentMCParticleListReplacement << std::endl;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ListPreparationAlgorithm::ProcessCaloHits()
{
    // Split input calo hit list into different views
    const CaloHitList *pCaloHitList = NULL;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputCaloHitListName, pCaloHitList));

    if (pCaloHitList->empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    CaloHitList filteredCaloHitListU, filteredCaloHitListV, filteredCaloHitListW;

    for (CaloHitList::const_iterator hitIter = pCaloHitList->begin(), hitIterEnd = pCaloHitList->end(); hitIter != hitIterEnd; ++hitIter)
    {
        CaloHit *pCaloHit = *hitIter;

        if (m_onlyAvailableCaloHits && !PandoraContentApi::IsAvailable(*this, pCaloHit))
            continue;

        if (pCaloHit->GetMipEquivalentEnergy() < m_mipEquivalentCut)
            continue;

        if (TPC_VIEW_U == (*hitIter)->GetHitType())
        {
            if (!filteredCaloHitListU.insert(*hitIter).second)
                throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
        }
        else if (TPC_VIEW_V == (*hitIter)->GetHitType())
        {
            if (!filteredCaloHitListV.insert(*hitIter).second)
                throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
        }
        else if (TPC_VIEW_W == (*hitIter)->GetHitType())
        {
            if (!filteredCaloHitListW.insert(*hitIter).second)
                throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
        }
    }

    // Group together views into overall lists
    CaloHitList filteredInputList;
    filteredInputList.insert(filteredCaloHitListU.begin(), filteredCaloHitListU.end());
    filteredInputList.insert(filteredCaloHitListV.begin(), filteredCaloHitListV.end());
    filteredInputList.insert(filteredCaloHitListW.begin(), filteredCaloHitListW.end());

    // Save the lists
    if (!filteredInputList.empty() && !m_filteredCaloHitListName.empty())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, filteredInputList, m_filteredCaloHitListName));

    if (!filteredCaloHitListU.empty() && !m_outputCaloHitListNameU.empty())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, filteredCaloHitListU, m_outputCaloHitListNameU));

    if (!filteredCaloHitListV.empty() && !m_outputCaloHitListNameV.empty())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, filteredCaloHitListV, m_outputCaloHitListNameV));

    if (!filteredCaloHitListW.empty() && !m_outputCaloHitListNameW.empty())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, filteredCaloHitListW, m_outputCaloHitListNameW));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ListPreparationAlgorithm::ProcessMCParticles()
{
    // Split input MC particles into different views
    const MCParticleList *pMCParticleList = NULL;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputMCParticleListName, pMCParticleList));

    if (pMCParticleList->empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    MCParticleList mcParticleListU, mcParticleListV, mcParticleListW, mcParticleList3D;

    for (MCParticleList::const_iterator mcIter = pMCParticleList->begin(), mcIterEnd = pMCParticleList->end(); mcIter != mcIterEnd; ++mcIter)
    {
        if (m_mcNeutrinoSelection && !LArMCParticleHelper::GetPrimaryNeutrino(*mcIter))
            continue;

        if (MC_VIEW_U == (*mcIter)->GetMCParticleType())
        {
            if (!mcParticleListU.insert(*mcIter).second)
                throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
        }
        else if (MC_VIEW_V == (*mcIter)->GetMCParticleType())
        {
            if (!mcParticleListV.insert(*mcIter).second)
                throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
        }
        else if (MC_VIEW_W == (*mcIter)->GetMCParticleType())
        {
            if (!mcParticleListW.insert(*mcIter).second)
                throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
        }
        else if (MC_STANDARD == (*mcIter)->GetMCParticleType())
        {
            if (!mcParticleList3D.insert(*mcIter).second)
                throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
        }
    }

    // Save the lists
    if (!mcParticleListU.empty() && !m_outputMCParticleListNameU.empty())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, mcParticleListU, m_outputMCParticleListNameU));

    if (!mcParticleListV.empty() && !m_outputMCParticleListNameV.empty())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, mcParticleListV, m_outputMCParticleListNameV));

    if (!mcParticleListW.empty() && !m_outputMCParticleListNameW.empty())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, mcParticleListW, m_outputMCParticleListNameW));

    if (!mcParticleList3D.empty() && !m_outputMCParticleListName3D.empty())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, mcParticleList3D, m_outputMCParticleListName3D));
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

    m_outputCaloHitListNameU.clear();
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputCaloHitListNameU", m_outputCaloHitListNameU));

    m_outputCaloHitListNameV.clear();
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputCaloHitListNameV", m_outputCaloHitListNameV));

    m_outputCaloHitListNameW.clear();
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputCaloHitListNameW", m_outputCaloHitListNameW));

    m_filteredCaloHitListName.clear();
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FilteredCaloHitListName", m_filteredCaloHitListName));

    m_currentCaloHitListReplacement.clear();
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CurrentCaloHitListReplacement", m_currentCaloHitListReplacement));

    m_mcNeutrinoSelection = false;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MCNeutrinoSelection", m_mcNeutrinoSelection));

    m_inputMCParticleListName = "Input";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "InputMCParticleListName", m_inputMCParticleListName));

    m_outputMCParticleListNameU.clear();
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputMCParticleListNameU", m_outputMCParticleListNameU));

    m_outputMCParticleListNameV.clear();
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputMCParticleListNameV", m_outputMCParticleListNameV));

    m_outputMCParticleListNameW.clear();
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputMCParticleListNameW", m_outputMCParticleListNameW));

    m_outputMCParticleListName3D.clear();
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputMCParticleListName3D", m_outputMCParticleListName3D));

    m_currentMCParticleListReplacement.clear();
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CurrentMCParticleListReplacement", m_currentMCParticleListReplacement));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
