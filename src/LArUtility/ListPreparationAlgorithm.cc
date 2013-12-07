/**
 *  @file   LArContent/src/LArUtility/ListPreparationAlgorithm.cc
 * 
 *  @brief  Implementation of the event preparation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArUtility/ListPreparationAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode ListPreparationAlgorithm::Run()
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, this->ProcessMCParticles());
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, this->ProcessCaloHits());

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ListPreparationAlgorithm::ProcessMCParticles()
{
    // Split input MC particles into different views
    const MCParticleList *pMCParticleList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentMCParticleList(*this, pMCParticleList));

    if (pMCParticleList->empty())
        return STATUS_CODE_NOT_INITIALIZED;

    MCParticleList mcParticleListU, mcParticleListV, mcParticleListW, mcParticleList3D;

    for (MCParticleList::const_iterator mcIter = pMCParticleList->begin(), mcIterEnd = pMCParticleList->end(); mcIter != mcIterEnd; ++mcIter)
    {

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
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveMCParticleList(*this, mcParticleListU, m_outputMCParticleListNameU));

    if (!mcParticleListV.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveMCParticleList(*this, mcParticleListV, m_outputMCParticleListNameV));

    if (!mcParticleListW.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveMCParticleList(*this, mcParticleListW, m_outputMCParticleListNameW));

    if (!mcParticleList3D.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveMCParticleList(*this, mcParticleList3D, m_outputMCParticleListName3D));

    if (!m_currentMCParticleListReplacement.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentMCParticleList(*this, m_currentMCParticleListReplacement));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ListPreparationAlgorithm::ProcessCaloHits()
{
    // Split input calo hit list into different views
    const CaloHitList *pCaloHitList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentCaloHitList(*this, pCaloHitList));

    if (pCaloHitList->empty())
        return STATUS_CODE_NOT_INITIALIZED;

    CaloHitList filteredCaloHitListU, filteredCaloHitListV, filteredCaloHitListW;

    for (CaloHitList::const_iterator hitIter = pCaloHitList->begin(), hitIterEnd = pCaloHitList->end(); hitIter != hitIterEnd; ++hitIter)
    {
        CaloHit *pCaloHit = *hitIter;

        const bool isFiltered(pCaloHit->GetMipEquivalentEnergy() > m_mipEquivalentCut);

        if (VIEW_U == (*hitIter)->GetHitType())
        {
	    if (isFiltered)
	    {
                if (!filteredCaloHitListU.insert(*hitIter).second)
                    return STATUS_CODE_ALREADY_PRESENT;
	    }
            else
	    {
                // TODO: Persist the low pulse height hits  
	    }
        }
        else if (VIEW_V == (*hitIter)->GetHitType())
        {
            if (isFiltered)
	    {
                if (!filteredCaloHitListV.insert(*hitIter).second)
                    return STATUS_CODE_ALREADY_PRESENT;
	    }
            else
	    {
                // TODO: Persist the low pulse height hits   
	    }
        }
        else if (VIEW_W == (*hitIter)->GetHitType())
        {   
            if (isFiltered)
	    {
                if (!filteredCaloHitListW.insert(*hitIter).second)
                    return STATUS_CODE_ALREADY_PRESENT;
	    }
            else
	    {
	        // TODO: Persist the low pulse height hits
	    }
        }
    }

    // Group together views into overall lists
    CaloHitList filteredInputList;
    filteredInputList.insert(filteredCaloHitListU.begin(), filteredCaloHitListU.end());
    filteredInputList.insert(filteredCaloHitListV.begin(), filteredCaloHitListV.end());
    filteredInputList.insert(filteredCaloHitListW.begin(), filteredCaloHitListW.end());

    // Save the lists
    if (!filteredInputList.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveCaloHitList(*this, filteredInputList, m_filteredCaloHitListName));

    if (!filteredCaloHitListU.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveCaloHitList(*this, filteredCaloHitListU, m_outputCaloHitListNameU));

    if (!filteredCaloHitListV.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveCaloHitList(*this, filteredCaloHitListV, m_outputCaloHitListNameV));

    if (!filteredCaloHitListW.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveCaloHitList(*this, filteredCaloHitListW, m_outputCaloHitListNameW));

    if (!m_currentCaloHitListReplacement.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentCaloHitList(*this, m_currentCaloHitListReplacement));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ListPreparationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_mipEquivalentCut = 0.25f; // mips
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MipEquivalentCut", m_mipEquivalentCut));

    m_outputCaloHitListNameU = "caloHitListU";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputCaloHitListNameU", m_outputCaloHitListNameU));

    m_outputCaloHitListNameV = "caloHitListV";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputCaloHitListNameV", m_outputCaloHitListNameV));

    m_outputCaloHitListNameW = "caloHitListW";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputCaloHitListNameW", m_outputCaloHitListNameW));

    m_filteredCaloHitListName = "caloHitListFiltered";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FilteredCaloHitListName", m_filteredCaloHitListName));

    m_currentCaloHitListReplacement = "";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CurrentCaloHitListReplacement", m_currentCaloHitListReplacement));

    m_outputMCParticleListNameU = "MCParticleListU";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputMCParticleListNameU", m_outputMCParticleListNameU));

    m_outputMCParticleListNameV = "MCParticleListV";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputMCParticleListNameV", m_outputMCParticleListNameV));

    m_outputMCParticleListNameW = "MCParticleListW";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputMCParticleListNameW", m_outputMCParticleListNameW));

    m_outputMCParticleListName3D = "MCParticleList3D";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputMCParticleListName3D", m_outputMCParticleListName3D));

    m_currentMCParticleListReplacement = "";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CurrentMCParticleListReplacement", m_currentMCParticleListReplacement));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
