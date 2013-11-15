/**
 *  @file   LArContent/src/LArUtility/EventPreparationAlgorithm.cc
 * 
 *  @brief  Implementation of the event preparation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArUtility/EventPreparationAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode EventPreparationAlgorithm::Run()
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, this->ProcessMCParticles());
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, this->ProcessCaloHits());

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventPreparationAlgorithm::ProcessMCParticles()
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

StatusCode EventPreparationAlgorithm::ProcessCaloHits()
{
    // Split input calo hit list into different views
    const CaloHitList *pCaloHitList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentCaloHitList(*this, pCaloHitList));

    if (pCaloHitList->empty())
        return STATUS_CODE_NOT_INITIALIZED;

    CaloHitList caloHitListU, caloHitListV, caloHitListW;

    for (CaloHitList::const_iterator hitIter = pCaloHitList->begin(), hitIterEnd = pCaloHitList->end(); hitIter != hitIterEnd; ++hitIter)
    {
        CaloHit *pCaloHit = *hitIter;

        if (pCaloHit->GetMipEquivalentEnergy() < m_mipEquivalentCut)
            continue;

        if (VIEW_U == (*hitIter)->GetHitType())
        {
            if (!caloHitListU.insert(*hitIter).second)
                return STATUS_CODE_ALREADY_PRESENT;
        }
        else if (VIEW_V == (*hitIter)->GetHitType())
        {
            if (!caloHitListV.insert(*hitIter).second)
                return STATUS_CODE_ALREADY_PRESENT;
        }
        else if (VIEW_W == (*hitIter)->GetHitType())
        {
            if (!caloHitListW.insert(*hitIter).second)
                return STATUS_CODE_ALREADY_PRESENT;
        }
    }

    // If multiple hits in same location, take most energetic
    CaloHitList finalCaloHitListU, finalCaloHitListV, finalCaloHitListW;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->FinalizeHitList(caloHitListU, finalCaloHitListU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->FinalizeHitList(caloHitListV, finalCaloHitListV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->FinalizeHitList(caloHitListW, finalCaloHitListW));

    CaloHitList filteredInputList;
    filteredInputList.insert(finalCaloHitListU.begin(), finalCaloHitListU.end());
    filteredInputList.insert(finalCaloHitListV.begin(), finalCaloHitListV.end());
    filteredInputList.insert(finalCaloHitListW.begin(), finalCaloHitListW.end());

    // Save the lists
    if (!filteredInputList.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveCaloHitList(*this, filteredInputList, m_filteredCaloHitListName));

    if (!finalCaloHitListU.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveCaloHitList(*this, finalCaloHitListU, m_outputCaloHitListNameU));

    if (!finalCaloHitListV.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveCaloHitList(*this, finalCaloHitListV, m_outputCaloHitListNameV));

    if (!finalCaloHitListW.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveCaloHitList(*this, finalCaloHitListW, m_outputCaloHitListNameW));

    if (!m_currentCaloHitListReplacement.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentCaloHitList(*this, m_currentCaloHitListReplacement));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventPreparationAlgorithm::FinalizeHitList(const CaloHitList &originalList, CaloHitList &finalList) const
{
    if (originalList.empty())
        return STATUS_CODE_SUCCESS;

    OrderedCaloHitList orderedCaloHitList;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, orderedCaloHitList.Add(originalList));

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        CaloHitList *pLayerHitList = iter->second;

        for (CaloHitList::const_iterator iterI = pLayerHitList->begin(), iterIEnd = pLayerHitList->end(); iterI != iterIEnd; ++iterI)
        {
            CaloHit *pCaloHitI = *iterI;
            bool useCaloHit(true);

            for (CaloHitList::const_iterator iterJ = pLayerHitList->begin(), iterJEnd = pLayerHitList->end(); iterJ != iterJEnd; ++iterJ)
            {
                CaloHit *pCaloHitJ = *iterJ;

                if ((pCaloHitI->GetMipEquivalentEnergy() < pCaloHitJ->GetMipEquivalentEnergy()) &&
                    ((pCaloHitI->GetPositionVector() - pCaloHitJ->GetPositionVector()).GetMagnitudeSquared() < m_minCaloHitSeparationSquared))
                {
                    useCaloHit = false;
                    break;
                }
            }

            if (useCaloHit)
                finalList.insert(pCaloHitI);
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventPreparationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_mipEquivalentCut = 0.25f; // mips
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MipEquivalentCut", m_mipEquivalentCut));

    float minCaloHitSeparation = 0.4f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCaloHitSeparation", minCaloHitSeparation));
    m_minCaloHitSeparationSquared = minCaloHitSeparation * minCaloHitSeparation;

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

    m_currentCaloHitListReplacement = "caloHitListW";
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

    m_currentMCParticleListReplacement = "MCParticleList3D";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CurrentMCParticleListReplacement", m_currentMCParticleListReplacement));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
