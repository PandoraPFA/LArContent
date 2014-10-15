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

namespace lar_content
{

ListPreparationAlgorithm::ListPreparationAlgorithm() :
    m_mipEquivalentCut(std::numeric_limits<float>::epsilon()),
    m_minCellLengthScale(std::numeric_limits<float>::epsilon()),
    m_maxCellLengthScale(2.f),
    m_onlyAvailableCaloHits(true),
    m_inputCaloHitListName("Input"),
    m_mcNeutrinoSelection(false),
    m_inputMCParticleListName("Input")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ListPreparationAlgorithm::Run()
{
    if (!this->GetPandora().GetSettings()->SingleHitTypeClusteringMode())
    {
        std::cout << "ListPreparationAlgorithm: expect Pandora to be configured in SingleHitTypeClusteringMode." << std::endl;
        return STATUS_CODE_FAILURE;
    }

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
            std::cout << "ListPreparationAlgorithm: could not replace current calo hit list with list named: " << m_currentCaloHitListReplacement << std::endl;
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
            std::cout << "ListPreparationAlgorithm: could not replace current MC particle list with list named: " << m_currentMCParticleListReplacement << std::endl;
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

    CaloHitList selectedCaloHitListU, selectedCaloHitListV, selectedCaloHitListW;

    for (CaloHitList::const_iterator hitIter = pCaloHitList->begin(), hitIterEnd = pCaloHitList->end(); hitIter != hitIterEnd; ++hitIter)
    {
        CaloHit *pCaloHit = *hitIter;

        if (m_onlyAvailableCaloHits && !PandoraContentApi::IsAvailable(*this, pCaloHit))
            continue;

        if (pCaloHit->GetMipEquivalentEnergy() < m_mipEquivalentCut)
            continue;

        if ((pCaloHit->GetCellLengthScale() < m_minCellLengthScale) || (pCaloHit->GetCellLengthScale() > m_maxCellLengthScale))
        {
            std::cout << "ListPreparationAlgorithm: found a hit with extent " << pCaloHit->GetCellLengthScale() << ", will remove it" << std::endl;
            continue;
        }

        if (pCaloHit->GetInputEnergy() < std::numeric_limits<float>::epsilon())
        {
            std::cout << "ListPreparationAlgorithm: found a hit with zero energy, will remove it" << std::endl;
            continue;
        }

        if (TPC_VIEW_U == (*hitIter)->GetHitType())
        {
            if (!selectedCaloHitListU.insert(*hitIter).second)
                throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
        }
        else if (TPC_VIEW_V == (*hitIter)->GetHitType())
        {
            if (!selectedCaloHitListV.insert(*hitIter).second)
                throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
        }
        else if (TPC_VIEW_W == (*hitIter)->GetHitType())
        {
            if (!selectedCaloHitListW.insert(*hitIter).second)
                throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
        }
    }

    // Filter the selected hits
    CaloHitList filteredCaloHitListU, filteredCaloHitListV, filteredCaloHitListW;
    this->GetFilteredCaloHitList(selectedCaloHitListU, filteredCaloHitListU);
    this->GetFilteredCaloHitList(selectedCaloHitListV, filteredCaloHitListV);
    this->GetFilteredCaloHitList(selectedCaloHitListW, filteredCaloHitListW);

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

void ListPreparationAlgorithm::GetFilteredCaloHitList(const CaloHitList &inputList, CaloHitList &outputList)
{
    // Remove hits that are in the same physical location!
    for (CaloHitList::const_iterator hIter1 = inputList.begin(), hIterEnd1 = inputList.end(); hIter1 != hIterEnd1; ++hIter1)
    {
        CaloHit *pCaloHit1 = *hIter1;

        bool isUnique(true);

        for (CaloHitList::const_iterator hIter2 = inputList.begin(), hIterEnd2 = inputList.end(); hIter2 != hIterEnd2; ++hIter2)
        {
            CaloHit *pCaloHit2 = *hIter2;

            if (pCaloHit1 == pCaloHit2)
                continue;

            if ((pCaloHit2->GetPositionVector() - pCaloHit1->GetPositionVector()).GetMagnitudeSquared() < std::numeric_limits<float>::epsilon() &&
                pCaloHit2->GetMipEquivalentEnergy() > pCaloHit1->GetMipEquivalentEnergy())
            {
                isUnique = false;
                break;
            }
        }

        if (isUnique)
        {
            outputList.insert(pCaloHit1);
        }
        else
        {   
            std::cout << "ListPreparationAlgorithm: found two hits in same location, will remove lowest pulse height" << std::endl;
        }
    }

    // TODO Could merge these hits instead. Could also chop up long hits around this point
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
        if (m_mcNeutrinoSelection && !LArMCParticleHelper::GetParentNeutrinoId(*mcIter))
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
        else if (MC_3D == (*mcIter)->GetMCParticleType())
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
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MipEquivalentCut", m_mipEquivalentCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCellLengthScale", m_minCellLengthScale));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxCellLengthScale", m_maxCellLengthScale));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OnlyAvailableCaloHits", m_onlyAvailableCaloHits));

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

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "CurrentCaloHitListReplacement", m_currentCaloHitListReplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MCNeutrinoSelection", m_mcNeutrinoSelection));

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

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CurrentMCParticleListReplacement", m_currentMCParticleListReplacement));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
