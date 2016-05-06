/**
 *  @file   larpandoracontent/LArUtility/ListPreparationAlgorithm.cc
 *
 *  @brief  Implementation of the list preparation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"
#include "larpandoracontent/LArUtility/ListPreparationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ListPreparationAlgorithm::ListPreparationAlgorithm() :
    m_mipEquivalentCut(std::numeric_limits<float>::epsilon()),
    m_minCellLengthScale(std::numeric_limits<float>::epsilon()),
    m_maxCellLengthScale(2.f),
    m_searchRegion1D(0.1f),
    m_onlyAvailableCaloHits(true),
    m_inputCaloHitListName("Input"),
    m_inputMCParticleListName("Input"),
    m_selectNeutrinos(true),
    m_selectCosmics(true)
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
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "ListPreparationAlgorithm: could not replace current calo hit list with list named: " << m_currentCaloHitListReplacement << std::endl;
        }
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
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "ListPreparationAlgorithm: could not replace current MC particle list with list named: " << m_currentMCParticleListReplacement << std::endl;
        }
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

    const bool checkMC(!m_selectNeutrinos || !m_selectCosmics);

    CaloHitList selectedCaloHitListU, selectedCaloHitListV, selectedCaloHitListW;

    for (CaloHitList::const_iterator hitIter = pCaloHitList->begin(), hitIterEnd = pCaloHitList->end(); hitIter != hitIterEnd; ++hitIter)
    {
        const CaloHit *const pCaloHit = *hitIter;

        if (m_onlyAvailableCaloHits && !PandoraContentApi::IsAvailable(*this, pCaloHit))
            continue;

        if (pCaloHit->GetMipEquivalentEnergy() < m_mipEquivalentCut)
            continue;

        if (pCaloHit->GetInputEnergy() < std::numeric_limits<float>::epsilon())
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "ListPreparationAlgorithm: found a hit with zero energy, will remove it" << std::endl;
            continue;
        }

        if (pCaloHit->GetCellLengthScale() < m_minCellLengthScale)
        {
            if (pCaloHit->GetCellLengthScale() < std::numeric_limits<float>::epsilon())
            {
                if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                    std::cout << "ListPreparationAlgorithm: found a hit with zero extent, will remove it" << std::endl;
            }
            else
            {
                if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                    std::cout << "ListPreparationAlgorithm: found a hit with extent " << pCaloHit->GetCellLengthScale()
                              << " (<" << m_minCellLengthScale << "), will remove it" << std::endl;
            }
            continue;
        }

        if (pCaloHit->GetCellLengthScale() > m_maxCellLengthScale)
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "ListPreparationAlgorithm: found a hit with extent " << pCaloHit->GetCellLengthScale()
                          << " (>" << m_maxCellLengthScale << "), will remove it" << std::endl;
            continue;
        }

        if (checkMC)
        {
            try
            {
                const bool isNeutrinoInduced(LArMCParticleHelper::IsNeutrinoInduced(pCaloHit));
                const bool isSelected(isNeutrinoInduced ? m_selectNeutrinos : m_selectCosmics);

                if (!isSelected)
                    continue;
            }
            catch (StatusCodeException &)
            {
            }
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
    // Initialize kd tree
    HitKDTree2D kdTree;
    HitKDNode2DList hitKDNode2DList;
    KDTreeBox hitsBoundingRegion2D = fill_and_bound_2d_kd_tree(this, inputList, hitKDNode2DList, true);
    kdTree.build(hitKDNode2DList, hitsBoundingRegion2D);

    // Remove hits that are in the same physical location!
    for (const CaloHit *const pCaloHit1 : inputList)
    {
        bool isUnique(true);

        // Get nearby hits from kd tree
        CaloHitList nearbyHits;
        KDTreeBox searchRegionHits = build_2d_kd_search_region(pCaloHit1, m_searchRegion1D, m_searchRegion1D);

        HitKDNode2DList found;
        kdTree.search(searchRegionHits, found);

        for (const auto &hit : found)
        {
            const CaloHit *const pCaloHit2 = hit.data;

            if (pCaloHit1 == pCaloHit2)
                continue;

            if (((pCaloHit2->GetPositionVector() - pCaloHit1->GetPositionVector()).GetMagnitudeSquared() < std::numeric_limits<float>::epsilon()) &&
                ((std::fabs(pCaloHit2->GetMipEquivalentEnergy() - pCaloHit1->GetMipEquivalentEnergy()) < std::numeric_limits<float>::epsilon()) ||
                (pCaloHit2->GetMipEquivalentEnergy() > pCaloHit1->GetMipEquivalentEnergy())) )
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
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
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

    const bool checkMC(!m_selectNeutrinos || !m_selectCosmics);

    MCParticleList mcParticleListU, mcParticleListV, mcParticleListW, mcParticleList3D;

    for (MCParticleList::const_iterator mcIter = pMCParticleList->begin(), mcIterEnd = pMCParticleList->end(); mcIter != mcIterEnd; ++mcIter)
    {
        if (checkMC)
        {
            const bool isNeutrinoInduced(LArMCParticleHelper::IsNeutrinoInduced(*mcIter));
            const bool isSelected(isNeutrinoInduced ? m_selectNeutrinos : m_selectCosmics);

            if (!isSelected)
                continue;
        }

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
        "SearchRegion1D", m_searchRegion1D));

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

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SelectNeutrinos", m_selectNeutrinos));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SelectCosmics", m_selectCosmics));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
