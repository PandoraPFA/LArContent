/**
 *  @file   larpandoracontent/LArUtility/ListPreparationAlgorithm.cc
 *
 *  @brief  Implementation of the list preparation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
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
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputCaloHitListName, pCaloHitList));

    if (pCaloHitList->empty())
        return;

    const bool checkMC(!m_selectNeutrinos || !m_selectCosmics);
    CaloHitList selectedCaloHitListU, selectedCaloHitListV, selectedCaloHitListW;

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        if (checkMC)
        {
            try
            {
                const bool isSelected(LArMCParticleHelper::IsNeutrinoInduced(pCaloHit) ? m_selectNeutrinos : m_selectCosmics);

                if (!isSelected)
                    continue;
            }
            catch (StatusCodeException &) {}
        }

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

        if ((pCaloHit->GetCellLengthScale() < m_minCellLengthScale) || (pCaloHit->GetCellLengthScale() > m_maxCellLengthScale))
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            {
                std::cout << "ListPreparationAlgorithm: found a hit with extent " << pCaloHit->GetCellLengthScale()
                          << ", require (" << m_minCellLengthScale << " - " << m_maxCellLengthScale << "), will remove it" << std::endl;
            }

            continue;
        }

        if (TPC_VIEW_U == pCaloHit->GetHitType())
        {
            selectedCaloHitListU.push_back(pCaloHit);
        }
        else if (TPC_VIEW_V == pCaloHit->GetHitType())
        {
            selectedCaloHitListV.push_back(pCaloHit);
        }
        else if (TPC_VIEW_W == pCaloHit->GetHitType())
        {
            selectedCaloHitListW.push_back(pCaloHit);
        }
    }

    CaloHitList filteredCaloHitListU, filteredCaloHitListV, filteredCaloHitListW;
    this->GetFilteredCaloHitList(selectedCaloHitListU, filteredCaloHitListU);
    this->GetFilteredCaloHitList(selectedCaloHitListV, filteredCaloHitListV);
    this->GetFilteredCaloHitList(selectedCaloHitListW, filteredCaloHitListW);

    CaloHitList filteredInputList;
    filteredInputList.insert(filteredInputList.end(), filteredCaloHitListU.begin(), filteredCaloHitListU.end());
    filteredInputList.insert(filteredInputList.end(), filteredCaloHitListV.begin(), filteredCaloHitListV.end());
    filteredInputList.insert(filteredInputList.end(), filteredCaloHitListW.begin(), filteredCaloHitListW.end());

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
    HitKDTree2D kdTree;
    HitKDNode2DList hitKDNode2DList;

    KDTreeBox hitsBoundingRegion2D = fill_and_bound_2d_kd_tree(inputList, hitKDNode2DList);
    kdTree.build(hitKDNode2DList, hitsBoundingRegion2D);

    // Remove hits that are in the same physical location!
    for (const CaloHit *const pCaloHit1 : inputList)
    {
        bool isUnique(true);
        KDTreeBox searchRegionHits(build_2d_kd_search_region(pCaloHit1, m_searchRegion1D, m_searchRegion1D));

        HitKDNode2DList found;
        kdTree.search(searchRegionHits, found);

        for (const auto &hit : found)
        {
            const CaloHit *const pCaloHit2(hit.data);

            if (pCaloHit1 == pCaloHit2)
                continue;

            const float displacementSquared((pCaloHit2->GetPositionVector() - pCaloHit1->GetPositionVector()).GetMagnitudeSquared());

            if (displacementSquared < std::numeric_limits<float>::epsilon())
            {
                const float deltaMip(pCaloHit2->GetMipEquivalentEnergy() > pCaloHit1->GetMipEquivalentEnergy());

                if ((deltaMip > std::numeric_limits<float>::epsilon()) ||
                    ((std::fabs(deltaMip) < std::numeric_limits<float>::epsilon()) && (outputList.end() != std::find(outputList.begin(), outputList.end(), pCaloHit2))))
                {
                    isUnique = false;
                    break;
                }
            }
        }

        if (isUnique)
        {
            outputList.push_back(pCaloHit1);
        }
        else
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "ListPreparationAlgorithm: found two hits in same location, will remove lowest pulse height" << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ListPreparationAlgorithm::ProcessMCParticles()
{
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputMCParticleListName, pMCParticleList));

    if (pMCParticleList->empty())
        return;

    const bool checkMC(!m_selectNeutrinos || !m_selectCosmics);
    MCParticleList mcParticleListU, mcParticleListV, mcParticleListW, mcParticleList3D;

    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        if (checkMC)
        {
            const bool isSelected(LArMCParticleHelper::IsNeutrinoInduced(pMCParticle) ? m_selectNeutrinos : m_selectCosmics);

            if (!isSelected)
                continue;
        }

        if (MC_VIEW_U == pMCParticle->GetMCParticleType())
        {
            mcParticleListU.push_back(pMCParticle);
        }
        else if (MC_VIEW_V == pMCParticle->GetMCParticleType())
        {
            mcParticleListV.push_back(pMCParticle);
        }
        else if (MC_VIEW_W == pMCParticle->GetMCParticleType())
        {
            mcParticleListW.push_back(pMCParticle);
        }
        else if (MC_3D == pMCParticle->GetMCParticleType())
        {
            mcParticleList3D.push_back(pMCParticle);
        }
    }

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

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputMCParticleListNameU", m_outputMCParticleListNameU));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputMCParticleListNameV", m_outputMCParticleListNameV));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputMCParticleListNameW", m_outputMCParticleListNameW));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
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
