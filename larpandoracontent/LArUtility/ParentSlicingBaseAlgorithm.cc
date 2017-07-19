/**
 *  @file   larpandoracontent/LArUtility/ParentSlicingBaseAlgorithm.cc
 *
 *  @brief  Implementation of the parent slicing base algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArUtility/ParentNeutrinoAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ParentSlicingBaseAlgorithm::ParentSlicingBaseAlgorithm() :
    m_shouldPerformSlicing(true),
    m_pSlicingTool(nullptr)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

ParentSlicingBaseAlgorithm::~ParentSlicingBaseAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParentSlicingBaseAlgorithm::CopyAllHitsToSingleSlice(SliceList &sliceList) const
{
    if (!sliceList.empty())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const CaloHitList *pCaloHitListU(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this,
        m_caloHitListNames.at(TPC_VIEW_U), pCaloHitListU));

    const CaloHitList *pCaloHitListV(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this,
        m_caloHitListNames.at(TPC_VIEW_V), pCaloHitListV));

    const CaloHitList *pCaloHitListW(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this,
        m_caloHitListNames.at(TPC_VIEW_W), pCaloHitListW));

    if (pCaloHitListU || pCaloHitListV || pCaloHitListW)
    {
        sliceList.push_back(Slice());
        Slice &slice(sliceList.at(0));

        if (pCaloHitListU) slice.m_caloHitListU = *pCaloHitListU;
        if (pCaloHitListV) slice.m_caloHitListV = *pCaloHitListV;
        if (pCaloHitListW) slice.m_caloHitListW = *pCaloHitListW;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ParentSlicingBaseAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    if (m_shouldPerformSlicing)
    {
        AlgorithmTool *pAlgorithmTool(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle,
            "Slicing", pAlgorithmTool));

        m_pSlicingTool = dynamic_cast<SlicingTool*>(pAlgorithmTool);

        if (!m_pSlicingTool)
            return STATUS_CODE_INVALID_PARAMETER;

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, xmlHandle,
            "ListDeletion", m_listDeletionAlgorithm));
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, xmlHandle,
        "ListMoving", m_listMovingAlgorithm));

    return ParentBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
