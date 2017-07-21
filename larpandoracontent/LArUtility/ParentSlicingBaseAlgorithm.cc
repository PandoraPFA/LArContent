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

void ParentSlicingBaseAlgorithm::PerformSlicing(SliceList &sliceList) const
{
    this->FastReconstruction();
    m_pSlicingTool->Slice(this, m_caloHitListNames, m_clusterListNames, sliceList);
    this->RunAlgorithm(m_slicingListDeletionAlgorithm);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParentSlicingBaseAlgorithm::SaveTwoDCaloHitLists(const Slice &slice, const std::string &sliceIndexString) const
{
    for (const HitType hitType : m_hitTypeList)
    {
        const CaloHitList &caloHitList((TPC_VIEW_U == hitType) ? slice.m_caloHitListU : (TPC_VIEW_V == hitType) ? slice.m_caloHitListV : slice.m_caloHitListW);

        const std::string workingCaloHitListName(m_caloHitListNames.at(hitType) + sliceIndexString);
        const CaloHitList *pWorkingCaloHitList(nullptr);

        if (!sliceIndexString.empty() && (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, workingCaloHitListName, pWorkingCaloHitList)))
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, caloHitList, workingCaloHitListName));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParentSlicingBaseAlgorithm::RunTwoDClustering(const std::string &sliceIndexString, const std::string &clusteringAlgName,
    const bool existingClusterList, const StringVector &additionalTwoDAlgorithms) const
{
    for (const HitType hitType : m_hitTypeList)
    {
        const StatusCode listStatusCode(PandoraContentApi::ReplaceCurrentList<CaloHit>(*this, m_caloHitListNames.at(hitType) + sliceIndexString));

        if (STATUS_CODE_NOT_FOUND == listStatusCode)
            continue;

        if (STATUS_CODE_SUCCESS != listStatusCode)
            throw StatusCodeException(listStatusCode);

        std::string clusterListName;
        const ClusterList *pClusterList(nullptr);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunClusteringAlgorithm(*this, clusteringAlgName, pClusterList, clusterListName));

        if (pClusterList->empty())
        {
            if (!existingClusterList || (STATUS_CODE_SUCCESS != PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_clusterListNames.at(hitType))))
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::DropCurrentList<Cluster>(*this));
        }
        else
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, m_clusterListNames.at(hitType)));
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_clusterListNames.at(hitType)));
        }

        this->RunAlgorithms(additionalTwoDAlgorithms);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ParentSlicingBaseAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "ShouldPerformSlicing", m_shouldPerformSlicing));

    if (m_shouldPerformSlicing)
    {
        AlgorithmTool *pAlgorithmTool(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle,
            "Slicing", pAlgorithmTool));

        m_pSlicingTool = dynamic_cast<SlicingTool*>(pAlgorithmTool);

        if (!m_pSlicingTool)
            return STATUS_CODE_INVALID_PARAMETER;

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, xmlHandle,
            "SlicingListDeletion", m_slicingListDeletionAlgorithm));
    }

    return ParentBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
