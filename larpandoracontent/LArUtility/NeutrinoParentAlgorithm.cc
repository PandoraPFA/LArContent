/**
 *  @file   LArContent/src/LArUtility/NeutrinoParentAlgorithm.cc
 * 
 *  @brief  Implementation of the neutrino parent algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArUtility/NeutrinoParentAlgorithm.h"

using namespace pandora;

namespace lar_content
{

NeutrinoParentAlgorithm::NeutrinoParentAlgorithm() :
    m_shouldPerformSlicing(true),
    m_pSlicingTool(nullptr)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode NeutrinoParentAlgorithm::Initialize()
{
    m_hitTypeList.push_back(TPC_VIEW_U);
    m_hitTypeList.push_back(TPC_VIEW_V);
    m_hitTypeList.push_back(TPC_VIEW_W);

    m_caloHitListNames[TPC_VIEW_U] = m_caloHitListNameU;
    m_caloHitListNames[TPC_VIEW_V] = m_caloHitListNameV;
    m_caloHitListNames[TPC_VIEW_W] = m_caloHitListNameW;

    m_clusterListNames[TPC_VIEW_U] = m_clusterListNameU;
    m_clusterListNames[TPC_VIEW_V] = m_clusterListNameV;
    m_clusterListNames[TPC_VIEW_W] = m_clusterListNameW;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode NeutrinoParentAlgorithm::Run()
{
    SliceList sliceList;

    if (m_shouldPerformSlicing)
    {
        this->PerformSlicing(sliceList);
    }
    else
    {
        this->CopyAllHitsToSingleSlice(sliceList);
    }

    unsigned int sliceCounter(0);

    for (const Slice &slice : sliceList)
    {
        for (const HitType hitType : m_hitTypeList)
        {
            const CaloHitList &caloHitList((TPC_VIEW_U == hitType) ? slice.m_caloHitListU : (TPC_VIEW_V == hitType) ? slice.m_caloHitListV : slice.m_caloHitListW);
            const std::string workingCaloHitListName(m_caloHitListNames.at(hitType) + TypeToString(sliceCounter++));

            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, caloHitList, workingCaloHitListName));
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<CaloHit>(*this, workingCaloHitListName));

            std::string clusterListName;
            const ClusterList *pClusterList(nullptr);
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunClusteringAlgorithm(*this, m_clusteringAlgorithm, pClusterList, clusterListName));

            if (pClusterList->empty())
            {
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::DropCurrentList<Cluster>(*this));
                continue;
            }

            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, m_clusterListNames.at(hitType)));
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_clusterListNames.at(hitType)));

            for (const std::string &algorithmName : m_twoDAlgorithms)
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunDaughterAlgorithm(*this, algorithmName));
        }

        StringVector algorithms;
        algorithms.insert(algorithms.end(), m_vertexAlgorithms.begin(), m_vertexAlgorithms.end());
        algorithms.insert(algorithms.end(), m_threeDAlgorithms.begin(), m_threeDAlgorithms.end());
        algorithms.insert(algorithms.end(), m_mopUpAlgorithms.begin(), m_mopUpAlgorithms.end());
        algorithms.insert(algorithms.end(), m_threeDHitAlgorithms.begin(), m_threeDHitAlgorithms.end());
        algorithms.insert(algorithms.end(), m_neutrinoAlgorithms.begin(), m_neutrinoAlgorithms.end());
        algorithms.insert(algorithms.end(), m_listMovingAlgorithm);

        for (const std::string &algorithmName : algorithms)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunDaughterAlgorithm(*this, algorithmName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoParentAlgorithm::PerformSlicing(SliceList &sliceList) const
{
    for (const HitType hitType : m_hitTypeList)
    {
        const StatusCode caloHitStatusCode(PandoraContentApi::ReplaceCurrentList<CaloHit>(*this, m_caloHitListNames.at(hitType)));

        if (STATUS_CODE_NOT_FOUND == caloHitStatusCode)
            continue;

        if (STATUS_CODE_SUCCESS != caloHitStatusCode)
            throw StatusCodeException(caloHitStatusCode);

        std::string clusterListName;
        const ClusterList *pClusterList(nullptr);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunClusteringAlgorithm(*this, m_clusteringAlgorithm, pClusterList, clusterListName));

        if (pClusterList->empty())
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::DropCurrentList<Cluster>(*this));
            continue;
        }

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, m_clusterListNames.at(hitType)));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_clusterListNames.at(hitType)));

        for (const std::string &algorithmName : m_twoDAlgorithms)
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunDaughterAlgorithm(*this, algorithmName));
    }

    StringVector preSlicingAlgorithms;
    preSlicingAlgorithms.insert(preSlicingAlgorithms.end(), m_threeDAlgorithms.begin(), m_threeDAlgorithms.end());
    preSlicingAlgorithms.insert(preSlicingAlgorithms.end(), m_threeDHitAlgorithms.begin(), m_threeDHitAlgorithms.end());

    for (const std::string &algorithmName : preSlicingAlgorithms)
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunDaughterAlgorithm(*this, algorithmName));

    // Slice the three dimensional clusters into separate, distinct interactions for reprocessing
    m_pSlicingTool->Slice(this, m_caloHitListNames, m_clusterListNames, sliceList);

    // Delete all existing algorithm objects and process each slice separately
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunDaughterAlgorithm(*this, m_listDeletionAlgorithm));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoParentAlgorithm::CopyAllHitsToSingleSlice(SliceList &sliceList) const
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

StatusCode NeutrinoParentAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "CaloHitListNameU", m_caloHitListNameU));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "CaloHitListNameV", m_caloHitListNameV));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "CaloHitListNameW", m_caloHitListNameW));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterListNameU", m_clusterListNameU));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterListNameV", m_clusterListNameV));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterListNameW", m_clusterListNameW));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, xmlHandle,
        "TwoDClustering", m_clusteringAlgorithm));

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
            "ListDeletion", m_listDeletionAlgorithm));
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, xmlHandle,
        "ListMoving", m_listMovingAlgorithm));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle,
        "TwoDAlgorithms", m_twoDAlgorithms));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle,
        "ThreeDAlgorithms", m_threeDAlgorithms));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle,
        "ThreeDHitAlgorithms", m_threeDHitAlgorithms));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle,
        "VertexAlgorithms", m_vertexAlgorithms));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle,
        "MopUpAlgorithms", m_mopUpAlgorithms));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle,
        "NeutrinoAlgorithms", m_neutrinoAlgorithms));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
