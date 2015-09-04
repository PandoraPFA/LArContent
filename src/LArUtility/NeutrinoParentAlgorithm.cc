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
    // Initial reconstruction pass
    for (const HitType hitType : m_hitTypeList)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<CaloHit>(*this, m_caloHitListNames.at(hitType)));

        std::string clusterListName;
        const ClusterList *pClusterList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunClusteringAlgorithm(*this, m_clusteringAlgorithm, pClusterList, clusterListName));

        if (pClusterList->empty())
            continue;

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, m_clusterListNames.at(hitType)));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_clusterListNames.at(hitType)));

        for (const std::string &algorithmName : m_twoDAlgorithms)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunDaughterAlgorithm(*this, algorithmName));
    }

    StringVector preSlicingAlgorithms;
    preSlicingAlgorithms.insert(preSlicingAlgorithms.end(), m_threeDAlgorithms.begin(), m_threeDAlgorithms.end());
    preSlicingAlgorithms.insert(preSlicingAlgorithms.end(), m_threeDHitAlgorithms.begin(), m_threeDHitAlgorithms.end());

    for (const std::string &algorithmName : preSlicingAlgorithms)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunDaughterAlgorithm(*this, algorithmName));

    // Slicing - TODO: work out interface
    SliceList sliceList;
    // PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunDaughterAlgorithm(*this, m_slicingAlgorithm));

// Temporary hack
sliceList.push_back(Slice());
Slice &sliceHack(sliceList.at(0));
const CaloHitList *pTempListU(nullptr);
PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListNameU, pTempListU));
sliceHack.m_caloHitListU = *pTempListU;
const CaloHitList *pTempListV(nullptr);
PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListNameV, pTempListV));
sliceHack.m_caloHitListV = *pTempListV;
const CaloHitList *pTempListW(nullptr);
PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListNameW, pTempListW));
sliceHack.m_caloHitListW = *pTempListW;

    // Delete all existing algorithm objects and process each slice separately
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunDaughterAlgorithm(*this, m_listDeletionAlgorithm));

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
                continue;

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

    //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, xmlHandle,
    //    "Slicing", m_slicingAlgorithm));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, xmlHandle,
        "ListDeletion", m_listDeletionAlgorithm));

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
