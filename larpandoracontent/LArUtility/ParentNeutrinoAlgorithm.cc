/**
 *  @file   larpandoracontent/LArUtility/ParentNeutrinoAlgorithm.cc
 *
 *  @brief  Implementation of the parent neutrino algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArUtility/ParentNeutrinoAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode ParentNeutrinoAlgorithm::Run()
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
        algorithms.insert(algorithms.end(), m_twoDMopUpAlgorithms.begin(), m_twoDMopUpAlgorithms.end());
        algorithms.insert(algorithms.end(), m_threeDHitAlgorithms.begin(), m_threeDHitAlgorithms.end());
        algorithms.insert(algorithms.end(), m_threeDMopUpAlgorithms.begin(), m_threeDMopUpAlgorithms.end());
        algorithms.insert(algorithms.end(), m_neutrinoAlgorithms.begin(), m_neutrinoAlgorithms.end());
        algorithms.insert(algorithms.end(), m_listMovingAlgorithm);

        for (const std::string &algorithmName : algorithms)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunDaughterAlgorithm(*this, algorithmName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParentNeutrinoAlgorithm::PerformSlicing(SliceList &sliceList) const
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

StatusCode ParentNeutrinoAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, xmlHandle,
        "TwoDClustering", m_clusteringAlgorithm));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle,
        "TwoDAlgorithms", m_twoDAlgorithms));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle,
        "ThreeDAlgorithms", m_threeDAlgorithms));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle,
        "ThreeDHitAlgorithms", m_threeDHitAlgorithms));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle,
        "VertexAlgorithms", m_vertexAlgorithms));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle,
        "TwoDMopUpAlgorithms", m_twoDMopUpAlgorithms));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle,
        "ThreeDMopUpAlgorithms", m_threeDMopUpAlgorithms));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle,
        "NeutrinoAlgorithms", m_neutrinoAlgorithms));

    return ParentSlicingBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
