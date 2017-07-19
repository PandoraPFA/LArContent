/**
 *  @file   larpandoracontent/LArUtility/ParentCosmicRayAlgorithm.cc
 *
 *  @brief  Implementation of the parent cosmic ray algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArUtility/ParentCosmicRayAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode ParentCosmicRayAlgorithm::Run()
{
    // Track oriented cosmic-ray muon reconstruction
    for (const HitType hitType : m_hitTypeList)
    {
        const StatusCode caloHitStatusCode(PandoraContentApi::ReplaceCurrentList<CaloHit>(*this, m_caloHitListNames.at(hitType)));

        if (STATUS_CODE_NOT_FOUND == caloHitStatusCode)
            continue;

        if (STATUS_CODE_SUCCESS != caloHitStatusCode)
            throw StatusCodeException(caloHitStatusCode);

        std::string clusterListName;
        const ClusterList *pClusterList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunClusteringAlgorithm(*this, m_trackClusteringAlgorithm, pClusterList, clusterListName));

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

    for (const std::string &algorithmName : m_threeDAlgorithms)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunDaughterAlgorithm(*this, algorithmName));

    // For clusters not in track particles, start again and form delta-ray clusters
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunDaughterAlgorithm(*this, m_listPruningAlgorithm));

    for (const HitType hitType : m_hitTypeList) // TODO refactor
    {
        const StatusCode caloHitStatusCode(PandoraContentApi::ReplaceCurrentList<CaloHit>(*this, m_caloHitListNames.at(hitType)));

        if (STATUS_CODE_NOT_FOUND == caloHitStatusCode)
            continue;

        if (STATUS_CODE_SUCCESS != caloHitStatusCode)
            throw StatusCodeException(caloHitStatusCode);

        std::string clusterListName;
        const ClusterList *pClusterList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunClusteringAlgorithm(*this, m_deltaRayClusteringAlgorithm, pClusterList, clusterListName));

        if (pClusterList->empty())
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_clusterListNames.at(hitType))); // NOTE, different to above
            continue;
        }

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, m_clusterListNames.at(hitType)));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_clusterListNames.at(hitType)));
    }

    // Complete the reconstruction
    StringVector algorithms;
    algorithms.insert(algorithms.end(), m_deltaRayAlgorithms.begin(), m_deltaRayAlgorithms.end());
    algorithms.insert(algorithms.end(), m_twoDRemnantAlgorithms.begin(), m_twoDRemnantAlgorithms.end());
    algorithms.insert(algorithms.end(), m_threeDRemnantAlgorithms.begin(), m_threeDRemnantAlgorithms.end());
    algorithms.insert(algorithms.end(), m_threeDHitAlgorithms.begin(), m_threeDHitAlgorithms.end());
    algorithms.insert(algorithms.end(), m_vertexAlgorithms.begin(), m_vertexAlgorithms.end());

    for (const std::string &algorithmName : algorithms)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunDaughterAlgorithm(*this, algorithmName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ParentCosmicRayAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, xmlHandle,
        "TwoDTrackClustering", m_trackClusteringAlgorithm));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, xmlHandle,
        "TwoDDeltaRayClustering", m_deltaRayClusteringAlgorithm));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, xmlHandle,
        "ListPruning", m_listPruningAlgorithm));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle,
        "TwoDAlgorithms", m_twoDAlgorithms));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle,
        "ThreeDAlgorithms", m_threeDAlgorithms));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle,
        "DeltaRayAlgorithms", m_deltaRayAlgorithms));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle,
        "TwoDRemnantAlgorithms", m_twoDRemnantAlgorithms));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle,
        "ThreeDRemnantAlgorithms", m_threeDRemnantAlgorithms));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle,
        "ThreeDHitAlgorithms", m_threeDHitAlgorithms));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle,
        "VertexAlgorithms", m_vertexAlgorithms));

    return ParentBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
