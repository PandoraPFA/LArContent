/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterCreation/ClusteringParentAlgorithm.cc
 *
 *  @brief  Implementation of the clustering parent algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArTwoDReco/LArClusterCreation/ClusteringParentAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ClusteringParentAlgorithm::ClusteringParentAlgorithm() :
    m_replaceCurrentCaloHitList(false),
    m_replaceCurrentClusterList(true)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusteringParentAlgorithm::Run()
{
    // If specified, change the current calo hit list, i.e. the input to the clustering algorithm
    std::string originalCaloHitListName;

    if (!m_inputCaloHitListName.empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<CaloHit>(*this, originalCaloHitListName));
        const StatusCode statusCode(PandoraContentApi::ReplaceCurrentList<CaloHit>(*this, m_inputCaloHitListName));

        if (STATUS_CODE_NOT_FOUND == statusCode)
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "ClusteringParentAlgorithm: calohit list not found " << m_inputCaloHitListName << std::endl;

            return STATUS_CODE_SUCCESS;
        }

        if (STATUS_CODE_SUCCESS != statusCode)
            return statusCode;
    }

    // Run the initial cluster formation algorithm
    const ClusterList *pClusterList = NULL;
    std::string newClusterListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        PandoraContentApi::RunClusteringAlgorithm(*this, m_clusteringAlgorithmName, pClusterList, newClusterListName));

    // Run the topological association algorithms to modify clusters
    if (!pClusterList->empty() && !m_associationAlgorithmName.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunDaughterAlgorithm(*this, m_associationAlgorithmName));

    // Save the new cluster list
    if (!pClusterList->empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, m_clusterListName));

        if (m_replaceCurrentClusterList)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_clusterListName));
    }

    // Unless specified, return current calo hit list to that when algorithm started
    if (!m_replaceCurrentCaloHitList && !m_inputCaloHitListName.empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<CaloHit>(*this, originalCaloHitListName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusteringParentAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    // Daughter algorithm parameters
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, xmlHandle, "ClusterFormation", m_clusteringAlgorithmName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ProcessAlgorithm(*this, xmlHandle, "ClusterAssociation", m_associationAlgorithmName));

    // Input parameters: name of input calo hit list and whether it should persist as the current list after algorithm has finished
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListName", m_inputCaloHitListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ReplaceCurrentCaloHitList", m_replaceCurrentCaloHitList));

    // Output parameters: name of output cluster list and whether it should subsequently be used as the current list
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListName", m_clusterListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ReplaceCurrentClusterList", m_replaceCurrentClusterList));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
