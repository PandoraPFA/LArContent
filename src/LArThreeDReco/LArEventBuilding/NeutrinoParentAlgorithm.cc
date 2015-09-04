/**
 *  @file   LArContent/src/LArThreeDReco/LArEventBuilding/NeutrinoParentAlgorithm.cc
 * 
 *  @brief  Implementation of the neutrino parent algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArThreeDReco/LArEventBuilding/NeutrinoParentAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode NeutrinoParentAlgorithm::Run()
{
    for (const std::string &inputClusterListName : m_inputClusterListNames)
    {
        const StatusCode statusCode(PandoraContentApi::SaveList<Cluster>(*this, inputClusterListName, m_outputClusterListName));

        if (STATUS_CODE_SUCCESS != statusCode)
            std::cout << "NeutrinoParentAlgorithm: input cluster list not available " << inputClusterListName << std::endl;
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_outputClusterListName));

    // TODO - run algorithm to collect together separate cluster lists for each neutrino interaction
    // TODO - consider the precise interface for this
    const ClusterList *pAllClusterList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pAllClusterList));

    if (!pAllClusterList || pAllClusterList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "NeutrinoParentAlgorithm: no input clusters available " << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    ClusterCollectionList clusterCollectionList;
    clusterCollectionList.push_back(*pAllClusterList);

    for (const ClusterList &singleNeutrinoClusters : clusterCollectionList)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_outputClusterListName, m_workingClusterListName, singleNeutrinoClusters));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_workingClusterListName));

        for (const std::string &reconstructionAlgName : m_reconstructionAlgNames)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunDaughterAlgorithm(*this, reconstructionAlgName));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunDaughterAlgorithm(*this, m_listManagementAlgName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode NeutrinoParentAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "InputClusterListNames", m_inputClusterListNames));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputClusterListName", m_outputClusterListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "WorkingClusterListName", m_workingClusterListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle,
        "ReconstructionAlgorithms", m_reconstructionAlgNames));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, xmlHandle,
        "ListManagement", m_listManagementAlgName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
