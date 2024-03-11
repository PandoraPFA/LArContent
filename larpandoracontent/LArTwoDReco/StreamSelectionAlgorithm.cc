/**
 *  @file   larpandoradlcontent/LArTwoDReco/StreamSelectionAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning track shower cluster streaming algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArTwoDReco/StreamSelectionAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include <numeric>

using namespace pandora;

namespace lar_content
{

StreamSelectionAlgorithm::StreamSelectionAlgorithm() :
    m_inputListName{""},
    m_listType{"cluster"}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode StreamSelectionAlgorithm::Run()
{
    if (m_listNames.empty())
    {
        std::cout << "StreamSelectionAlgorithm::Run - Error: No output lists found" << std::endl;
        return STATUS_CODE_NOT_FOUND;
    }

    const ClusterList *pClusterList{nullptr};
    std::string originalClusterListName;
    if (m_inputListName.empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*this, originalClusterListName));
    }
    else
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputListName, pClusterList));
        originalClusterListName = m_inputListName;
    }

    for (std::string listName : m_listNames)
        m_clusterListMap[listName] = ClusterList();

    for (const Cluster *pCluster : *pClusterList)
    {
        StatusCode status{this->AllocateToStreams(pCluster)};
        if (status != STATUS_CODE_SUCCESS)
            return status;
    }

    // ATTN - We're ok with saving empty lists here and allowing future algorithms to simply do nothing if there are no clusters
    // Moves the subset of clusters in the cluster list from the old list to the new list
    for (std::string listName : m_listNames)
        PandoraContentApi::SaveList<ClusterList>(*this, originalClusterListName, listName, m_clusterListMap.at(listName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode StreamSelectionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InputListName", m_inputListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ListType", m_listType));
    std::transform(m_listType.begin(), m_listType.end(), m_listType.begin(), ::tolower);
    if (m_listType != "cluster")
    {
        std::cout << "StreamingAlgorithm::ReadSettings - Error: Only Cluster list type is supported at this time" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    // ATTN - Classes implementing this interface should ensure that m_listNames gets populated from the list names provided to the
    // implementing class
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "ListNames", m_listNames));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
