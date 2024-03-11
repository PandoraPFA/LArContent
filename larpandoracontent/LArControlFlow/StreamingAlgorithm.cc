/**
 *  @file   larpandoradlcontent/LArTwoDReco/StreamingAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning track shower cluster streaming algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArControlFlow/StreamingAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include <algorithm>
#include <numeric>

using namespace pandora;

namespace lar_content
{

StreamingAlgorithm::StreamingAlgorithm() :
    m_listType{"cluster"}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StreamingAlgorithm::~StreamingAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode StreamingAlgorithm::Run()
{
    unsigned int i{0};
    for (std::string listName : m_inputListNames)
    {
        std::string algStreamName{"Algorithms" + listName};
        const ClusterList *pClusterList{nullptr};
        // Set the input list as current
        PandoraContentApi::ReplaceCurrentList<Cluster>(*this, listName);
        StatusCode code{PandoraContentApi::GetCurrentList(*this, pClusterList)};
        if (code == STATUS_CODE_SUCCESS)
        {
            for (const auto &alg : m_streamAlgorithmMap.at(algStreamName))
            { // ATTN - The algorithms replace the current list as they go
                PandoraContentApi::GetCurrentList(*this, pClusterList);
                if (!pClusterList->empty())
                {
                    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunDaughterAlgorithm(*this, alg));
                }
            }
            // Save the current list to the target output list
            if (!m_outputListName.empty())
                PandoraContentApi::SaveList<Cluster>(*this, m_outputListName);
            else
                PandoraContentApi::SaveList<Cluster>(*this, m_outputListNames.at(i));
        }
        else if (code != STATUS_CODE_NOT_INITIALIZED)
        {
            return code;
        }
        ++i;
    }

    // If we have a single output list specified, set that list as the current list
    if (!m_outputListName.empty())
        PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_outputListName);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode StreamingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ListType", m_listType));
    std::transform(m_listType.begin(), m_listType.end(), m_listType.begin(), ::tolower);
    if (m_listType != "cluster")
    {
        std::cout << "StreamingAlgorithm::ReadSettings - Error: Only Cluster list type is supported at this time" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputListNames", m_inputListNames));
    if (m_inputListNames.empty())
    {
        std::cout << "StreamingAlgorithm::ReadSettings - Error: No input lists found" << std::endl;
        return STATUS_CODE_NOT_FOUND;
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "OutputListName", m_outputListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "OutputListNames", m_outputListNames));
    if ((m_outputListName.empty() && m_outputListNames.empty()) || (!m_outputListName.empty() && !m_outputListNames.empty()))
    {
        std::cout << "StreamingAlgorithm::ReadSettings - Error: You must provide either a single output list name OR a list of output list names"
                  << std::endl;
        return STATUS_CODE_NOT_FOUND;
    }
    if (!m_outputListNames.empty() && m_inputListNames.size() != m_outputListNames.size())
    {
        std::cout << "StreamingAlgorithm::ReadSettings - Error: When providing a list of output lists, there should be a one-to-one "
                  << "correspondence with the list of input lists" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    for (std::string listName : m_inputListNames)
    {
        std::string algStreamName{"Algorithms" + listName};
        if (m_streamAlgorithmMap.find(algStreamName) != m_streamAlgorithmMap.end())
        {
            std::cout << "StreamingAlgorithm::ReadSettings - Error: Duplicate stream name found" << std::endl;
            return STATUS_CODE_INVALID_PARAMETER;
        }
        PANDORA_RETURN_RESULT_IF(
            STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle, algStreamName, m_streamAlgorithmMap[algStreamName]));
        if (m_streamAlgorithmMap.at(algStreamName).empty())
        {
            std::cout << "StreamingAlgorithm::ReadSettings - Error: Found no algorithms for \'" << algStreamName << "\'" << std::endl;
            return STATUS_CODE_NOT_FOUND;
        }
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
