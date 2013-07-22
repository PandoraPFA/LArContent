/**
 *  @file   LArContent/src/LArReclustering/TrackSplittingAlgorithm.cc
 * 
 *  @brief  Implementation of the track splitting algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArReclustering/TrackSplittingAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode TrackSplittingAlgorithm::Run()
{
    // Reclustering operations require input list to be current list
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentClusterList(*this, m_seedClusterListName));

    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentClusterList(*this, pClusterList));

    ClusterVector clusterVector(pClusterList->begin(), pClusterList->end());
    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNOccupiedLayers);

    ClusterList clusterRelegations;

    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        // Selection here
        //if (true)
        //    continue;

        // Perform the reclustering
        ClusterList inputClusterList;
        inputClusterList.insert(pCluster);
        std::string inputClusterListName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeReclustering(*this, TrackList(), inputClusterList, inputClusterListName));

        const ClusterList *pReclusterList = NULL;
        std::string reclusterListName;

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunClusteringAlgorithm(*this, m_clusteringAlgorithmName, pReclusterList, reclusterListName));

        for (StringVector::const_iterator iter = m_clusterAssociationAlgorithms.begin(), iterEnd = m_clusterAssociationAlgorithms.end(); iter != iterEnd; ++iter)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunDaughterAlgorithm(*this, *iter));

        // Now examine resulting cluster list
        std::string chosenListName(inputClusterListName);
        const unsigned int nProtoClusters(pReclusterList->size());

        if (nProtoClusters > 1)
        {
            ClusterVector subClusterVector(pReclusterList->begin(), pReclusterList->end());
            std::sort(subClusterVector.begin(), subClusterVector.end(), LArClusterHelper::SortByNHits);

            for (ClusterVector::const_iterator subIter = subClusterVector.begin(), subIterEnd = subClusterVector.end(); subIter != subIterEnd; ++subIter)
            {
                Cluster *pSubCluster = *subIter;

                if (pSubCluster->GetNCaloHits() < 10)
                {
                    clusterRelegations.insert(pSubCluster);
                    continue;
                }

                if (LArClusterHelper::LArTrackWidth(pSubCluster) < 1.5f)
                {
                    chosenListName = reclusterListName;
                    ClusterList subClusterList;
                    subClusterList.insert(pSubCluster);
                    PandoraMonitoringApi::VisualizeClusters(&subClusterList, "SubCluster", BLUE);
                    PandoraMonitoringApi::ViewEvent();
                }
            }

            PandoraMonitoringApi::VisualizeClusters(&inputClusterList, "InputClusterList", BLUE);
            PandoraMonitoringApi::VisualizeClusters(pReclusterList, "ReclusterList", RED);
            PandoraMonitoringApi::ViewEvent();
        }

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndReclustering(*this, chosenListName));
    }

    // Relegate fragments to non-seed list
    if (!clusterRelegations.empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveClusterList(*this, m_seedClusterListName,
            m_nonSeedClusterListName, clusterRelegations));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "SeedClusterListName", m_seedClusterListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NonSeedClusterListName", m_nonSeedClusterListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, xmlHandle, "Clustering",
        m_clusteringAlgorithmName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle, "associationAlgorithms",
        m_clusterAssociationAlgorithms));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
