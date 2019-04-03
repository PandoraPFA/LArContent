/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterMopUp/ExampleClusterMopUpAlgorithm.cc
 * 
 *  @brief  Implementation of the example cluster mop up algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterMopUp/ExampleClusterMopUpAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ExampleClusterMopUpAlgorithm::ExampleClusterMopUpAlgorithm() :
    m_maxMergeDistance(5.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ExampleClusterMopUpAlgorithm::Run()
{
    const ClusterList *pClusterList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

    // Need to be very careful with cluster list iterators here, as we are deleting elements from the std::list owned by the manager.
    // If user chooses to iterate over that same list, must adhere to rule that iterators pointing at the deleted element will be invalidated.

    // Here, iterate over a local copy of the cluster list, and keep track of dangling pointers
    ClusterSet deletedClusters;
    const ClusterList localClusterList(*pClusterList);

    for (const Cluster *const pParentCluster : localClusterList)
    {
        // Check to see whether parent cluster (address stored in local copy of list) still exists in manager-owned list, and hasn't been
        // removed by the cluster merging operations in this algorithm. Many alternative methods to check this, of course.
        if (deletedClusters.count(pParentCluster))
            continue;

        std::cout << "ExampleClusterMopUpAlgorithm: new candidate parent cluster " << std::endl;

        const Cluster *pBestDaughterCluster(nullptr);
        float bestDistance(std::numeric_limits<float>::max());

        for (const Cluster *const pDaughterCluster : *pClusterList)
        {
            if (pParentCluster == pDaughterCluster)
                continue;

            // TODO Calculate association details, using pParentCluster and pDaughterCluster. Just a simple example here.
            const float distance(LArClusterHelper::GetClosestDistance(pParentCluster, pDaughterCluster));

            if ((distance < m_maxMergeDistance) && (distance < bestDistance))
            {
                bestDistance = distance;
                pBestDaughterCluster = pDaughterCluster;
            }
        }

        if (pBestDaughterCluster)
        {
            std::cout << "Make cluster merge, with distance " << bestDistance << std::endl;
            const ClusterList parentClusterList(1, pParentCluster);
            const ClusterList bestDaughterClusterList(1, pBestDaughterCluster);
            PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f);
            PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &parentClusterList, "parentClusterList", RED);
            PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &bestDaughterClusterList, "bestDaughterClusterList", BLUE);
            PandoraMonitoringApi::ViewEvent(this->GetPandora());

            // The API implementation will enforce the availability of the daughter cluster and ensure that the parent and daughter are not one and the same
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pParentCluster, pBestDaughterCluster));

            const ClusterList afterMergeList(1, pParentCluster);
            PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &afterMergeList, "afterMergeList", GREEN);
            PandoraMonitoringApi::ViewEvent(this->GetPandora());

            // pBestDaughterCluster is now a dangling pointer, which exists only in the local cluster list - do not deference!
            deletedClusters.insert(pBestDaughterCluster);
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ExampleClusterMopUpAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxMergeDistance", m_maxMergeDistance));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
