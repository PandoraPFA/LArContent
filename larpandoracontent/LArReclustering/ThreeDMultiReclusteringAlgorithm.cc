/**
 *  @file   larpandoracontent/LArReclustering/ThreeDMultiReclusteringAlgorithm.cc
 *
 *  @brief  Implementation file for the reclustering algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArReclustering/ThreeDMultiReclusteringAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ThreeDMultiReclusteringAlgorithm::ThreeDMultiReclusteringAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDMultiReclusteringAlgorithm::~ThreeDMultiReclusteringAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDMultiReclusteringAlgorithm::Run()
{
    std::cout << "\n";
    const PfoList *pPfos {nullptr};
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=,
        PandoraContentApi::GetList(*this, m_pfoListName, pPfos));
    if (!pPfos || pPfos->empty())
    {
        std::cout << "pPfos empty or not initialised\n";
        return STATUS_CODE_SUCCESS;
    }
    std::cout << "pPfos->size()=" << pPfos->size() << "\n";

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_clusterListName));

    PfoList pfosToRecluster;
    for (const Pfo *const pPfo : *pPfos)
    {
        // To be replaced by an appropriate figure of merit calculation with:
        // - Configurable threshold
        // - Probably with implemented in an alg tool
        if (rand() % 2 == 1)
            pfosToRecluster.push_back(pPfo);
    }
    std::cout << "pfosToRecluster.size()=" << pfosToRecluster.size() << "\n";
    if (pfosToRecluster.empty())
        return STATUS_CODE_SUCCESS;

    ClusterList freedClusters3D;
    for (const Pfo *const pPfo : pfosToRecluster)
    {
        ClusterList clusters3D;
        LArPfoHelper::GetThreeDClusterList(pPfo, clusters3D);
        if (clusters3D.empty())
        {
            std::cout << "clusters3D.empty()\n";
            continue;
        }
        if (clusters3D.size() != 1)
        {
            std::cout << "clusters3D.size() != 1, clusters.size()=" << clusters3D.size() << "\n";
            continue;
        }
        for (const Cluster *const pCluster3D : clusters3D)
        {
            freedClusters3D.push_back(pCluster3D);
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                PandoraContentApi::RemoveFromPfo(*this, pPfo, pCluster3D));
        }
    }
    std::cout << "freedClusters3D.size()=" << freedClusters3D.size() << "\n";

    std::string originalClusterListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        PandoraContentApi::InitializeReclustering(*this, TrackList(), freedClusters3D, originalClusterListName));
    // std::string name;
    // PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
    //     PandoraContentApi::GetCurrentListName<Cluster>(*this, name));
    // const ClusterList *pCurrClusters;
    // PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
    //     PandoraContentApi::GetCurrentList(*this, pCurrClusters));
    // std::cout << "InitializeReclustering: current list is " << name << " with " << pCurrClusters->size() << " clusters.\n";

    for (const std::string &clusteringAlg : m_clusteringAlgs)
    {
        std::string reclusterListName;
        const ClusterList *pReclusterList {nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
            PandoraContentApi::RunClusteringAlgorithm(*this, clusteringAlg, pReclusterList, reclusterListName));
        if (!pReclusterList)
        {
            std::cout << "!pReclusterList\n";
            continue;
        }
        if (pReclusterList->empty())
        {
            std::cout << "pReclusterList->empty()\n";
            continue;
        }
        std::cout << "pReclusterList->size()=" << pReclusterList->size() << "\n";
    }

    // 1. InitialiseReclustering
    // 2. Loop reclustering algs
    // 3. Run reclutering alg
    // 4. Compute reclustering metric, track the best one
    // 5. Finish loop, choose the best new list (if the better than the original clustering metric)
    // 6. EndReclustering passing the chosen best list to it
    // 7. A load of mess to make the new Pfos from the new 3d clusters, ensuring that the original 2D clusters go to the correct place

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDMultiReclusteringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListName", m_clusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle, "ClusteringAlgorithms", m_clusteringAlgs));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
