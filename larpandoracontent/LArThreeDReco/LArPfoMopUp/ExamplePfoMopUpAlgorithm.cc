/**
 *  @file   larpandoracontent/LArThreeDReco/LArPfoMopUp/ExamplePfoMopUpAlgorithm.cc
 * 
 *  @brief  Implementation of the example pfo mop up algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArThreeDReco/LArPfoMopUp/ExamplePfoMopUpAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ExamplePfoMopUpAlgorithm::ExamplePfoMopUpAlgorithm() :
    m_maxMergeDistance(5.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ExamplePfoMopUpAlgorithm::Run()
{
    ClusterList clusterList3D;
    ClusterToPfoMap clusterToPfoMap;
    this->GetThreeDClusters(clusterList3D, clusterToPfoMap);

    ClusterMergeMap clusterMergeMap;
    this->GetClusterMergeMap(clusterList3D, clusterMergeMap);

    this->MakePfoMerges(clusterList3D, clusterToPfoMap, clusterMergeMap);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ExamplePfoMopUpAlgorithm::GetThreeDClusters(ClusterList &clusterList3D, ClusterToPfoMap &clusterToPfoMap) const
{
    for (const std::string &pfoListName : m_inputPfoListNames)
    {
        const PfoList *pPfoList(nullptr);

        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, pfoListName, pPfoList))
        {
            std::cout << "ExamplePfoMopUpAlgorithm: could not find pfo list with name " << pfoListName << std::endl;
            continue;
        }

        for (const Pfo *const pPfo : *pPfoList)
        {
            ClusterList pfoClusters3D;
            LArPfoHelper::GetThreeDClusterList(pPfo, pfoClusters3D);

            for (const Cluster *const pCluster3D : pfoClusters3D)
            {
                if (!clusterToPfoMap.insert(ClusterToPfoMap::value_type(pCluster3D, pPfo)).second)
                    throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);

                clusterList3D.push_back(pCluster3D);
            }
        }
    }

    clusterList3D.sort(LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ExamplePfoMopUpAlgorithm::GetClusterMergeMap(const ClusterList &clusters3D, ClusterMergeMap &clusterMergeMap) const
{
    unsigned int displayNumber(0);

    for (const Cluster *const pParentCluster : clusters3D)
    {
        std::cout << "ExamplePfoMopUpAlgorithm: new candidate parent cluster " << std::endl;

        for (const Cluster *const pDaughterCluster : clusters3D)
        {
            if (pDaughterCluster == pParentCluster)
                continue;

            // TODO Calculate association details, using pParentCluster and pDaughterCluster. Just a simple example here.
            const float distance(LArClusterHelper::GetClosestDistance(pParentCluster, pDaughterCluster));

            std::cout << "ExamplePfoMopUpAlgorithm: displayNumber " << displayNumber++ << ", distance " << distance << std::endl;
            const ClusterList parentClusterList(1, pParentCluster);
            const ClusterList daughterClusterList(1, pDaughterCluster);
            PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_DEFAULT, -1.f, -1.f, 1.f);
            PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &parentClusterList, "parentClusterList", RED);
            PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &daughterClusterList, "daughterClusterList", BLUE);
            PandoraMonitoringApi::ViewEvent(this->GetPandora());

            if (distance < m_maxMergeDistance)
                clusterMergeMap[pParentCluster].emplace_back(pDaughterCluster, distance);
        }
    }

    // Rank merges by figure of merit
    for (ClusterMergeMap::value_type &mapEntry : clusterMergeMap)
        std::sort(mapEntry.second.begin(), mapEntry.second.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ExamplePfoMopUpAlgorithm::MakePfoMerges(const ClusterList &clusters3D, const ClusterToPfoMap &clusterToPfoMap, const ClusterMergeMap &clusterMergeMap) const
{
    // ATTN All the complex steps here are to address cases where Pfos are merged, then merged again, etc.
    ClusterReplacementMap clusterReplacementMap;

    for (const Cluster *pParentCluster : clusters3D)
    {
        if (!clusterMergeMap.count(pParentCluster))
            continue;

        // Have these clusters already been affected by prior merges in this algorithm?
        const Cluster *pBestDaughterCluster(clusterMergeMap.at(pParentCluster).at(0).GetDaughterCluster());
        pParentCluster = (clusterReplacementMap.count(pParentCluster)) ? clusterReplacementMap.at(pParentCluster) : pParentCluster;
        pBestDaughterCluster = (clusterReplacementMap.count(pBestDaughterCluster)) ? clusterReplacementMap.at(pBestDaughterCluster) : pBestDaughterCluster;

        if (pParentCluster == pBestDaughterCluster)
            continue;

        std::cout << "Make pfo merge, with distance " << clusterMergeMap.at(pParentCluster).at(0).GetDistance() << std::endl;
        const PfoList parentPfoList(1, clusterToPfoMap.at(pParentCluster));
        const PfoList daughterPfoList(1, clusterToPfoMap.at(pBestDaughterCluster));
        PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &parentPfoList, "parentPfoList", RED);
        PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &daughterPfoList, "daughterPfoList", BLUE);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());

        // Make the actual merge
        this->MergeAndDeletePfos(clusterToPfoMap.at(pParentCluster), clusterToPfoMap.at(pBestDaughterCluster));

        const PfoList afterMergeList(1, clusterToPfoMap.at(pParentCluster));
        PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &afterMergeList, "afterMergeList", GREEN);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());

        // Book-keeping for reciprocal relationships and progressive merges
        clusterReplacementMap[pBestDaughterCluster] = pParentCluster;

        for (ClusterReplacementMap::value_type &mapEntry : clusterReplacementMap)
        {
            if (pBestDaughterCluster == mapEntry.second)
                mapEntry.second = pParentCluster;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ExamplePfoMopUpAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "InputPfoListNames", m_inputPfoListNames));

    m_daughterListNames.insert(m_daughterListNames.end(), m_inputPfoListNames.begin(), m_inputPfoListNames.end());

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxMergeDistance", m_maxMergeDistance));

    return PfoMopUpBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
