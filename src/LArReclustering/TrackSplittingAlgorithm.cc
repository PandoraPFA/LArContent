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
    const ClusterList *pInputClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentClusterList(*this, m_inputSeedClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentClusterList(*this, pInputClusterList));

    ClusterVector clusterVector(pInputClusterList->begin(), pInputClusterList->end());
    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNOccupiedLayers);
PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
PandoraMonitoringApi::VisualizeClusters(pInputClusterList, "inputClusterList", BLUE);
PandoraMonitoringApi::ViewEvent();

    // The target lists
    ClusterList allShowerSeedClusters, allTrackSeedClusters, allNonSeedClusters;

    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;
        ClusterList trackSeedClusters, nonSeedClusters;

        // Obvious showers
        const unsigned int nHits(pCluster->GetNCaloHits());

        if (0 == nHits)
            return STATUS_CODE_FAILURE;

        const float length((pCluster->GetCentroid(pCluster->GetOuterPseudoLayer()) - pCluster->GetCentroid(pCluster->GetInnerPseudoLayer())).GetMagnitude());
        const float hitsPerUnitLength((length > 0.f) ? static_cast<float>(nHits) / length : std::numeric_limits<float>::max());

        if (hitsPerUnitLength > 10.f) // TODO
        {
            allShowerSeedClusters.insert(pCluster);
std::cout << " hitsPerUnitLength " << hitsPerUnitLength << std::endl;
PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
ClusterList tempList; tempList.insert(pCluster);
PandoraMonitoringApi::VisualizeClusters(&tempList, "ObviousShower", BLUE);
PandoraMonitoringApi::ViewEvent();
            continue;
        }

        // Reduce to subclusters
        ClusterList inputClusterList; inputClusterList.insert(pCluster); std::string inputClusterListName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeReclustering(*this, TrackList(), inputClusterList, inputClusterListName));

        const ClusterList *pReclusterList = NULL; std::string reclusterListName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunClusteringAlgorithm(*this, m_clusteringAlgorithmName, pReclusterList, reclusterListName));

        for (StringVector::const_iterator iter = m_clusterAssociationAlgorithms.begin(), iterEnd = m_clusterAssociationAlgorithms.end(); iter != iterEnd; ++iter)
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunDaughterAlgorithm(*this, *iter));
        }

        // Examine subclusters
        ClusterVector subClusterVector;

        for (ClusterList::const_iterator subIter = pReclusterList->begin(), subIterEnd = pReclusterList->end(); subIter != subIterEnd; ++subIter)
        {
            Cluster *pSubCluster = *subIter;

            if (pSubCluster->GetNCaloHits() >= 10)
            {
                subClusterVector.push_back(pSubCluster);
            }
            else
            {
                nonSeedClusters.insert(pSubCluster);
            }
        }

        std::sort(subClusterVector.begin(), subClusterVector.end(), LArClusterHelper::SortByNHits);

        // Identify clusters which are branches or prongs
        ClusterList branches, prongs, parents;

        for (ClusterVector::iterator subIter = subClusterVector.begin(), subIterEnd = subClusterVector.end(); subIter != subIterEnd; ++subIter)
        {
            Cluster *pParent = *subIter;

            for (ClusterVector::iterator subIter2 = subIter; subIter2 != subIterEnd; ++subIter2)
            {
                Cluster *pDaughter = *subIter2;

                if (pDaughter == pParent)
                    continue;

                const CartesianVector pInnerCentroid(pParent->GetCentroid(pParent->GetInnerPseudoLayer()));
                const CartesianVector pOuterCentroid(pParent->GetCentroid(pParent->GetOuterPseudoLayer()));
                const CartesianVector dInnerCentroid(pDaughter->GetCentroid(pDaughter->GetInnerPseudoLayer()));
                const CartesianVector dOuterCentroid(pDaughter->GetCentroid(pDaughter->GetOuterPseudoLayer()));

                // prongs
                const float rsqInnerInner((pInnerCentroid - dInnerCentroid).GetMagnitudeSquared());
                const float rsqInnerOuter((pInnerCentroid - dOuterCentroid).GetMagnitudeSquared());
                const float rsqOuterInner((pOuterCentroid - dInnerCentroid).GetMagnitudeSquared());
                const float rsqOuterOuter((pOuterCentroid - dOuterCentroid).GetMagnitudeSquared()); 

                if (rsqInnerInner < 2.5 * 2.5 || rsqInnerOuter < 2.5 * 2.5 ||
                    rsqOuterInner < 2.5 * 2.5 || rsqOuterOuter < 2.5 * 2.5) 
                {
                    parents.insert(pParent);
                    prongs.insert(pDaughter);
                    continue;
                }

                // branches
                const float pOuter(LArClusterHelper::GetClosestDistance(pOuterCentroid, pDaughter));
                const float pInner(LArClusterHelper::GetClosestDistance(pInnerCentroid, pDaughter));
                const float dOuter(LArClusterHelper::GetClosestDistance(dOuterCentroid, pParent));
                const float dInner(LArClusterHelper::GetClosestDistance(dInnerCentroid, pParent));

                if ((pInner > 2.5) && (pOuter > 2.5) && ((dInner < 2.5) || (dOuter < 2.5)))
                {
                    parents.insert(pParent);
                    branches.insert(pDaughter);
                    continue;
                }
            }
        }

        // Select clean clusters
        unsigned int nParents(0), nProngs(0), nBranches(0), nIsolated(0);

        for (ClusterVector::iterator subIter = subClusterVector.begin(), subIterEnd = subClusterVector.end(); subIter != subIterEnd; ++subIter)
        {
            Cluster *pSubCluster = *subIter;

            const bool isParent(parents.count(pSubCluster));
            const bool isBranch(branches.count(pSubCluster));
            const bool isProng(prongs.count(pSubCluster));
            bool isIsolated(false);

            if (isBranch)
            {
                ++nBranches;
            }
            else if (isProng)
            {
                ++nProngs;
            }
            else if (isParent)
            {
                ++nParents;
            }
            else
            {
                ++nIsolated;
                isIsolated = true;
            }

            const float lengthSquared(LArClusterHelper::GetLengthSquared(pSubCluster));
            const float hitFraction(static_cast<float>(pSubCluster->GetNCaloHits()) / static_cast<float>(nHits));

            if ((lengthSquared > 200.f) || (!isBranch && (lengthSquared > 100.f)) || (isIsolated && (hitFraction > 0.8f)))
            {
                trackSeedClusters.insert(pSubCluster);
            }
            else
            {
                nonSeedClusters.insert(pSubCluster);
            }
        }

std::cout << " nProtoClusters " << pReclusterList->size() << ", nProtoClusters10 " << subClusterVector.size() << std::endl;
std::cout << " nParents " << nParents << " nBranches " << nBranches << " nProngs " << nProngs << " nIsolated " << nIsolated << std::endl;
std::cout << " seedClusters " << trackSeedClusters.size() << " nonSeedClusters " << nonSeedClusters.size() << std::endl;

const unsigned int nProtoClusters10(subClusterVector.size());
const float clustersPerUnitLength(static_cast<float>(subClusterVector.size()) / length);
const float branchesPerUnitLength(static_cast<float>(nBranches) / length);
unsigned int nSelectedHits(0);
for (ClusterList::const_iterator tIter = trackSeedClusters.begin(), tIterEnd = trackSeedClusters.end(); tIter != tIterEnd; ++tIter)
    nSelectedHits += (*tIter)->GetNCaloHits();
const float selectedHitFraction(static_cast<float>(nSelectedHits) / static_cast<float>(nHits));
std::cout << " clustersPerUnitLength " << clustersPerUnitLength << " branchesPerUnitLength " << branchesPerUnitLength << " selectedHitFraction " << selectedHitFraction << std::endl;

        std::string chosenListName;

        if ((selectedHitFraction > 0.8f) && ((nProtoClusters10 == 1) || (clustersPerUnitLength < 0.1f)))
        {
std::cout << "TRACK-LIKE - Use new fragments " << std::endl;
            chosenListName = reclusterListName;
            allTrackSeedClusters.insert(trackSeedClusters.begin(), trackSeedClusters.end());
            allNonSeedClusters.insert(nonSeedClusters.begin(), nonSeedClusters.end());
        }
        else
        {
std::cout << "SHOWER - Revert to input cluster " << std::endl;
            chosenListName = inputClusterListName;
            allShowerSeedClusters.insert(pCluster);
        }
//PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
//PandoraMonitoringApi::VisualizeClusters(&inputClusterList, "inputClusterList", BLUE);
//PandoraMonitoringApi::VisualizeClusters(&trackSeedClusters, "trackSeedClusters", RED);
//PandoraMonitoringApi::VisualizeClusters(&nonSeedClusters, "nonSeedClusters", GREEN);
//PandoraMonitoringApi::VisualizeClusters(&branches, "branches", ORANGE);
//PandoraMonitoringApi::VisualizeClusters(&prongs, "prongs", MAGENTA);
//PandoraMonitoringApi::VisualizeClusters(pReclusterList, "reclusterList", GRAY);
//PandoraMonitoringApi::ViewEvent();

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndReclustering(*this, chosenListName));
    }

    // Cluster output
    if (!allShowerSeedClusters.empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveClusterList(*this, m_inputSeedClusterListName,
            m_showerSeedClusterListName, allShowerSeedClusters));
    }

    if (!allTrackSeedClusters.empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveClusterList(*this, m_inputSeedClusterListName,
            m_trackSeedClusterListName, allTrackSeedClusters));
    }

    if (!allNonSeedClusters.empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveClusterList(*this, m_inputSeedClusterListName,
            m_nonSeedClusterListName, allNonSeedClusters));
    }

    if (!pInputClusterList->empty())
        return STATUS_CODE_FAILURE;
PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
PandoraMonitoringApi::VisualizeClusters(pInputClusterList, "InputClusterList", BLUE);
PandoraMonitoringApi::VisualizeClusters(&allShowerSeedClusters, "allShowerSeedClusters", RED);
PandoraMonitoringApi::VisualizeClusters(&allTrackSeedClusters, "allTrackSeedClusters", GREEN);
PandoraMonitoringApi::VisualizeClusters(&allNonSeedClusters, "allNonSeedClusters", ORANGE);
PandoraMonitoringApi::ViewEvent();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputSeedClusterListName",
        m_inputSeedClusterListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ShowerSeedClusterListName",
        m_showerSeedClusterListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrackSeedClusterListName",
        m_trackSeedClusterListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NonSeedClusterListName",
        m_nonSeedClusterListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, xmlHandle, "Clustering",
        m_clusteringAlgorithmName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle, "associationAlgorithms",
        m_clusterAssociationAlgorithms));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
