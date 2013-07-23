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
    ClusterList allClusterRelegations;

    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

std::cout << " Select cluster for reclustering? " << std::endl;
const unsigned int nHits(pCluster->GetNCaloHits());
const unsigned int nLayers(pCluster->GetOrderedCaloHitList().size());
const float length((pCluster->GetCentroid(pCluster->GetOuterPseudoLayer()) - pCluster->GetCentroid(pCluster->GetInnerPseudoLayer())).GetMagnitude());
const float slidingFitWidth(LArClusterHelper::LArTrackWidth(pCluster));
std::cout << " nHits " << nHits << " nLayers " << nLayers << " length " << length << std::endl;
std::cout << " nHits/Layer " << static_cast<float>(nHits)/static_cast<float>(nLayers) << " nHits/Length " << static_cast<float>(nHits)/length << std::endl;
std::cout << " slidingFitWidth " << slidingFitWidth << std::endl;
ClusterList tempList; tempList.insert(pCluster);
PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
PandoraMonitoringApi::VisualizeClusters(&tempList, "ClusterListU", BLUE);
PandoraMonitoringApi::ViewEvent();

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
            ClusterVector subClusterVector;
            ClusterList relegations;

            for (ClusterList::const_iterator subIter = pReclusterList->begin(), subIterEnd = pReclusterList->end(); subIter != subIterEnd; ++subIter)
            {
                if ((*subIter)->GetNCaloHits() > 9)
                {
                    subClusterVector.push_back(*subIter);
                }
                else
                {
                    relegations.insert(*subIter);
                }
            }

            std::sort(subClusterVector.begin(), subClusterVector.end(), LArClusterHelper::SortByNHits);
            ClusterList parents, branches, prongs;

            for (ClusterVector::iterator subIter = subClusterVector.begin(), subIterEnd = subClusterVector.end(); subIter != subIterEnd; ++subIter)
            {
                Cluster *pParent = *subIter;

                if (NULL == pParent)
                    continue;

                for (ClusterVector::iterator subIter2 = subIter; subIter2 != subIterEnd; ++subIter2)
                {
                    Cluster *pDaughter = *subIter2;

                    if (NULL == pDaughter)
                        continue;

                    const float pOuter(LArClusterHelper::GetClosestDistance(pParent->GetCentroid(pParent->GetOuterPseudoLayer()), pDaughter));
                    const float dOuter(LArClusterHelper::GetClosestDistance(pDaughter->GetCentroid(pDaughter->GetOuterPseudoLayer()), pParent));
                    const float pInner(LArClusterHelper::GetClosestDistance(pParent->GetCentroid(pParent->GetInnerPseudoLayer()), pDaughter));
                    const float dInner(LArClusterHelper::GetClosestDistance(pDaughter->GetCentroid(pDaughter->GetInnerPseudoLayer()), pParent));

                    // Branches
                    if ((pInner > 10.) && (pOuter > 10.) && ((dInner < 2.5) || (dOuter < 2.5)))
                    {
                        parents.insert(pParent);
                        branches.insert(pDaughter);
                        *subIter2 = NULL;
                    }

                    // Prongs
                    if ((pInner > 10.) && (dOuter > 10.) && (pOuter < 2.5) && (dInner < 2.5))
                    {
                        parents.insert(pParent);
                        prongs.insert(pDaughter);
                        *subIter2 = NULL;
                    }
                }
            }

std::cout << " nProtoClusters " << nProtoClusters << ", nProtoClusters10 " << subClusterVector.size() << std::endl;
std::cout << " nParents " << parents.size() << " nBranches " << branches.size() << " nProngs " << prongs.size() << std::endl;
PandoraMonitoringApi::VisualizeClusters(&inputClusterList, "InputClusterList", BLUE);
PandoraMonitoringApi::VisualizeClusters(&parents, "parents", GREEN);
PandoraMonitoringApi::VisualizeClusters(&branches, "branches", ORANGE);
PandoraMonitoringApi::VisualizeClusters(&prongs, "prongs", MAGENTA);
PandoraMonitoringApi::VisualizeClusters(&relegations, "junk", YELLOW);
PandoraMonitoringApi::VisualizeClusters(pReclusterList, "ReclusterList", RED);
PandoraMonitoringApi::ViewEvent();

            if (false)
            {
                chosenListName = reclusterListName;
                //allClusterRelegations.insert();
            }
        }

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndReclustering(*this, chosenListName));
    }

    // Relegate fragments to non-seed list
    if (!allClusterRelegations.empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveClusterList(*this, m_seedClusterListName,
            m_nonSeedClusterListName, allClusterRelegations));
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
