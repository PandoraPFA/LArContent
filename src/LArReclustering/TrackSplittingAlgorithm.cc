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
	    // Move small sub-clusters into relegated list
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

            // Identify clusters which are branches or prongs
            ClusterList branches, prongs, parents, daughters;  // TODO: Which of these lists is superfluous ?

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
                        daughters.insert(pDaughter);
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
                        daughters.insert(pDaughter);
                        branches.insert(pDaughter);
			continue;
                    }
                }
            }

            // Select clean clusters
            ClusterList theSeedClusters, theNonSeedClusters;

            unsigned int nParents(0), nProngs(0), nBranches(0), nIsolated(0);

            for (ClusterVector::iterator subIter = subClusterVector.begin(), subIterEnd = subClusterVector.end(); subIter != subIterEnd; ++subIter)
            {
                Cluster* pCluster = *subIter;

                bool isParent(false), isDaughter(false), isBranch(false), isProng(false);

                ClusterList::const_iterator parentIter = parents.find(pCluster);
                if (parentIter != parents.end()) isParent = true;

                ClusterList::const_iterator daughterIter = daughters.find(pCluster);
                if (daughterIter != daughters.end()) isDaughter = true;

                ClusterList::const_iterator branchIter = branches.find(pCluster);
                if (branchIter != branches.end()) isBranch = true;

                ClusterList::const_iterator prongIter = prongs.find(pCluster);
                if (prongIter != prongs.end()) isProng = true;
		
                if (true == isBranch) ++nBranches;
                else if (true == isProng) ++nProngs;
                else if (false == isParent) ++nIsolated;
                else ++nParents;

		const float lengthSquared(LArClusterHelper::GetLengthSquared(pCluster));

		if ( lengthSquared > 200.f || (false == isBranch && lengthSquared > 100.f) )
		{
		    theSeedClusters.insert(pCluster);
		}
                else 
                {
                    theNonSeedClusters.insert(pCluster);
		}
	    }

std::cout << " nProtoClusters " << nProtoClusters << ", nProtoClusters10 " << subClusterVector.size() << std::endl;
std::cout << " nParents " << nParents << " nBranches " << nBranches << " nProngs " << nProngs << " nIsolated " << nIsolated << std::endl;
std::cout << " seedClusters " << theSeedClusters.size() << " nonSeedClusters " << theNonSeedClusters.size() << std::endl;
PandoraMonitoringApi::VisualizeClusters(&inputClusterList, "InputClusterList", BLUE);
PandoraMonitoringApi::VisualizeClusters(&theSeedClusters, "Seeds", RED);
PandoraMonitoringApi::VisualizeClusters(&theNonSeedClusters, "NonSeeds", GREEN);
PandoraMonitoringApi::VisualizeClusters(&branches, "branches", ORANGE);
PandoraMonitoringApi::VisualizeClusters(&prongs, "prongs", MAGENTA);
PandoraMonitoringApi::VisualizeClusters(&relegations, "junk", YELLOW);
PandoraMonitoringApi::VisualizeClusters(pReclusterList, "ReclusterList", GRAY);
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
