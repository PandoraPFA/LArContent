/**
 *  @file   larpandoracontent/LArTrackShowerId/BranchGrowingAlgorithm.cc
 *
 *  @brief  Implementation of the branch growing algorithm base class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArTrackShowerId/BranchGrowingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

void BranchGrowingAlgorithm::FindAssociatedClusters(const Cluster *const pParticleSeed, ClusterVector &candidateClusters,
    ClusterUsageMap &forwardUsageMap, ClusterUsageMap &backwardUsageMap) const
{
    ClusterVector currentSeedAssociations, newSeedAssociations;
    currentSeedAssociations.push_back(pParticleSeed);

    unsigned int associationOrder(1);
    std::map<std::pair<const Cluster *, const Cluster *>, Association> associationCache;

    while (!currentSeedAssociations.empty())
    {
        for (ClusterVector::iterator iterI = candidateClusters.begin(), iterIEnd = candidateClusters.end(); iterI != iterIEnd; ++iterI)
        {
            const Cluster *const pCandidateCluster = *iterI;

            if (NULL == pCandidateCluster)
                continue;

            for (ClusterVector::iterator iterJ = currentSeedAssociations.begin(), iterJEnd = currentSeedAssociations.end(); iterJ != iterJEnd; ++iterJ)
            {
                const Cluster *const pAssociatedCluster = *iterJ;
                const auto clusterPair = pAssociatedCluster->GetNCaloHits() > pCandidateCluster->GetNCaloHits()
                    ? std::make_pair(pCandidateCluster, pAssociatedCluster)
                    : std::make_pair(pAssociatedCluster, pCandidateCluster);
                Association association;

                if (associationCache.count(clusterPair) > 0)
                    association = associationCache[clusterPair];
                else
                {
                    const AssociationType associationType(this->AreClustersAssociated(pAssociatedCluster, pCandidateCluster));
                    association = Association(associationOrder, associationType);
                    associationCache[clusterPair] = association;
                }

                if (association.GetType() == NONE)
                    continue;

                // Check we store best association between this seed and candidate
                const Association &existingAssociation = forwardUsageMap[pParticleSeed][pCandidateCluster];

                if (association.GetType() > existingAssociation.GetType())
                {
                    // If not first association, check strength of previous association in chain
                    if (pParticleSeed != pAssociatedCluster)
                        association.SetType(std::min(association.GetType(), backwardUsageMap[pAssociatedCluster][pParticleSeed].GetType()));

                    forwardUsageMap[pParticleSeed][pCandidateCluster] = association;
                    backwardUsageMap[pCandidateCluster][pParticleSeed] = association;
                }

                newSeedAssociations.push_back(pCandidateCluster);
                *iterI = NULL;
            }
        }

        currentSeedAssociations = newSeedAssociations;
        newSeedAssociations.clear();
        ++associationOrder;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BranchGrowingAlgorithm::IdentifyClusterMerges(
    const ClusterVector &particleSeedVector, const ClusterUsageMap &backwardUsageMap, SeedAssociationList &seedAssociationList) const
{
    ClusterVector sortedCandidates;
    for (const auto &mapEntry : backwardUsageMap)
        sortedCandidates.push_back(mapEntry.first);
    std::sort(sortedCandidates.begin(), sortedCandidates.end(), LArClusterHelper::SortByNHits);

    for (const Cluster *const pCluster : sortedCandidates)
    {
        const ClusterAssociationMap &particleSeedUsageMap(backwardUsageMap.at(pCluster));

        if (particleSeedUsageMap.empty())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        ClusterVector sortedSeeds;
        for (const auto &mapEntry : particleSeedUsageMap)
            sortedSeeds.push_back(mapEntry.first);
        std::sort(sortedSeeds.begin(), sortedSeeds.end(), LArClusterHelper::SortByNHits);

        const Cluster *pBestParticleSeed = NULL;
        AssociationType bestType(NONE);
        unsigned int bestOrder(std::numeric_limits<unsigned int>::max());

        for (const Cluster *const pParticleSeed : sortedSeeds)
        {
            const Association &association(particleSeedUsageMap.at(pParticleSeed));

            if ((association.GetType() > bestType) || ((association.GetType() == bestType) && (association.GetOrder() < bestOrder)))
            {
                // Break-out condition for single order associations
                if ((SINGLE_ORDER == association.GetType()) && (association.GetOrder() > 1))
                    continue;

                // Type is primary consideration; order breaks ties
                pBestParticleSeed = pParticleSeed;
                bestType = association.GetType();
                bestOrder = association.GetOrder();
            }
            else if ((association.GetType() == bestType) && (association.GetOrder() == bestOrder))
            {
                // Remove ambiguous cluster from algorithm
                pBestParticleSeed = NULL;
            }
        }

        if (NULL == pBestParticleSeed)
            continue;

        seedAssociationList[pBestParticleSeed].push_back(pCluster);
    }

    // Now deal with seeds that have no associations
    for (ClusterVector::const_iterator iter = particleSeedVector.begin(), iterEnd = particleSeedVector.end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pParticleSeed = *iter;

        if (seedAssociationList.end() == seedAssociationList.find(pParticleSeed))
            seedAssociationList[pParticleSeed] = ClusterVector();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode BranchGrowingAlgorithm::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
