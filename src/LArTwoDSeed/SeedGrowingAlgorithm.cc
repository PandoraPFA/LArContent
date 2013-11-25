/**
 *  @file   LArContent/src/LArTwoDSeed/SeedGrowingAlgorithm.cc
 * 
 *  @brief  Implementation of the seed growing algorithm base class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArParticleIdHelper.h"
#include "LArHelpers/LArVertexHelper.h"

#include "LArTwoDSeed/SeedGrowingAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode SeedGrowingAlgorithm::Run()
{
    const ClusterList *pSeedClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetClusterList(*this, m_seedClusterListName, pSeedClusterList));

    ClusterVector particleSeedVector(pSeedClusterList->begin(), pSeedClusterList->end());
    std::sort(particleSeedVector.begin(), particleSeedVector.end(), LArClusterHelper::SortByNOccupiedLayers);

    const ClusterList *pNonSeedClusterList = NULL;
    const StatusCode statusCode(PandoraContentApi::GetClusterList(*this, m_nonSeedClusterListName, pNonSeedClusterList));

    if ((STATUS_CODE_SUCCESS != statusCode) && (STATUS_CODE_NOT_INITIALIZED != statusCode))
        return statusCode;

    if (STATUS_CODE_NOT_INITIALIZED == statusCode)
        return STATUS_CODE_SUCCESS;

    ClusterUsageMap forwardUsageMap, backwardUsageMap;

    for (ClusterVector::const_iterator iter = particleSeedVector.begin(), iterEnd = particleSeedVector.end(); iter != iterEnd; ++iter)
    {
        ClusterVector candidateClusters;
        this->GetCandidateClusters(pNonSeedClusterList, candidateClusters);
        std::sort(candidateClusters.begin(), candidateClusters.end(), LArClusterHelper::SortByNOccupiedLayers);
        this->FindAssociatedClusters(*iter, candidateClusters, forwardUsageMap, backwardUsageMap);
    }

    SeedAssociationList seedAssociationList;
    this->IdentifyClusterMerges(particleSeedVector, backwardUsageMap, seedAssociationList);
    this->MakeClusterMerges(seedAssociationList);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SeedGrowingAlgorithm::FindAssociatedClusters(Cluster *const pParticleSeed, ClusterVector &candidateClusters,
    ClusterUsageMap &forwardUsageMap, ClusterUsageMap &backwardUsageMap) const
{
    ClusterList currentSeedAssociations, newSeedAssociations;
    currentSeedAssociations.insert(pParticleSeed);

    unsigned int associationOrder(1);

    while (!currentSeedAssociations.empty())
    {
        for (ClusterVector::iterator iterI = candidateClusters.begin(), iterIEnd = candidateClusters.end(); iterI != iterIEnd; ++iterI)
        {
            Cluster *pCandidateCluster = *iterI;

            if (NULL == pCandidateCluster)
                continue;

            for (ClusterList::iterator iterJ = currentSeedAssociations.begin(), iterJEnd = currentSeedAssociations.end(); iterJ != iterJEnd; ++iterJ)
            {
                Cluster *pAssociatedCluster = *iterJ;

                const AssociationType associationType(this->AreClustersAssociated(pAssociatedCluster, pCandidateCluster));

                if (NONE == associationType)
                    continue;

                // Check we store best association between this seed and candidate
                Association association(associationOrder, associationType);
                const Association &existingAssociation = forwardUsageMap[pParticleSeed][pCandidateCluster];

                if (association.GetType() > existingAssociation.GetType())
                {
                    // If not first association, check strength of previous association in chain
                    if (pParticleSeed != pAssociatedCluster)
                        association.SetType(std::min(association.GetType(), backwardUsageMap[pAssociatedCluster][pParticleSeed].GetType()));

                    forwardUsageMap[pParticleSeed][pCandidateCluster] = association;
                    backwardUsageMap[pCandidateCluster][pParticleSeed] = association;
                }

                (void) newSeedAssociations.insert(pCandidateCluster);
                *iterI = NULL;
            }
        }

        currentSeedAssociations = newSeedAssociations;
        newSeedAssociations.clear();
        ++associationOrder;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SeedGrowingAlgorithm::IdentifyClusterMerges(const ClusterVector &particleSeedVector, const ClusterUsageMap &backwardUsageMap,
    SeedAssociationList &seedAssociationList) const
{
    for (ClusterUsageMap::const_iterator iterB = backwardUsageMap.begin(), iterBEnd = backwardUsageMap.end(); iterB != iterBEnd; ++iterB)
    {
        Cluster *pCluster = iterB->first;
        const ClusterAssociationMap &particleSeedUsageMap = iterB->second;

        if (particleSeedUsageMap.empty())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        Cluster *pBestParticleSeed = NULL;
        AssociationType bestType(NONE);
        unsigned int bestOrder(std::numeric_limits<unsigned int>::max());

        for (ClusterAssociationMap::const_iterator iterP = particleSeedUsageMap.begin(), iterPEnd = particleSeedUsageMap.end(); iterP != iterPEnd; ++iterP)
        {
            Cluster *const pParticleSeed = iterP->first;
            const Association &association(iterP->second);

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
        Cluster *const pParticleSeed = *iter;

        if (seedAssociationList.end() == seedAssociationList.find(pParticleSeed))
            seedAssociationList[pParticleSeed] = ClusterVector();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SeedGrowingAlgorithm::MakeClusterMerges(const SeedAssociationList &seedAssociationList) const
{
    for (SeedAssociationList::const_iterator iter = seedAssociationList.begin(), iterEnd = seedAssociationList.end(); iter != iterEnd; ++iter)
    {
        Cluster *pParticleSeed = iter->first;

        for (ClusterVector::const_iterator iterM = iter->second.begin(), iterMEnd = iter->second.end(); iterM != iterMEnd; ++iterM)
        {
            Cluster *pDaughterCluster = *iterM;

            if (pDaughterCluster == pParticleSeed)
                continue;

            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pParticleSeed, pDaughterCluster,
                m_seedClusterListName, m_nonSeedClusterListName));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SeedGrowingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "SeedClusterListName", m_seedClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NonSeedClusterListName", m_nonSeedClusterListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
