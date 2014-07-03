/**
 *  @file   LArContent/src/LArTwoDReco/LArSeedFinding/ClusterCharacterisationAlgorithm.cc
 * 
 *  @brief  Implementation of the cluster characterisation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArTwoDReco/LArSeedFinding/ClusterCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode ClusterCharacterisationAlgorithm::Run()
{
    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputClusterListName, pClusterList));

    ClusterList usedClusters;
    Cluster *pSeedCluster = NULL;

    while (this->GetNextSeedCandidate(pClusterList, usedClusters, pSeedCluster))
    {
        SeedAssociationList seedAssociationList;
        this->GetSeedAssociationList(ClusterVector(1, pSeedCluster), pClusterList, seedAssociationList);

        if (seedAssociationList.size() != 1)
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        SeedAssociationList finalSeedAssociationList;
        this->CheckSeedAssociationList(seedAssociationList.begin(), finalSeedAssociationList);

        for (SeedAssociationList::const_iterator iter = finalSeedAssociationList.begin(), iterEnd = finalSeedAssociationList.end(); iter != iterEnd; ++iter)
        {
            usedClusters.insert(iter->first);
            usedClusters.insert(iter->second.begin(), iter->second.end());

            for (ClusterVector::const_iterator iter2 = iter->second.begin(), iter2End = iter->second.end(); iter2 != iter2End; ++iter2)
            {
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, iter->first, *iter2, m_inputClusterListName, m_inputClusterListName));
            }
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ClusterCharacterisationAlgorithm::GetNextSeedCandidate(const ClusterList *const pClusterList, const ClusterList &usedClusters,
    Cluster *&pSeedCluster) const
{
    pSeedCluster = NULL;

    ClusterVector clusterVector;
    clusterVector.insert(clusterVector.end(), pClusterList->begin(), pClusterList->end());
    std::sort(clusterVector.begin(), clusterVector.end(), ClusterCharacterisationAlgorithm::SortClusters);

    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (usedClusters.count(pCluster))
            continue;

        if (pCluster->IsAvailable() && (pCluster->GetNCaloHits() < 5)) // TODO
            continue;

        pSeedCluster = pCluster;
        return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterCharacterisationAlgorithm::GetSeedAssociationList(const ClusterVector &particleSeedVector, const ClusterList *const pClusterList,
    SeedAssociationList &seedAssociationList) const
{
    ClusterVector candidateClusters;
    ClusterList seedClusters(particleSeedVector.begin(), particleSeedVector.end());

    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster *pCandidateCluster = *iter;

        if (!pCandidateCluster->IsAvailable() || seedClusters.count(pCandidateCluster) || (pCandidateCluster->GetNCaloHits() < 5)) // TODO
            continue;

        candidateClusters.push_back(pCandidateCluster);
    }

    std::sort(candidateClusters.begin(), candidateClusters.end(), ClusterCharacterisationAlgorithm::SortClusters);
    ClusterUsageMap forwardUsageMap, backwardUsageMap;

    for (ClusterVector::const_iterator iter = particleSeedVector.begin(), iterEnd = particleSeedVector.end(); iter != iterEnd; ++iter)
    {
        this->FindAssociatedClusters(*iter, candidateClusters, forwardUsageMap, backwardUsageMap);
    }

    this->IdentifyClusterMerges(particleSeedVector, backwardUsageMap, seedAssociationList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterCharacterisationAlgorithm::CheckSeedAssociationList(SeedAssociationList::const_iterator seedIter, SeedAssociationList &finalSeedAssociationList) const
{
    ClusterList usedClusters, availableClusters;
    usedClusters.insert(seedIter->first);
    availableClusters.insert(seedIter->second.begin(), seedIter->second.end());

    SeedAssociationList seedAssociationList;
    seedAssociationList.insert(SeedAssociationList::value_type(seedIter->first, seedIter->second));
    const float fom1(std::max(1.f, this->GetFigureOfMerit(seedAssociationList)));

    Cluster *pTrialSeedCluster = NULL;
    bool betterConfigurationFound(false);

    while (this->GetNextSeedCandidate(&availableClusters, usedClusters, pTrialSeedCluster))
    {
        usedClusters.insert(pTrialSeedCluster);

        ClusterVector trialParticleSeedVector;
        trialParticleSeedVector.push_back(seedIter->first);
        trialParticleSeedVector.push_back(pTrialSeedCluster);

        SeedAssociationList trialSeedAssociationList;
        this->GetSeedAssociationList(trialParticleSeedVector, &availableClusters, trialSeedAssociationList);

        const float fom2(this->GetFigureOfMerit(trialSeedAssociationList));

        if (fom2 > fom1)
        {
            for (SeedAssociationList::const_iterator iter = trialSeedAssociationList.begin(), iterEnd = trialSeedAssociationList.end(); iter != iterEnd; ++iter)
            {
                this->CheckSeedAssociationList(iter, finalSeedAssociationList);
            }

            betterConfigurationFound = true;
            break;
        }
    }

    if (!betterConfigurationFound)
        finalSeedAssociationList.insert(SeedAssociationList::value_type(seedIter->first, seedIter->second));
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ClusterCharacterisationAlgorithm::GetFigureOfMerit(const SeedAssociationList &seedAssociationList) const
{
    const unsigned int nSeeds(seedAssociationList.size());

    if (!((nSeeds == 1) || (nSeeds == 2)))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    typedef std::set<const MCParticle*> ConstMCParticleList;
    ConstMCParticleList mcParticleList;

    for (SeedAssociationList::const_iterator iter1 = seedAssociationList.begin(), iter1End = seedAssociationList.end(); iter1 != iter1End; ++iter1)
    {
        Cluster *pParent(iter1->first);

        try // TODO don't cheat this bit!
        {
            const MCParticle *pMCParticle(MCParticleHelper::GetMainMCParticle(pParent));
            mcParticleList.insert(pMCParticle);
        }
        catch (StatusCodeException &)
        {
            mcParticleList.insert(NULL);
        }
    }

    if (mcParticleList.size() == nSeeds)
        return nSeeds;

    return 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ClusterCharacterisationAlgorithm::SortClusters(const Cluster *const pLhs, const Cluster *const pRhs)
{
    CartesianVector innerCoordinateLhs(0.f, 0.f, 0.f), outerCoordinateLhs(0.f, 0.f, 0.f);
    LArClusterHelper::GetExtremalCoordinatesXZ(pLhs, innerCoordinateLhs, outerCoordinateLhs);
    const float dLhs2((outerCoordinateLhs - innerCoordinateLhs).GetMagnitudeSquared());

    CartesianVector innerCoordinateRhs(0.f, 0.f, 0.f), outerCoordinateRhs(0.f, 0.f, 0.f);
    LArClusterHelper::GetExtremalCoordinatesXZ(pRhs, innerCoordinateRhs, outerCoordinateRhs);
    const float dRhs2((outerCoordinateRhs - innerCoordinateRhs).GetMagnitudeSquared());

    return (dLhs2 > dRhs2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListName", m_inputClusterListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
