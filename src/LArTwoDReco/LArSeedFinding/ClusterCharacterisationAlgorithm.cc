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

    PfoList pfoList;
    this->GetInputPfoList(pfoList);
    ClusterInfoMap nCaloHitsPerCluster, nBranchesPerCluster;

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

            this->StoreNCaloHitsPerCluster(iter->first, nCaloHitsPerCluster);
            this->StoreNBranchesPerCluster(iter->first, iter->second, nBranchesPerCluster);

            for (ClusterVector::const_iterator iter2 = iter->second.begin(), iter2End = iter->second.end(); iter2 != iter2End; ++iter2)
            {
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, iter->first, *iter2, m_inputClusterListName, m_inputClusterListName));
            }
        }
    }

    if (m_shouldRemoveShowerPfos)
    {
        this->RemoveShowerPfos(pClusterList, pfoList, nCaloHitsPerCluster, nBranchesPerCluster);
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

        if (pCluster->IsAvailable() && (pCluster->GetNCaloHits() < m_minCaloHitsPerCluster))
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

        if (!pCandidateCluster->IsAvailable() || seedClusters.count(pCandidateCluster) || (pCandidateCluster->GetNCaloHits() < m_minCaloHitsPerCluster))
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
    const float originalFigureOfMerit(this->GetFigureOfMerit(seedAssociationList));
    float bestFigureOfMerit(originalFigureOfMerit);

    Cluster *pTrialSeedCluster = NULL;
    bool betterConfigurationFound(false);
    SeedAssociationList bestSeedAssociationList;

    while (this->GetNextSeedCandidate(&availableClusters, usedClusters, pTrialSeedCluster))
    {
        usedClusters.insert(pTrialSeedCluster);

        ClusterVector trialParticleSeedVector;
        trialParticleSeedVector.push_back(seedIter->first);
        trialParticleSeedVector.push_back(pTrialSeedCluster);

        SeedAssociationList trialSeedAssociationList;
        this->GetSeedAssociationList(trialParticleSeedVector, &availableClusters, trialSeedAssociationList);
        const float trialFigureOfMerit(this->GetFigureOfMerit(trialSeedAssociationList));

        if (trialFigureOfMerit > bestFigureOfMerit)
        {
            betterConfigurationFound = true;
            bestFigureOfMerit = trialFigureOfMerit;
            bestSeedAssociationList = trialSeedAssociationList;
        }
    }

    if (betterConfigurationFound)
    {
        for (SeedAssociationList::const_iterator iter = bestSeedAssociationList.begin(), iterEnd = bestSeedAssociationList.end(); iter != iterEnd; ++iter)
        {
            this->CheckSeedAssociationList(iter, finalSeedAssociationList);
        }
    }
    else
    {
        finalSeedAssociationList.insert(SeedAssociationList::value_type(seedIter->first, seedIter->second));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ClusterCharacterisationAlgorithm::GetFigureOfMerit(const SeedAssociationList &seedAssociationList) const
{
    const unsigned int nSeeds(seedAssociationList.size());

    if (!((nSeeds == 1) || (nSeeds == 2)))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    unsigned int nMatchedClusters(0), nClusters(0);

    for (SeedAssociationList::const_iterator iter1 = seedAssociationList.begin(), iter1End = seedAssociationList.end(); iter1 != iter1End; ++iter1)
    {
        Cluster *pParentCluster(iter1->first);
        ++nClusters;

        // TODO don't cheat this bit!
        const MCParticle *pParentMCParticle(NULL);
        try {pParentMCParticle = MCParticleHelper::GetMainMCParticle(pParentCluster);} catch (StatusCodeException &) {}

        ++nMatchedClusters;
        const ClusterVector &associatedClusters(iter1->second);

        for (ClusterVector::const_iterator iter2 = associatedClusters.begin(), iter2End = associatedClusters.end(); iter2 != iter2End; ++iter2)
        {
            Cluster *pBranchCluster = *iter2;
            ++nClusters;

            const MCParticle *pDaughterMCParticle(NULL);
            try {pDaughterMCParticle = MCParticleHelper::GetMainMCParticle(pBranchCluster);} catch (StatusCodeException &) {}

            if (pParentMCParticle == pDaughterMCParticle)
            {
                ++nMatchedClusters;
            }
        }
    }

    if (0 == nClusters)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const float figureOfMerit(static_cast<float>(nMatchedClusters) / static_cast<float>(nClusters));
    return figureOfMerit;
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

void ClusterCharacterisationAlgorithm::GetInputPfoList(PfoList &pfoList) const
{
    for (StringVector::const_iterator iter = m_inputPfoListNames.begin(), iterEnd = m_inputPfoListNames.end(); iter != iterEnd; ++iter)
    {
        const PfoList *pPfoList = NULL;

        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, *iter, pPfoList))
        {
            std::cout << "ClusterCharacterisationAlgorithm : Could not find input pfo list with name " << *iter << std::endl;
            continue;
        }

        pfoList.insert(pPfoList->begin(), pPfoList->end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterCharacterisationAlgorithm::StoreNCaloHitsPerCluster(const Cluster *const pCluster, ClusterInfoMap &clusterInfoMap) const
{
    // ATTN: Stores only first value provided per cluster
    if (clusterInfoMap.count(pCluster))
        return;

    if (!clusterInfoMap.insert(ClusterInfoMap::value_type(pCluster, pCluster->GetNCaloHits())).second)
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterCharacterisationAlgorithm::StoreNBranchesPerCluster(const Cluster *const pCluster, const ClusterVector &branchList,
    ClusterInfoMap &clusterInfoMap) const
{
    // ATTN: Stores total number of branches linked to eventual parent cluster
    ClusterInfoMap::const_iterator iIter = clusterInfoMap.find(pCluster);
    unsigned int nBranchesSum((clusterInfoMap.end() == iIter) ? 0 : iIter->second);

    for (ClusterVector::const_iterator iter = branchList.begin(), iterEnd = branchList.end(); iter != iterEnd; ++iter)
    {
        Cluster *pBranchCluster = *iter;
        ClusterInfoMap::const_iterator bIter = clusterInfoMap.find(pBranchCluster);
        const unsigned int nBranches((clusterInfoMap.end() == bIter) ? 0 : bIter->second);
        nBranchesSum += (1 + nBranches);
    }

    clusterInfoMap[pCluster] = nBranchesSum;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterCharacterisationAlgorithm::RemoveShowerPfos(const ClusterList *const pClusterList, PfoList &pfoList,
    const ClusterInfoMap &nCaloHitsPerCluster, const ClusterInfoMap &nBranchesPerCluster) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        try
        {
            Cluster *pCluster = *iter;

            if (pCluster->IsAvailable())
                continue;

            ClusterInfoMap::const_iterator nCaloHitsIter = nCaloHitsPerCluster.find(pCluster);
            ClusterInfoMap::const_iterator nBranchesIter = nBranchesPerCluster.find(pCluster);

            if ((nCaloHitsPerCluster.end() == nCaloHitsIter) || (nBranchesPerCluster.end() == nBranchesIter))
                continue;

            if (0 == nCaloHitsIter->second)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            const float nCaloHitsRatio(static_cast<float>(pCluster->GetNCaloHits()) / static_cast<float>(nCaloHitsIter->second));
            const unsigned int nBranches(nBranchesIter->second);

            if ((nBranches < m_showerLikeNBranches) && (nCaloHitsRatio < m_showerLikeCaloHitRatio))
                continue;

            Pfo *pTargetPfo = NULL;
            this->FindTargetPfo(pCluster, pfoList, pTargetPfo);

            pfoList.erase(pTargetPfo);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pTargetPfo));
        }
        catch (StatusCodeException &)
        {
            std::cout << "ClusterCharacterisationAlgorithm: Unable to remove shower-like particle flow object " << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterCharacterisationAlgorithm::FindTargetPfo(Cluster *const pCluster, const PfoList &pfoList, Pfo *&pTargetPfo) const
{
    pTargetPfo = NULL;

    for (PfoList::const_iterator iter = pfoList.begin(), iterEnd = pfoList.end(); iter != iterEnd; ++iter)
    {
        Pfo *pPfo = *iter;

        if (pPfo->GetClusterList().count(pCluster))
        {
            pTargetPfo = pPfo;
            return;
        }
    }

    if (NULL == pTargetPfo)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListName", m_inputClusterListName));

    m_inputPfoListNames.clear();
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "InputPfoListNames", m_inputPfoListNames));

    m_minCaloHitsPerCluster = 5;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCaloHitsPerCluster", m_minCaloHitsPerCluster));

    m_shouldRemoveShowerPfos = true;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShouldRemoveShowerPfos", m_shouldRemoveShowerPfos));

    m_showerLikeNBranches = 5;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowerLikeNBranches", m_showerLikeNBranches));

    m_showerLikeCaloHitRatio = 2.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowerLikeCaloHitRatio", m_showerLikeCaloHitRatio));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
