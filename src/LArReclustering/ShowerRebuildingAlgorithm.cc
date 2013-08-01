/**
 *  @file   LArContent/src/LArReclustering/ShowerRebuildingAlgorithm.cc
 * 
 *  @brief  Implementation of the shower rebuilding algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArThreeDHelper.h"

#include "LArReclustering/ShowerRebuildingAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode ShowerRebuildingAlgorithm::Run()
{
    // Rebuild 2d showers relating to seeds that have been added to 3d particles
    this->RebuildThreeDShowers();

    // Rebuild any 2d showers for which 2d->3d matching was not possible
    this->RebuildTwoDShowers(m_seedClusterListNameU, m_nonSeedClusterListNameU);
    this->RebuildTwoDShowers(m_seedClusterListNameV, m_nonSeedClusterListNameV);
    this->RebuildTwoDShowers(m_seedClusterListNameW, m_nonSeedClusterListNameW);

    // Move any seed showers that did not have promising spines back into the shower seed list
    this->RestoreLoneShowers();

    LArThreeDHelper::RemoveAllStoredClusters();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerRebuildingAlgorithm::RebuildThreeDShowers() const
{
    try
    {
        const PfoList *pPfoList = NULL;
        PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetPfoList(*this, m_particleListName, pPfoList));

        if (NULL == pPfoList)
        {
            std::cout << "ShowerRebuildingAlgorithm: pfo list " << m_particleListName << " unavailable." << std::endl;
            throw StatusCodeException(STATUS_CODE_SUCCESS);
        }

        PfoVector pfoVector(pPfoList->begin(), pPfoList->end());
        std::sort(pfoVector.begin(), pfoVector.end(), ShowerRebuildingAlgorithm::SortByNHits);

        for (PfoVector::iterator iter = pfoVector.begin(), iterEnd = pfoVector.end(); iter != iterEnd; ++iter)
        {
            ParticleFlowObject *pPfo = *iter;

            if (NULL == pPfo)
                continue;

            *iter = NULL;

            try
            {
                Cluster *pSeedClusterU(NULL), *pSeedClusterV(NULL), *pSeedClusterW(NULL);
                this->GetSeedClusters(pPfo, pSeedClusterU, pSeedClusterV, pSeedClusterW);

                ClusterList seedMergesU, seedMergesV, seedMergesW, nonSeedMergesU, nonSeedMergesV, nonSeedMergesW;
                this->RebuildThreeDShower(pfoVector, pSeedClusterU, pSeedClusterV, pSeedClusterW, seedMergesU, nonSeedMergesU);
                this->RebuildThreeDShower(pfoVector, pSeedClusterV, pSeedClusterW, pSeedClusterU, seedMergesV, nonSeedMergesV);
                this->RebuildThreeDShower(pfoVector, pSeedClusterW, pSeedClusterU, pSeedClusterV, seedMergesW, nonSeedMergesW);

                this->PerformClusterMerges(pSeedClusterU, seedMergesU, m_seedClusterListNameU, m_seedClusterListNameU);
                this->PerformClusterMerges(pSeedClusterV, seedMergesV, m_seedClusterListNameV, m_seedClusterListNameV);
                this->PerformClusterMerges(pSeedClusterW, seedMergesW, m_seedClusterListNameW, m_seedClusterListNameW);
                this->PerformClusterMerges(pSeedClusterU, nonSeedMergesU, m_seedClusterListNameU, m_nonSeedClusterListNameU);
                this->PerformClusterMerges(pSeedClusterV, nonSeedMergesV, m_seedClusterListNameV, m_nonSeedClusterListNameV);
                this->PerformClusterMerges(pSeedClusterW, nonSeedMergesW, m_seedClusterListNameW, m_nonSeedClusterListNameW);
            }
            catch (StatusCodeException &statusCodeException)
            {
                if (STATUS_CODE_NOT_ALLOWED != statusCodeException.GetStatusCode())
                    throw statusCodeException;

                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::DeletePfo(*this, pPfo));
            }
        }
    }
    catch (StatusCodeException &statusCodeException)
    {
        if (STATUS_CODE_SUCCESS != statusCodeException.GetStatusCode())
            std::cout << "ShowerRebuildingAlgorithm::RebuildThreeDShowers, " << statusCodeException.ToString() << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerRebuildingAlgorithm::RebuildThreeDShower(PfoVector &pfoVector, Cluster *pSeedCluster, Cluster *pSeedCluster2, Cluster *pSeedCluster3,
    ClusterList &seedMerges, ClusterList &nonSeedMerges) const
{
    ClusterList seeds, nonSeeds;
    LArThreeDHelper::GetAllSeedComponents(pSeedCluster, seeds);
    LArThreeDHelper::GetAllNonSeedComponents(pSeedCluster, nonSeeds);

    const unsigned int seedId2(LArThreeDHelper::GetClusterIdFromSeed(pSeedCluster2));
    const unsigned int seedId3(LArThreeDHelper::GetClusterIdFromSeed(pSeedCluster3));

    // Merge-in the associated seeds
    for (ClusterList::const_iterator iter = seeds.begin(), iterEnd = seeds.end(); iter != iterEnd; ++iter)
    {
        Cluster *pSiblingSeed = *iter;

        if (pSeedCluster == pSiblingSeed)
            continue;

        if (!pSiblingSeed->IsAvailable())
        {
            for (PfoVector::iterator pIter = pfoVector.begin(), pIterEnd = pfoVector.end(); pIter != pIterEnd; ++pIter)
            {
                ParticleFlowObject *pPfo = *pIter;

                if (NULL == pPfo)
                    continue;

                Cluster *pSeedClusterU(NULL), *pSeedClusterV(NULL), *pSeedClusterW(NULL);
                this->GetSeedClusters(pPfo, pSeedClusterU, pSeedClusterV, pSeedClusterW);

                const bool isMatchU(pSiblingSeed == pSeedClusterU);
                const bool isMatchV(pSiblingSeed == pSeedClusterV);
                const bool isMatchW(pSiblingSeed == pSeedClusterW);

                if (!isMatchU && !isMatchV && !isMatchW)
                    continue;

                const unsigned int clusterIdU(LArThreeDHelper::GetClusterIdFromSeed(pSeedClusterU));
                const unsigned int clusterIdV(LArThreeDHelper::GetClusterIdFromSeed(pSeedClusterV));
                const unsigned int clusterIdW(LArThreeDHelper::GetClusterIdFromSeed(pSeedClusterW));

                const bool isSameShower(
                    (isMatchU && ((seedId2 == clusterIdV && seedId3 == clusterIdW) || (seedId2 == clusterIdW && seedId3 == clusterIdV))) ||
                    (isMatchV && ((seedId2 == clusterIdU && seedId3 == clusterIdW) || (seedId2 == clusterIdW && seedId3 == clusterIdU))) ||
                    (isMatchW && ((seedId2 == clusterIdU && seedId3 == clusterIdV) || (seedId2 == clusterIdV && seedId3 == clusterIdU))) );

                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::DeletePfo(*this, pPfo));
                *pIter = NULL;

                if (!isSameShower)
                    throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);
             }
        }

        seedMerges.insert(pSiblingSeed);
    }

    // Merge-in the associated non-seeds
    for (ClusterList::const_iterator iter = nonSeeds.begin(), iterEnd = nonSeeds.end(); iter != iterEnd; ++iter)
    {
        nonSeedMerges.insert(*iter);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerRebuildingAlgorithm::GetSeedClusters(const ParticleFlowObject *const pPfo, Cluster *&pClusterU, Cluster *&pClusterV, Cluster *&pClusterW) const
{
    if (3 != pPfo->GetClusterList().size())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    pClusterU = NULL; pClusterV = NULL; pClusterW = NULL;

    for (ClusterList::const_iterator cIter = pPfo->GetClusterList().begin(), cIterEnd = pPfo->GetClusterList().end(); cIter != cIterEnd; ++cIter)
    {
        Cluster *pCluster = *cIter;
        const bool containsU(pCluster->ContainsHitType(VIEW_U));
        const bool containsV(pCluster->ContainsHitType(VIEW_V));
        const bool containsW(pCluster->ContainsHitType(VIEW_W));

        if (containsU && !containsV && !containsW)
        {
            pClusterU = pCluster;
        }
        else if (!containsU && containsV && !containsW)
        {
            pClusterV = pCluster;
        }
        else if (!containsU && !containsV && containsW)
        {
            pClusterW = pCluster;
        }
        else
        {
            throw StatusCodeException(STATUS_CODE_FAILURE);
        }
    }

    if ((NULL == pClusterU) || (NULL == pClusterV) || (NULL == pClusterW))
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerRebuildingAlgorithm::PerformClusterMerges(Cluster *pCluster, const ClusterList &clusterList, const std::string &clusterListName1,
    const std::string &clusterListName2) const
{
    for (ClusterList::const_iterator iter = clusterList.begin(), iterEnd = clusterList.end(); iter != iterEnd; ++iter)
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pCluster, *iter, clusterListName1, clusterListName2));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerRebuildingAlgorithm::RebuildTwoDShowers(const std::string &seedClusterListName, const std::string &nonSeedClusterListName) const
{
    try
    {
        const ClusterList *pClusterList = NULL;
        PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetClusterList(*this, seedClusterListName, pClusterList));

        if (NULL == pClusterList)
        {
            std::cout << "ShowerRebuildingAlgorithm: cluster list " << seedClusterListName << " unavailable." << std::endl;
            throw StatusCodeException(STATUS_CODE_SUCCESS);
        }

        for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
        {
            try
            {
                Cluster *pSeedCluster = *iter;

                if (!pSeedCluster->IsAvailable())
                    continue;

                ClusterList seeds, nonSeeds;
                LArThreeDHelper::GetAllSeedComponents(pSeedCluster, seeds);
                LArThreeDHelper::GetAllNonSeedComponents(pSeedCluster, nonSeeds);

                // Merge-in the associated seeds
                for (ClusterList::const_iterator cIter = seeds.begin(), cIterEnd = seeds.end(); cIter != cIterEnd; ++cIter)
                {
                    if (pSeedCluster == *cIter)
                        continue;

                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pSeedCluster, *cIter, seedClusterListName, seedClusterListName));
                }

                // Merge-in the associated non-seeds
                for (ClusterList::const_iterator cIter = nonSeeds.begin(), cIterEnd = nonSeeds.end(); cIter != cIterEnd; ++cIter)
                {
                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pSeedCluster, *cIter, seedClusterListName, nonSeedClusterListName));
                }
            }
            catch (StatusCodeException &statusCodeException)
            {
                if (STATUS_CODE_NOT_FOUND != statusCodeException.GetStatusCode())
                    throw statusCodeException;
            }
        }
    }
    catch (StatusCodeException &statusCodeException)
    {
        if (STATUS_CODE_SUCCESS != statusCodeException.GetStatusCode())
            std::cout << "ShowerRebuildingAlgorithm::RebuildTwoDShowers, " << statusCodeException.ToString() << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerRebuildingAlgorithm::RestoreLoneShowers() const
{
    try
    {
        const ClusterList &loneClusters(LArThreeDHelper::GetLoneClusterList());

        if (loneClusters.empty())
            return;

        ClusterList loneClustersU, loneClustersV, loneClustersW;

        for (ClusterList::const_iterator iter = loneClusters.begin(), iterEnd = loneClusters.end(); iter != iterEnd; ++iter)
        {
            Cluster *pCluster = *iter;
            const bool containsU(pCluster->ContainsHitType(VIEW_U));
            const bool containsV(pCluster->ContainsHitType(VIEW_V));
            const bool containsW(pCluster->ContainsHitType(VIEW_W));

            if (containsU && !containsV && !containsW)
            {
                loneClustersU.insert(pCluster);
            }
            else if (!containsU && containsV && !containsW)
            {
                loneClustersV.insert(pCluster);
            }
            else if (!containsU && !containsV && containsW)
            {
                loneClustersW.insert(pCluster);
            }
            else
            {
                throw StatusCodeException(STATUS_CODE_FAILURE);
            }
        }

        if (!loneClustersU.empty())
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveClusterList(*this, m_nonSeedClusterListNameU, m_seedClusterListNameU, loneClustersU));

        if (!loneClustersV.empty())
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveClusterList(*this, m_nonSeedClusterListNameV, m_seedClusterListNameV, loneClustersV));

        if (!loneClustersW.empty())
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveClusterList(*this, m_nonSeedClusterListNameW, m_seedClusterListNameW, loneClustersW));
    }
    catch (StatusCodeException &statusCodeException)
    {
        if (STATUS_CODE_SUCCESS != statusCodeException.GetStatusCode())
            std::cout << "ShowerRebuildingAlgorithm::RestoreLoneShowers, " << statusCodeException.ToString() << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ShowerRebuildingAlgorithm::SortByNHits(const ParticleFlowObject *const pLhs, const ParticleFlowObject *const pRhs)
{
    unsigned int nHitsLhs(0);
    const ClusterList &clusterListLhs(pLhs->GetClusterList());

    for (ClusterList::const_iterator iter = clusterListLhs.begin(), iterEnd = clusterListLhs.end(); iter != iterEnd; ++iter)
        nHitsLhs += (*iter)->GetNCaloHits();

    unsigned int nHitsRhs(0);
    const ClusterList &clusterListRhs(pRhs->GetClusterList());

    for (ClusterList::const_iterator iter = clusterListRhs.begin(), iterEnd = clusterListRhs.end(); iter != iterEnd; ++iter)
        nHitsRhs += (*iter)->GetNCaloHits();

    return (nHitsLhs > nHitsRhs);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerRebuildingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ParticleListName", m_particleListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "SeedClusterListNameU", m_seedClusterListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "SeedClusterListNameV", m_seedClusterListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "SeedClusterListNameW", m_seedClusterListNameW));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NonSeedClusterListNameU", m_nonSeedClusterListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NonSeedClusterListNameV", m_nonSeedClusterListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NonSeedClusterListNameW", m_nonSeedClusterListNameW));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
