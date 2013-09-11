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

    LArThreeDHelper::Reset();

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
        std::sort(pfoVector.begin(), pfoVector.end(), ShowerRebuildingAlgorithm::SortPfosByNHits);

        for (PfoVector::iterator iter = pfoVector.begin(), iterEnd = pfoVector.end(); iter != iterEnd; ++iter)
        {
            ParticleFlowObject *pPfo = *iter;

            if (NULL == pPfo)
                continue;

            try
            {
                Cluster *pClusterU(NULL), *pClusterV(NULL), *pClusterW(NULL);
                this->GetSeedClusters(pPfo, pClusterU, pClusterV, pClusterW);

                ClusterList seedMergesU, seedMergesV, seedMergesW, nonSeedMergesU, nonSeedMergesV, nonSeedMergesW;
                if (pClusterU) this->RebuildThreeDShower(pPfo, pClusterU, pfoVector, seedMergesU, nonSeedMergesU);
                if (pClusterV) this->RebuildThreeDShower(pPfo, pClusterV, pfoVector, seedMergesV, nonSeedMergesV);
                if (pClusterW) this->RebuildThreeDShower(pPfo, pClusterW, pfoVector, seedMergesW, nonSeedMergesW);

                if (pClusterU) this->PerformClusterMerges(pClusterU, seedMergesU, m_seedClusterListNameU, m_seedClusterListNameU);
                if (pClusterV) this->PerformClusterMerges(pClusterV, seedMergesV, m_seedClusterListNameV, m_seedClusterListNameV);
                if (pClusterW) this->PerformClusterMerges(pClusterW, seedMergesW, m_seedClusterListNameW, m_seedClusterListNameW);

                if (pClusterU) this->PerformClusterMerges(pClusterU, nonSeedMergesU, m_seedClusterListNameU, m_nonSeedClusterListNameU);
                if (pClusterV) this->PerformClusterMerges(pClusterV, nonSeedMergesV, m_seedClusterListNameV, m_nonSeedClusterListNameV);
                if (pClusterW) this->PerformClusterMerges(pClusterW, nonSeedMergesW, m_seedClusterListNameW, m_nonSeedClusterListNameW);
            }
            catch (StatusCodeException &statusCodeException)
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::DeletePfo(*this, pPfo));
            }

            *iter = NULL;
        }
    }
    catch (StatusCodeException &statusCodeException)
    {
        if (STATUS_CODE_SUCCESS != statusCodeException.GetStatusCode())
            std::cout << "ShowerRebuildingAlgorithm::RebuildThreeDShowers, " << statusCodeException.ToString() << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerRebuildingAlgorithm::GetSeedClusters(const ParticleFlowObject *const pPfo, Cluster *&pSeedClusterU, Cluster *&pSeedClusterV,
    Cluster *&pSeedClusterW) const
{
    ClusterVector seedClustersU, seedClustersV, seedClustersW;

    for (ClusterList::const_iterator cIter = pPfo->GetClusterList().begin(), cIterEnd = pPfo->GetClusterList().end(); cIter != cIterEnd; ++cIter)
    {
        Cluster *pCluster = *cIter;
        const bool containsU(pCluster->ContainsHitType(VIEW_U));
        const bool containsV(pCluster->ContainsHitType(VIEW_V));
        const bool containsW(pCluster->ContainsHitType(VIEW_W));

        if (containsU && !containsV && !containsW)
        {
            seedClustersU.push_back(pCluster);
        }
        else if (!containsU && containsV && !containsW)
        {
            seedClustersV.push_back(pCluster);
        }
        else if (!containsU && !containsV && containsW)
        {
            seedClustersW.push_back(pCluster);
        }
        else
        {
            throw StatusCodeException(STATUS_CODE_FAILURE);
        }
    }

    if ((seedClustersU.empty() && seedClustersV.empty()) ||
        (seedClustersV.empty() && seedClustersW.empty()) ||
        (seedClustersW.empty() && seedClustersU.empty()) )
    {
        throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);
    }

    pSeedClusterU = NULL;
    pSeedClusterV = NULL;
    pSeedClusterW = NULL;

    if (!seedClustersU.empty())
    {
        std::sort(seedClustersU.begin(), seedClustersU.end(), ShowerRebuildingAlgorithm::SortClustersByNHits);
        pSeedClusterU = seedClustersU.at(0);
    }

    if (!seedClustersV.empty())
    {
        std::sort(seedClustersV.begin(), seedClustersV.end(), ShowerRebuildingAlgorithm::SortClustersByNHits);
        pSeedClusterV = seedClustersV.at(0);
    }

    if (!seedClustersW.empty())
    {
        std::sort(seedClustersW.begin(), seedClustersW.end(), ShowerRebuildingAlgorithm::SortClustersByNHits);
        pSeedClusterW = seedClustersW.at(0);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerRebuildingAlgorithm::RebuildThreeDShower(const ParticleFlowObject *const pPfo, Cluster *pSeedCluster, PfoVector &pfoVector,
    ClusterList &seedMerges, ClusterList &nonSeedMerges) const
{
    this->GetAssociatedClusters(pSeedCluster, seedMerges, nonSeedMerges);

    for (ClusterList::const_iterator iter = seedMerges.begin(), iterEnd = seedMerges.end(); iter != iterEnd; ++iter)
    {
        Cluster *pSiblingSeed = *iter;

        if (pSeedCluster == pSiblingSeed)
            continue;

        if (pSiblingSeed->IsAvailable())
            continue;

        // First pfo in sorted container wins in the event of any mismatches with initial 2D reconstruction
        for (PfoVector::iterator pfoIter = pfoVector.begin(), pfoIterEnd = pfoVector.end(); pfoIter != pfoIterEnd; ++pfoIter)
        {
            ParticleFlowObject *pAlternativePfo = *pfoIter;

            if (NULL == pAlternativePfo)
                continue;

            if (pAlternativePfo->GetClusterList().count(pSiblingSeed) == 0)
                continue;

            if (pAlternativePfo->GetNClusters() > 1)
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveClusterFromPfo(*this, pAlternativePfo, pSiblingSeed));
            }
            else
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::DeletePfo(*this, pAlternativePfo));
                *pfoIter = NULL;
            }

            /*
            if (pPfo == pAlternativePfo)
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveClusterFromPfo(*this, pAlternativePfo, pSiblingSeed));
            }
            else
            {
                std::cout << "Shower spines matched between multiple clusters - delete PFO with fewest hits " << std::endl;
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::DeletePfo(*this, pAlternativePfo));
                *pfoIter = NULL;
            }
            */
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerRebuildingAlgorithm::GetAssociatedClusters(const Cluster *const pCluster, ClusterList &seedMerges, ClusterList &nonSeedMerges) const
{
    ClusterList siblingSeeds, nonSeeds;
    LArThreeDHelper::GetAllSeedComponents(pCluster, siblingSeeds);
    LArThreeDHelper::GetAllNonSeedComponents(pCluster, nonSeeds);

    for (ClusterList::const_iterator iter = siblingSeeds.begin(), iterEnd = siblingSeeds.end(); iter != iterEnd; ++iter)
    {
        if ((pCluster != *iter) && (seedMerges.end() == seedMerges.find(*iter)))
        {
            seedMerges.insert(*iter);
            this->GetAssociatedClusters(*iter, seedMerges, nonSeedMerges);
        }
    }

    for (ClusterList::const_iterator iter = nonSeeds.begin(), iterEnd = nonSeeds.end(); iter != iterEnd; ++iter)
    {
        nonSeedMerges.insert(*iter);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerRebuildingAlgorithm::PerformClusterMerges(Cluster *pCluster, const ClusterList &clusterList, const std::string &clusterListName1,
    const std::string &clusterListName2) const
{
    for (ClusterList::const_iterator iter = clusterList.begin(), iterEnd = clusterList.end(); iter != iterEnd; ++iter)
    {
        if (*iter == pCluster)
            continue;

        if (!(*iter)->IsAvailable())
            continue;

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pCluster, *iter, clusterListName1, clusterListName2));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerRebuildingAlgorithm::RebuildTwoDShowers(const std::string &seedClusterListName, const std::string &nonSeedClusterListName) const
{
    try
    {
        const ClusterList *pSeedClusterList = NULL;
        const ClusterList *pNonSeedClusterList = NULL;
        PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetClusterList(*this, seedClusterListName, pSeedClusterList));
        PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetClusterList(*this, nonSeedClusterListName, pNonSeedClusterList));

        if (NULL == pSeedClusterList)
        {
            std::cout << "ShowerRebuildingAlgorithm: cluster list " << seedClusterListName << " unavailable." << std::endl;
            throw StatusCodeException(STATUS_CODE_SUCCESS);
        }

        for (ClusterList::const_iterator iter = pSeedClusterList->begin(), iterEnd = pSeedClusterList->end(); iter != iterEnd; ++iter)
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

                    if (pSeedClusterList->count(*cIter) == 0)
                        continue;

                    if (!(*cIter)->IsAvailable())
                        continue;

                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pSeedCluster, *cIter, seedClusterListName, seedClusterListName));
                }

                if ((NULL == pNonSeedClusterList) && !nonSeeds.empty())
                {
                    std::cout << "ShowerRebuildingAlgorithm: cluster list " << nonSeedClusterListName << " unavailable." << std::endl;
                    continue;
                }

                // Merge-in the associated non-seeds
                for (ClusterList::const_iterator cIter = nonSeeds.begin(), cIterEnd = nonSeeds.end(); cIter != cIterEnd; ++cIter)
                {
                    if (pSeedCluster == *cIter)
                        continue;

                    if (pNonSeedClusterList->count(*cIter) == 0)
                        continue;

                    if (!(*cIter)->IsAvailable())
                        continue;

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

bool ShowerRebuildingAlgorithm::SortPfosByNHits(const ParticleFlowObject *const pLhs, const ParticleFlowObject *const pRhs)
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

bool ShowerRebuildingAlgorithm::SortClustersByNHits(const Cluster *const pLhs, const Cluster *const pRhs)
{
    return (pLhs->GetNCaloHits() > pRhs->GetNCaloHits());
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
