/**
 * @file   larpandoracontent/LArCheating/CheatingClusterMergingAlgorithm.cc
 *
 * @brief  Implementation file for the cheating cluster merging algorithm.
 *
 * $Log:
 */

#include "Helpers/MCParticleHelper.h"
#include "Objects/MCParticle.h"
#include "Pandora/AlgorithmHeaders.h"

#include <cmath>

#include "larpandoracontent/LArCheating/CheatingClusterMergingAlgorithm.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"
#include "larpandoracontent/LArTrackShowerId/ShowerGrowingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CheatingClusterMergingAlgorithm::CheatingClusterMergingAlgorithm() :
    m_minNCaloHits(1)
{
}
//------------------------------------------------------------------------------------------//

StatusCode CheatingClusterMergingAlgorithm::Run()
{
    for (const std::string &clusterListName : m_inputClusterListNames)
    {
        try
        {
            const ClusterList *pClusterList{nullptr};
            PANDORA_RETURN_RESULT_IF_AND_IF(
                STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));

            if (!pClusterList || pClusterList->empty())
            {
                if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                    std::cout << "ClusterMergingAlgorithm: unable to find cluster list " << clusterListName << std::endl;

                continue;
            }
            this->CheatedClusterMerging(pClusterList, clusterListName);
        }
        catch (StatusCodeException &statusCodeException)
        {
            throw statusCodeException;
        }
    }
    return STATUS_CODE_SUCCESS;
}
//--------------------------------------------------------------------------------------------//
const MCParticle *CheatingClusterMergingAlgorithm::GetMCForCluster(
    const Cluster *const cluster, std::map<const Cluster *, const MCParticle *> &clusterToMCMap) const
{
    const MCParticle *clusterMC{nullptr};

    if (clusterToMCMap.count(cluster) > 0)
    {
        clusterMC = clusterToMCMap.at(cluster);
    }
    else
    {
        try
        {
            clusterMC = MCParticleHelper::GetMainMCParticle(cluster);
            clusterToMCMap[cluster] = clusterMC;
        }
        catch (StatusCodeException e)
        {
            std::cout << "Failed to get MC particle for cluster of " << cluster->GetOrderedCaloHitList().size() << " : " << e.ToString() << std::endl;
        }
    }
    return clusterMC;
}

//---------------------------------------------------------------------------------------------//

bool CheatingClusterMergingAlgorithm::IsValidToUse(const Cluster *const cluster, std::map<const Cluster *, bool> &clusterIsUsed) const
{
    if (!cluster->IsAvailable())
        return false;

    if (cluster->GetNCaloHits() < m_minNCaloHits)
        return false;

    if (clusterIsUsed.count(cluster) > 0)
        return false;

    return true;
}

//-----------------------------------------------------------------------------------------------//

void CheatingClusterMergingAlgorithm::CheatedClusterMerging(const pandora::ClusterList *const pClusterList, const std::string &listName) const
{
    std::map<const Cluster *, const MCParticle *> clusterToMCParticleMap;
    std::map<const Cluster *, bool> clusterIsUsed;
    std::map<const Cluster *, ClusterVector> clustersToMerge;

    for (auto it = pClusterList->begin(); it != pClusterList->end(); ++it)
    {
        const Cluster *cluster(*it);

        if (!this->IsValidToUse(cluster, clusterIsUsed))
            continue;

        const MCParticle *clusterMC(this->GetMCForCluster(cluster, clusterToMCParticleMap));

        for (auto it2 = std::next(it); it2 != pClusterList->end(); ++it2)
        {
            const Cluster *const otherCluster(*it2);

            if (!this->IsValidToUse(otherCluster, clusterIsUsed))
                continue;

            const MCParticle *otherClusterMC(this->GetMCForCluster(otherCluster, clusterToMCParticleMap));

            if (clusterMC == otherClusterMC)
            {
                clusterIsUsed[cluster] = true;
                clusterIsUsed[otherCluster] = true;
                clustersToMerge[cluster].emplace_back(otherCluster);
            }
        }
    }

    for (const auto &[pCurrentCluster, clusters] : clustersToMerge)
    {
        for (const Cluster *const pClusterToMerge : clusters)
        {
            if (!pClusterToMerge->IsAvailable())
                continue;

            try
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                    PandoraContentApi::MergeAndDeleteClusters(*this, pCurrentCluster, pClusterToMerge, listName, listName));
            }
            catch (StatusCodeException)
            {
            }
        }
    }
}

//-------------------------------------------------------------------------------------------------

StatusCode CheatingClusterMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputClusterListNames", m_inputClusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinNCaloHits", m_minNCaloHits));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
