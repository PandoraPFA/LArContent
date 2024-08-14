/**
 * @file   larpandoracontent/LArCheating/CheatingClusterMergingAlgorithm.cc
 *
 * @brief  Implementation file for the cheating cluster merging algorithm.
 *
 * $Log: 
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Helpers/MCParticleHelper.h"
#include "Objects/MCParticle.h"

#include <cmath>

#include "larpandoracontent/LArCheating/CheatingClusterMergingAlgorithm.h"
#include "larpandoracontent/LArTrackShowerId/ShowerGrowingAlgorithm.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"

using namespace pandora;

namespace lar_content
{

CheatingClusterMergingAlgorithm::CheatingClusterMergingAlgorithm() :
    m_maxClusterFraction(0.25f),
    m_minNCaloHits(1),
    m_writeTree{true}
{
    if (m_maxClusterFraction < 0)
        std::cout << "WARN: m_maxClusterFraction can not be below 0!" << std::endl;
}
//------------------------------------------------------------------------------------------//

StatusCode CheatingClusterMergingAlgorithm::Run()
{
    for (const std::string &clusterListName : m_inputClusterListNames)
    {
        try
        {
            const ClusterList *pClusterList = nullptr;
            PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));

            if (!pClusterList || pClusterList->empty())
            {
                if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                    std::cout << "ClusterMergingAlgorithm: unable to find cluster list " << clusterListName << std::endl;

                continue;
            }

            this->CheatedClusterMerging(pClusterList, clusterListName);
            std::cout << "---------------------------------------------------------------------------------" << std::endl;   
     
        }
        
        catch (StatusCodeException &statusCodeException)
        {
            throw statusCodeException;
        }
        
        if (m_writeTree)
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "RECREATE"))
         }

    }

    return STATUS_CODE_SUCCESS;
}

//--------------------------------------------------------------------------------------------//
const MCParticle* CheatingClusterMergingAlgorithm::GetMCForCluster(const Cluster *const cluster, std::map<const Cluster*,
    const MCParticle*> &clusterToMCMap) const
{
    const MCParticle* clusterMC = nullptr;

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
            std::cout << "Failed to get MC particle for cluster of " << cluster->GetOrderedCaloHitList().size()
                      << " : " << e.ToString() << std::endl;
        }
    }

    return clusterMC;
}

//---------------------------------------------------------------------------------------------//

bool CheatingClusterMergingAlgorithm::IsValidToUse(const Cluster *const cluster, std::map<const Cluster*, bool> &clusterIsUsed) const
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
    std::map<const Cluster*, const MCParticle*> clusterToMCParticleMap;

    std::map<const Cluster*, bool> clusterIsUsed;
    std::map<const Cluster*, ClusterVector> clustersToMerge;

    for (auto it = pClusterList->begin(); it != pClusterList->end(); ++it)
    {
        const Cluster *cluster(*it);

        if (!this->IsValidToUse(cluster, clusterIsUsed))
            continue;

        const MCParticle *clusterMC(this->GetMCForCluster(cluster, clusterToMCParticleMap));

        for (auto it2 = std::next(it); it2 != pClusterList->end(); ++it2) {
            const Cluster *const otherCluster(*it2);

            if (!this->IsValidToUse(otherCluster, clusterIsUsed))
                continue;

            const MCParticle *otherClusterMC(this->GetMCForCluster(otherCluster, clusterToMCParticleMap));

            if (clusterMC == otherClusterMC)
            {
                clusterIsUsed[cluster] = true;
                clusterIsUsed[otherCluster] = true;
                clustersToMerge[cluster].push_back(otherCluster);
            }    
        }
    }

 // What to do with clusters that are never merged?
 //     - Have a fallback check for second main MC, if over X%?
 
 for (auto clusterToMergePair : clustersToMerge)
    {
        const Cluster *currentCluster = clusterToMergePair.first;
        const auto clusters = clusterToMergePair.second;

        for (auto clusterToMerge : clusters)
        {
            if (! clusterToMerge->IsAvailable())
                continue;

            try
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, currentCluster, clusterToMerge, listName, listName));
            } catch (StatusCodeException) {}
        }
    }

    return;
}

//-------------------------------------------------------------------------------------------------

float CheatingClusterMergingAlgorithm::Distance(const CartesianVector vector1, const CartesianVector vector2) const
{
    const float dx{vector1.GetX() - vector2.GetX()};
    const float dz{vector1.GetZ() - vector2.GetZ()};
    float distance{std::sqrt(dx * dx + dz * dz)};
 
    return distance;
}

//-------------------------------------------------------------------------------------------------

double CheatingClusterMergingAlgorithm::Angle(const CartesianVector vector1, const CartesianVector vector2) const
{
    const double dx{vector1.GetX() - vector2.GetX()};
    const double dz{vector1.GetZ() - vector2.GetZ()};
    double angle{tan(dx / dz)};

    return angle;
}
 
//----------------------------------------------------------------------------------------------------

StatusCode CheatingClusterMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputClusterListNames", m_inputClusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MCParticleListName", m_mcParticleListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxClusterFraction", m_maxClusterFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinNCaloHits", m_minNCaloHits));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TreeName", m_treeName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "FileName", m_fileName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteTree", m_writeTree));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
