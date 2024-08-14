/**
 * @file   larpandoracontent/LArMonitoring/TopologyCalorimetryMonitoringAlgorithm.cc
 *
 * @brief  Implementation file for the topology and calorimetry monitoring algorithm.
 *
 * $Log: 
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Helpers/MCParticleHelper.h"
#include "Objects/MCParticle.h"

#include "larpandoracontent/LArMonitoring/TopologyCalorimetryMonitoringAlgorithm.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"

using namespace pandora;

namespace lar_content
{

TopologyCalorimetryMonitoringAlgorithm::TopologyCalorimetryMonitoringAlgorithm() : 
    m_writeTree{true}
{
}

//------------------------------------------------------------------------------------------//

TopologyCalorimetryMonitoringAlgorithm::~TopologyCalorimetryMonitoringAlgorithm()
{
    if (m_writeTree)
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "RECREATE"));
    }
}

//------------------------------------------------------------------------------------------//

StatusCode TopologyCalorimetryMonitoringAlgorithm::Run()
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

            this->OutputAxis(pClusterList);
//            this->ClusterMerging(pClusterList, &clusterListName);
        }
        catch (StatusCodeException &statusCodeException)
        {
            throw statusCodeException;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//--------------------------------------------------------------------------------------------//

/* const MCParticle* TopologyCalorimetryMonitoringAlgorithm::GetMCForCluster(const Cluster *const cluster, std::map<const Cluster*,
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

//--------------------------------------------------------------------------------------------//

bool TopologyCalorimetryMonitoringAlgorithm::IsValidToUse(const Cluster *const cluster, std::map<const Cluster*, bool> &clusterIsUsed) const
{
    const int m_minNCaloHits{3};
    
    if (!cluster->IsAvailable())
        return false;

    if (cluster->GetNCaloHits() < m_minNCaloHits)
        return false;

    if (clusterIsUsed.count(cluster) > 0)
        return false;

    return true;
}
*/
//--------------------------------------------------------------------------------------------//

/* void TopologyCalorimetryMonitoringAlgorithm::ClusterMerging(const pandora::ClusterList *const pClusterList, const std::string &listName) const
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

*/
//--------------------------------------------------------------------------------------------//

void TopologyCalorimetryMonitoringAlgorithm::OutputAxis(const pandora::ClusterList *const pClusterList) const 
{
    for(const Cluster *const pCluster : *pClusterList)
    {
        CaloHitList clusterCaloHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(clusterCaloHitList);

        CartesianVector centroid(0.f, 0.f, 0.f);
        LArPcaHelper::EigenVectors eigenVecs;
        LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
        LArPcaHelper::RunPca(clusterCaloHitList, centroid, eigenValues, eigenVecs);
    
        const float clusterLength(LArClusterHelper::GetLength(pCluster));           
 
        if(!eigenVecs.empty())
        { 
            const float pcaX{eigenVecs.front().GetX()};
            const float pcaY{eigenVecs.front().GetY()};
            const float pcaZ{eigenVecs.front().GetZ()};
            const int   pdg{std::abs(MCParticleHelper::GetMainMCParticle(pCluster)->GetParticleId())};
            
            const double clusterAngle{std::tan( pcaZ / pcaX)};


            if(m_printToTerminal)
            {

                for(const CartesianVector &axis : eigenVecs)
                {
                    std::cout << axis  << std::endl;
                }
            }
            
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ClusterLength", clusterLength));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPDG", pdg));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PcaX", pcaX));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PcaY", pcaY));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PcaZ", pcaZ));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ClusterAngle", clusterAngle));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str())); 
       }
    }

    return;
}

//-------------------------------------------------------------------------------------------------

StatusCode TopologyCalorimetryMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputClusterListNames", m_inputClusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
 
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TreeName", m_treeName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "FileName", m_fileName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteTree", m_writeTree));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PrintToTerminal", m_printToTerminal));
 
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
