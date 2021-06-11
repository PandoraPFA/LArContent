/**
 *  @file   larpandoracontent/LArUtility/PfoHitCleaningAlgorithm.cc
 *
 *  @brief  Implementation of the pfo hit cleaning algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArUtility/PfoHitCleaningAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode PfoHitCleaningAlgorithm::Run()
{
    for (unsigned int i = 0; i < m_pfoListNames.size(); ++i)
    {
        // ATTN - one-to-one correspondance required between PFO and cluster lists
        const std::string &pfoListName{m_pfoListNames.at(i)};
        const std::string &clusterListName{m_clusterListNames.at(i)};
        const PfoList *pList(nullptr);
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, pfoListName, pList));

        if (pList && !pList->empty())
        {
            for (const ParticleFlowObject *pPfo : *pList)
            {
                ClusterList clustersToRemove;
                LArPfoHelper::GetClusters(pPfo, TPC_3D, clustersToRemove);

                for (const Cluster *pCluster : clustersToRemove)
                {
                    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromPfo(*this, pPfo, pCluster));
                    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                        PandoraContentApi::Delete<Cluster>(*this, pCluster, clusterListName));
                }
            }
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoHitCleaningAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "PfoListNames", m_pfoListNames));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "ClusterListNames", m_clusterListNames));

    if (m_pfoListNames.size() != m_clusterListNames.size())
    {
        std::cout << "PfoHitCleaningAlgorithm: Mismatch between PFO and Cluster list sizes" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
