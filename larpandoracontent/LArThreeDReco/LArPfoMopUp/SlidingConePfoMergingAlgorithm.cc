/**
 *  @file   larpandoracontent/LArThreeDReco/LArPfoMopUp/SlidingConePfoMergingAlgorithm.cc
 * 
 *  @brief  Implementation of the sliding cone pfo merging algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingConeFitResult.h"

#include "larpandoracontent/LArThreeDReco/LArPfoMopUp/SlidingConePfoMergingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

SlidingConePfoMergingAlgorithm::SlidingConePfoMergingAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SlidingConePfoMergingAlgorithm::Run()
{
    ClusterToPfoMap clusterToPfoMap;
    ClusterVector trackClusters3D, showerClusters3D;
    this->GetThreeDClusters(trackClusters3D, showerClusters3D, clusterToPfoMap);

    

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SlidingConePfoMergingAlgorithm::GetThreeDClusters(ClusterVector &trackClusters3D, ClusterVector &showerClusters3D, ClusterToPfoMap &clusterToPfoMap) const
{
    for (const std::string &pfoListName : m_inputPfoListNames)
    {
        const PfoList *pPfoList(nullptr);

        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, pfoListName, pPfoList))
            continue;

        for (const Pfo *const pPfo : *pPfoList)
        {
            ClusterList pfoClusters3D;
            LArPfoHelper::GetThreeDClusterList(pPfo, pfoClusters3D);

            for (const Cluster *const pCluster3D : pfoClusters3D)
            {
                if (!clusterToPfoMap.insert(ClusterToPfoMap::value_type(pCluster3D, pPfo)).second)
                    throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);

                ClusterVector &targetVector((MU_MINUS == std::abs(pPfo->GetParticleId())) ? trackClusters3D : showerClusters3D);
                targetVector.push_back(pCluster3D);
            }
        }
    }

    std::sort(trackClusters3D.begin(), trackClusters3D.end(), LArClusterHelper::SortByNHits);
    std::sort(showerClusters3D.begin(), showerClusters3D.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SlidingConePfoMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "InputPfoListNames", m_inputPfoListNames));

    m_daughterListNames.insert(m_daughterListNames.end(), m_inputPfoListNames.begin(), m_inputPfoListNames.end());

    return PfoMergingBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
