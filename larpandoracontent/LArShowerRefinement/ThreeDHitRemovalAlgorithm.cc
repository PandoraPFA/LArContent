/**
 *  @file   larpandoracontent/LArShowerRefinement/ThreeDHitRemovalAlgorithm.cc
 *
 *  @brief  Implementation of the 3D hit removal algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArShowerRefinement/ThreeDHitRemovalAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ThreeDHitRemovalAlgorithm::ThreeDHitRemovalAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDHitRemovalAlgorithm::Run()
{
    for (const std::string &pfoListName : m_pfoListNames)
    {
        const PfoList *pPfoList(nullptr);
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, pfoListName, pPfoList));

        if (!pPfoList || pPfoList->empty())
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "ThreeDHitRemovalAlgorithm: unable to find pfo list " << pfoListName << std::endl;

            continue;
        }

        for (const ParticleFlowObject *const pPfo : *pPfoList)
            this->RemoveSpacePoints(pPfo);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDHitRemovalAlgorithm::RemoveSpacePoints(const ParticleFlowObject *const pPfo) const
{   
    ClusterList threeDClusterList;
    LArPfoHelper::GetClusters(pPfo, TPC_3D, threeDClusterList);

    for (const Cluster *const pCluster : threeDClusterList)
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromPfo(*this, pPfo, pCluster));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDHitRemovalAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadVectorOfValues(xmlHandle, "PfoListNames", m_pfoListNames));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
