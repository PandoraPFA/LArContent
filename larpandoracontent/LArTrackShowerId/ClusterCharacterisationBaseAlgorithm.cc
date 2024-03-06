/**
 *  @file   larpandoracontent/LArTrackShowerId/ClusterCharacterisationBaseAlgorithm.cc
 *
 *  @brief  Implementation of the cluster characterisation base algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArTrackShowerId/ClusterCharacterisationBaseAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ClusterCharacterisationBaseAlgorithm::ClusterCharacterisationBaseAlgorithm() :
    m_zeroMode(false), m_overwriteExistingId(false), m_useUnavailableClusters(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

ClusterCharacterisationBaseAlgorithm::~ClusterCharacterisationBaseAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterCharacterisationBaseAlgorithm::Run()
{
    for (const std::string &clusterListName : m_inputClusterListNames)
    {
        const ClusterList *pClusterList = NULL;
        PANDORA_RETURN_RESULT_IF_AND_IF(
            STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));

        if (!pClusterList || pClusterList->empty())
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "ClusterCharacterisationBaseAlgorithm: unable to find cluster list " << clusterListName << std::endl;

            continue;
        }

        if (m_zeroMode)
        {
            for (const Cluster *const pCluster : *pClusterList)
            {
                PandoraContentApi::Cluster::Metadata metadata;
                metadata.m_particleId = UNKNOWN_PARTICLE_TYPE;
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::AlterMetadata(*this, pCluster, metadata));
            }

            return STATUS_CODE_SUCCESS;
        }

        for (const Cluster *const pCluster : *pClusterList)
        {
            if (!m_overwriteExistingId && (UNKNOWN_PARTICLE_TYPE != pCluster->GetParticleId()))
                continue;

            if (!m_useUnavailableClusters && !PandoraContentApi::IsAvailable(*this, pCluster))
                continue;

            PandoraContentApi::Cluster::Metadata metadata;

            if (this->IsClearTrack(pCluster))
            {
                metadata.m_particleId = MU_MINUS;
            }
            else
            {
                metadata.m_particleId = E_MINUS;
            }

            if (pCluster->GetParticleId() != metadata.m_particleId.Get())
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::AlterMetadata(*this, pCluster, metadata));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterCharacterisationBaseAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputClusterListNames", m_inputClusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ZeroMode", m_zeroMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "OverwriteExistingId", m_overwriteExistingId));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UseUnavailableClusters", m_useUnavailableClusters));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
