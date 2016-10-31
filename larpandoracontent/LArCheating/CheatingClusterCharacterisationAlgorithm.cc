/**
 *  @file   larpandoracontent/LArCheating/CheatingClusterCharacterisationAlgorithm.cc
 * 
 *  @brief  Implementation of the cheating cluster characterisation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArCheating/CheatingClusterCharacterisationAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

using namespace pandora;

namespace lar_content
{

CheatingClusterCharacterisationAlgorithm::CheatingClusterCharacterisationAlgorithm() :
    m_overwriteExistingId(false),
    m_useUnavailableClusters(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingClusterCharacterisationAlgorithm::Run()
{
    for (const std::string &clusterListName : m_inputClusterListNames)
    {
        const ClusterList *pClusterList(nullptr);
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));

        if (!pClusterList || pClusterList->empty())
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "CheatingClusterCharacterisationAlgorithm: unable to find cluster list " << clusterListName << std::endl;

            continue;
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

bool CheatingClusterCharacterisationAlgorithm::IsClearTrack(const Cluster *const pCluster) const
{
    try
    {
        // ATTN Slightly curious definition of a clear track, but this is most-likely what is needed for shower-growing
        const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCluster));

        if ((PHOTON != pMCParticle->GetParticleId()) && (E_MINUS != std::abs(pMCParticle->GetParticleId())))
            return true;
    }
    catch (StatusCodeException &)
    {
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingClusterCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputClusterListNames", m_inputClusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OverwriteExistingId", m_overwriteExistingId));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseUnavailableClusters", m_useUnavailableClusters));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
