/**
 *  @file   larpandoracontent/LArCheating/CosmicRayShowerMatchingAlg.cc
 *
 *  @brief  Implementation of the cheater for the cosmic ray shower matching algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArCheating/CheatingCosmicRayShowerMatchingAlg.h"

using namespace pandora;

namespace lar_content
{

StatusCode CheatingCosmicRayShowerMatchingAlg::Run()
{
    ClusterList candidateClusterList;
    this->GetCandidateClusters(candidateClusterList);

    const PfoList *pPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputPfoListName, pPfoList));

    for (const ParticleFlowObject *const pPfo : *pPfoList)
    {
        ClusterList twoDClusters;
        LArPfoHelper::GetTwoDClusterList(pPfo, twoDClusters);

        for (const Cluster *const pPfoCluster : twoDClusters)
            this->CosmicRayShowerMatching(pPfo, pPfoCluster, candidateClusterList);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingCosmicRayShowerMatchingAlg::GetCandidateClusters(ClusterList &candidateClusterList) const
{
    for (const std::string &clusterListName : m_inputClusterListNames)
    {
        const ClusterList *pClusterList(nullptr);

        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, clusterListName, pClusterList))
        {
            std::cout << "CheatingCosmicRayShowerMatchingAlg - Could not access cluster list with name " << clusterListName << std::endl;
            continue;
        }

        candidateClusterList.insert(candidateClusterList.end(), pClusterList->begin(), pClusterList->end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingCosmicRayShowerMatchingAlg::CosmicRayShowerMatching(
    const ParticleFlowObject *const pPfo, const Cluster *const pPfoCluster, const ClusterList &candidateClusterList) const
{
    try
    {
        const HitType pfoClusterHitType(LArClusterHelper::GetClusterHitType(pPfoCluster));
        const MCParticle *const pPfoMCParticle(MCParticleHelper::GetMainMCParticle(pPfoCluster));
        const MCParticle *const pPfoParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pPfoMCParticle));

        for (const Cluster *const pCandidateCluster : candidateClusterList)
        {
            if (pfoClusterHitType != LArClusterHelper::GetClusterHitType(pCandidateCluster))
                continue;

            if (!PandoraContentApi::IsAvailable(*this, pCandidateCluster))
                continue;

            if (pPfoCluster == pCandidateCluster)
                continue;

            try
            {
                const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCandidateCluster));
                const MCParticle *const pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pMCParticle));

                if (!LArMCParticleHelper::IsNeutrino(pParentMCParticle) && (pPfoParentMCParticle == pParentMCParticle))
                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pPfo, pCandidateCluster));
            }
            catch (const StatusCodeException &)
            {
            }
        }
    }
    catch (const StatusCodeException &)
    {
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingCosmicRayShowerMatchingAlg::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputPfoListName", m_inputPfoListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputClusterListNames", m_inputClusterListNames));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
