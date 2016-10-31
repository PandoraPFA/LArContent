/**
 *  @file   larpandoracontent/LArCheating/CheatingPfoCharacterisationAlgorithm.cc
 * 
 *  @brief  Implementation of the cheating pfo characterisation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArCheating/CheatingPfoCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CheatingPfoCharacterisationAlgorithm::CheatingPfoCharacterisationAlgorithm() :
    m_updateClusterIds(true)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingPfoCharacterisationAlgorithm::Run()
{
    PfoList tracksToShowers, showersToTracks;

    for (const std::string &pfoListName : m_inputPfoListNames)
    {
        const PfoList *pPfoList(nullptr);
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, pfoListName, pPfoList));

        if (!pPfoList || pPfoList->empty())
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "CheatingPfoCharacterisationAlgorithm: unable to find pfo list " << pfoListName << std::endl;

            continue;
        }

        for (const ParticleFlowObject *const pPfo : *pPfoList)
        {
            PandoraContentApi::ParticleFlowObject::Metadata pfoMetadata;

            if (this->IsClearTrack(pPfo))
            {
                pfoMetadata.m_particleId = MU_MINUS;

                if (m_showerPfoListName == pfoListName)
                    showersToTracks.push_back(pPfo);
            }
            else
            {
                pfoMetadata.m_particleId = E_MINUS;

                if (m_trackPfoListName == pfoListName)
                    tracksToShowers.push_back(pPfo);
            }

            if (pPfo->GetParticleId() != pfoMetadata.m_particleId.Get())
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, pfoMetadata));

            if (!m_updateClusterIds)
                continue;

            ClusterList twoDClusterList;
            LArPfoHelper::GetTwoDClusterList(pPfo, twoDClusterList);

            for (const Cluster *const pCluster : twoDClusterList)
            {
                if (pCluster->GetParticleId() == pfoMetadata.m_particleId.Get())
                    continue;

                PandoraContentApi::Cluster::Metadata clusterMetadata;
                clusterMetadata.m_particleId = pfoMetadata.m_particleId.Get();
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::AlterMetadata(*this, pCluster, clusterMetadata));
            }
        }
    }

    if (!tracksToShowers.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_trackPfoListName, m_showerPfoListName, tracksToShowers));

    if (!showersToTracks.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_showerPfoListName, m_trackPfoListName, showersToTracks));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CheatingPfoCharacterisationAlgorithm::IsClearTrack(const ParticleFlowObject *const pPfo) const
{
    CaloHitList caloHitList;
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_U, caloHitList);
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_V, caloHitList);
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, caloHitList);

    MCParticleWeightMap mcParticleWeightMap;

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        for (const MCParticleWeightMap::value_type &mapEntry : pCaloHit->GetMCParticleWeightMap())
            mcParticleWeightMap[mapEntry.first] += mapEntry.second;
    }

    float bestWeight(0.f);
    const MCParticle *pBestMCParticle(nullptr);

    MCParticleList mcParticleList;
    for (const auto &mapEntry : mcParticleWeightMap) mcParticleList.push_back(mapEntry.first);
    mcParticleList.sort(LArMCParticleHelper::SortByMomentum);

    for (const MCParticle *const pMCParticle : mcParticleList)
    {
        const float weight(mcParticleWeightMap.at(pMCParticle));

        if (weight > bestWeight)
        {
            pBestMCParticle = pMCParticle;
            bestWeight = weight;
        }
    }

    if (!pBestMCParticle)
        return false;

    const MCParticle *const pPrimaryMCParticle(LArMCParticleHelper::GetPrimaryMCParticle(pBestMCParticle));
    const int absParticleId(std::abs(pPrimaryMCParticle->GetParticleId()));
    return ((MU_MINUS == absParticleId) || (PROTON == absParticleId) || (PI_PLUS == absParticleId));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingPfoCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrackPfoListName", m_trackPfoListName));
    m_inputPfoListNames.push_back(m_trackPfoListName);

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ShowerPfoListName", m_showerPfoListName));
    m_inputPfoListNames.push_back(m_showerPfoListName);

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UpdateClusterIds", m_updateClusterIds));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
