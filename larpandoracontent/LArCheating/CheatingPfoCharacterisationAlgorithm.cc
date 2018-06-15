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

bool CheatingPfoCharacterisationAlgorithm::IsClearTrack(const ParticleFlowObject *const pPfo)
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

    return ((PHOTON != pBestMCParticle->GetParticleId()) && (E_MINUS != std::abs(pBestMCParticle->GetParticleId())));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CheatingPfoCharacterisationAlgorithm::IsClearTrack(const Cluster *const /*pCluster*/) const
{
    throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);
}

} // namespace lar_content
