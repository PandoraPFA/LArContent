/**
 *  @file   larpandoracontent/LArCheating/CheatingSliceIdBaseTool.cc
 *
 *  @brief  Implementation of the cheating slice id base tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArCheating/CheatingSliceIdBaseTool.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

using namespace pandora;

namespace lar_content
{

void CheatingSliceIdBaseTool::GetTargetParticleWeight(
    const PfoList *const pPfoList, float &targetParticleWeight, float &totalWeight, std::function<bool(const MCParticle *const)> fCriteria)
{
    targetParticleWeight = 0.f;
    totalWeight = 0.f;

    for (const ParticleFlowObject *const pPfo : *pPfoList)
    {
        ClusterList twoDClusters;
        LArPfoHelper::GetTwoDClusterList(pPfo, twoDClusters);

        CaloHitList caloHitList;

        for (const Cluster *const pCluster : twoDClusters)
        {
            const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));

            if ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType))
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

            pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);
        }

        for (const CaloHit *const pCaloHit : caloHitList)
        {
            float thisTargetParticleWeight = 0.f, thisTotalWeight = 0.f;
            CheatingSliceIdBaseTool::GetTargetParticleWeight(pCaloHit, thisTargetParticleWeight, thisTotalWeight, fCriteria);

            targetParticleWeight += thisTargetParticleWeight;
            totalWeight += thisTotalWeight;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingSliceIdBaseTool::GetTargetParticleWeight(
    const CaloHit *const pCaloHit, float &targetParticleWeight, float &totalWeight, std::function<bool(const MCParticle *const)> fCriteria)
{
    targetParticleWeight = 0.f;
    totalWeight = 0.f;
    const MCParticleWeightMap &hitMCParticleWeightMap(pCaloHit->GetMCParticleWeightMap());

    if (hitMCParticleWeightMap.empty())
        return;

    MCParticleList mcParticleList;
    for (const auto &mapEntry : hitMCParticleWeightMap)
        mcParticleList.push_back(mapEntry.first);
    mcParticleList.sort(LArMCParticleHelper::SortByMomentum);

    for (const MCParticle *const pMCParticle : mcParticleList)
    {
        const float weight(hitMCParticleWeightMap.at(pMCParticle));
        const MCParticle *const pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pMCParticle));

        if (fCriteria(pParentMCParticle))
            targetParticleWeight += weight;

        totalWeight += weight;
    }

    // ATTN normalise arbitrary input weights at this point
    if (totalWeight > std::numeric_limits<float>::epsilon())
    {
        targetParticleWeight *= 1.f / totalWeight;
        totalWeight = 1.f;
    }
    else
    {
        targetParticleWeight = 0.f;
        totalWeight = 0.f;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingSliceIdBaseTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
