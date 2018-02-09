/**
 *  @file   larpandoracontent/LArCheating/CheatingNeutrinoIdTool.cc
 *
 *  @brief  Implementation of the cheating neutrino id tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArCheating/CheatingNeutrinoIdTool.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

using namespace pandora;

namespace lar_content
{

void CheatingNeutrinoIdTool::SelectOutputPfos(const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, PfoList &selectedPfos)
{
    if (nuSliceHypotheses.size() != crSliceHypotheses.size())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    float bestNeutrinoWeight(0.f);
    unsigned int bestSliceIndex(std::numeric_limits<unsigned int>::max());

    for (unsigned int sliceIndex = 0, nSlices = nuSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
    {
        float neutrinoWeight(0.f);
        const PfoList &neutrinoPfoList(nuSliceHypotheses.at(sliceIndex));

        for (const Pfo *const pNeutrinoPfo : neutrinoPfoList)
        {
            if (!LArPfoHelper::IsNeutrino(pNeutrinoPfo))
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

            PfoList downstreamPfos;
            LArPfoHelper::GetAllDownstreamPfos(pNeutrinoPfo, downstreamPfos);

            float thisNeutrinoWeight(0.f), thisTotalWeight(0.f);
            this->GetNeutrinoWeight(&downstreamPfos, thisNeutrinoWeight, thisTotalWeight);
            neutrinoWeight += thisNeutrinoWeight;
        }

        if (neutrinoWeight > bestNeutrinoWeight)
        {
            bestNeutrinoWeight = neutrinoWeight;
            bestSliceIndex = sliceIndex;
        }
    }

    for (unsigned int sliceIndex = 0, nSlices = nuSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
    {
        const PfoList &sliceOutput((bestSliceIndex == sliceIndex) ? nuSliceHypotheses.at(sliceIndex) : crSliceHypotheses.at(sliceIndex));
        selectedPfos.insert(selectedPfos.end(), sliceOutput.begin(), sliceOutput.end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingNeutrinoIdTool::GetNeutrinoWeight(const CaloHit *const pCaloHit, float &neutrinoWeight, float &totalWeight) const
{
    neutrinoWeight = 0.f; totalWeight = 0.f;

    const CaloHit *const pCaloHitMaster(reinterpret_cast<const CaloHit *>(pCaloHit->GetParentAddress()));

    if (!pCaloHitMaster)
    {
        std::cout << "CheatingNeutrinoIdTool: failed to cast CaloHit parent address as CaloHit" << std::endl;
        return;
    }

    const MCParticleWeightMap &hitMCParticleWeightMap(pCaloHitMaster->GetMCParticleWeightMap());

    if (hitMCParticleWeightMap.empty())
        return;

    MCParticleList mcParticleList;
    for (const auto &mapEntry : hitMCParticleWeightMap) mcParticleList.push_back(mapEntry.first);
    mcParticleList.sort(LArMCParticleHelper::SortByMomentum);

    for (const MCParticle *const pMCParticle : mcParticleList)
    {
        const float weight(hitMCParticleWeightMap.at(pMCParticle));

        if (LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCParticle))
            neutrinoWeight += weight;

        totalWeight += weight;
    }

    // ATTN normalise arbitrary input weights at this point
    if (totalWeight > std::numeric_limits<float>::epsilon())
    {
        neutrinoWeight *= 1.f / totalWeight;
        totalWeight = 1.f;
    }
    else
    {
        neutrinoWeight = 0.f;
        totalWeight = 0.f;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingNeutrinoIdTool::GetNeutrinoWeight(const PfoList *const pPfoList, float &neutrinoWeight, float &totalWeight) const
{
    neutrinoWeight = 0.f; totalWeight = 0.f;

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
            float thisNeutrinoWeight = 0.f, thisTotalWeight = 0.f;
            this->GetNeutrinoWeight(pCaloHit, thisNeutrinoWeight, thisTotalWeight);

            neutrinoWeight += thisNeutrinoWeight;
            totalWeight += thisTotalWeight;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingNeutrinoIdTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
