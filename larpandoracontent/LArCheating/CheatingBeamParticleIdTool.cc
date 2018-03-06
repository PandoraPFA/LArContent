/**
 *  @file   larpandoracontent/LArCheating/CheatingBeamParticleIdTool.cc
 *
 *  @brief  Implementation of the cheating beam particle id tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArCheating/CheatingBeamParticleIdTool.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

using namespace pandora;

namespace lar_content
{

CheatingBeamParticleIdTool::CheatingBeamParticleIdTool() : 
    m_minWeightFraction(0.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingBeamParticleIdTool::GetBeamParticleWeight(const PfoList *const pPfoList, const bool objectOwnedByMaster, float &beamParticleWeight, float &totalWeight)
{
    beamParticleWeight = 0.f; totalWeight = 0.f;

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
            float thisBeamParticleWeight = 0.f, thisTotalWeight = 0.f;
            CheatingBeamParticleIdTool::GetBeamParticleWeight(pCaloHit, objectOwnedByMaster, thisBeamParticleWeight, thisTotalWeight);

            beamParticleWeight += thisBeamParticleWeight;
            totalWeight += thisTotalWeight;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingBeamParticleIdTool::GetBeamParticleWeight(const CaloHit *const pCaloHit, const bool objectOwnedByMaster, float &beamParticleWeight, float &totalWeight)
{
    beamParticleWeight = 0.f; totalWeight = 0.f;

    const CaloHit *const pCaloHitMaster(objectOwnedByMaster ? pCaloHit : static_cast<const CaloHit *>(pCaloHit->GetParentAddress()));
    const MCParticleWeightMap &hitMCParticleWeightMap(pCaloHitMaster->GetMCParticleWeightMap());

    if (hitMCParticleWeightMap.empty())
        return;

    MCParticleList mcParticleList;
    for (const auto &mapEntry : hitMCParticleWeightMap) mcParticleList.push_back(mapEntry.first);
    mcParticleList.sort(LArMCParticleHelper::SortByMomentum);

    for (const MCParticle *const pMCParticle : mcParticleList)
    {
        const float weight(hitMCParticleWeightMap.at(pMCParticle));
        const MCParticle *const pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pMCParticle));

        if (LArMCParticleHelper::IsBeamParticle(pParentMCParticle))
            beamParticleWeight += weight;

        totalWeight += weight;
    }

    // ATTN normalise arbitrary input weights at this point
    if (totalWeight > std::numeric_limits<float>::epsilon())
    {
        beamParticleWeight *= 1.f / totalWeight;
        totalWeight = 1.f;
    }
    else
    {
        beamParticleWeight = 0.f;
        totalWeight = 0.f;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingBeamParticleIdTool::SelectOutputPfos(const pandora::Algorithm *const /*pAlgorithm*/, const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, PfoList &selectedPfos)
{
    if (nuSliceHypotheses.size() != crSliceHypotheses.size())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    IntVector beamSliceIndices;

    for (unsigned int sliceIndex = 0, nSlices = nuSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
    {
        float beamParticleWeight(0.f), totalWeight(0.f);
        const PfoList &neutrinoPfoList(nuSliceHypotheses.at(sliceIndex));

        for (const Pfo *const pNeutrinoPfo : neutrinoPfoList)
        {
            if (!LArPfoHelper::IsNeutrino(pNeutrinoPfo))
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

            PfoList downstreamPfos;
            LArPfoHelper::GetAllDownstreamPfos(pNeutrinoPfo, downstreamPfos);

            float thisBeamParticleWeight(0.f), thisTotalWeight(0.f);
            CheatingBeamParticleIdTool::GetBeamParticleWeight(&downstreamPfos, false, thisBeamParticleWeight, thisTotalWeight);

            beamParticleWeight += thisBeamParticleWeight;
            totalWeight += thisTotalWeight;
        }

        if (beamParticleWeight/totalWeight > m_minWeightFraction)
        {
            beamSliceIndices.push_back(sliceIndex);
        }
    }

    for (unsigned int sliceIndex = 0, nSlices = nuSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
    {
        bool isBeam(false);

        if (std::find(beamSliceIndices.begin(), beamSliceIndices.end(), sliceIndex) != beamSliceIndices.end())
            isBeam = true;

        const PfoList &sliceOutput(isBeam ? nuSliceHypotheses.at(sliceIndex) : crSliceHypotheses.at(sliceIndex));
        selectedPfos.insert(selectedPfos.end(), sliceOutput.begin(), sliceOutput.end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingBeamParticleIdTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MinimumWeightFraction", m_minWeightFraction));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
