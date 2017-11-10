/**
 *  @file   larpandoracontent/LArThreeDReco/LArEventBuilding/NeutrinoIdTool.cc
 *
 *  @brief  Implementation of the neutrino id tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include "larpandoracontent/LArThreeDReco/LArEventBuilding/NeutrinoIdTool.h"

using namespace pandora;

namespace lar_content
{

NeutrinoIdTool::NeutrinoIdTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoIdTool::FillNeutrinoProperties(const PfoList *const pPfoList, SliceProperties &sliceProperties) const
{
    if (1 != pPfoList->size())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const ParticleFlowObject *const pNeutrinoPfo(pPfoList->front());

    if (!LArPfoHelper::IsNeutrino(pNeutrinoPfo))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const CartesianVector vertexPosition(LArPfoHelper::GetVertex(pNeutrinoPfo)->GetPosition());
    sliceProperties.m_nuVtxY = vertexPosition.GetY();
    sliceProperties.m_nuVtxZ = vertexPosition.GetZ();

    CartesianVector weightedDirection(0.f, 0.f, 0.f);

    for (const ParticleFlowObject *const pDaughterPfo : pNeutrinoPfo->GetDaughterPfoList())
    {
        ClusterList clusters3D;
        LArPfoHelper::GetThreeDClusterList(pDaughterPfo, clusters3D);

        if (1 != clusters3D.size())
            continue;

        const ThreeDSlidingFitResult slidingFitResult(clusters3D.front(), 100, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const CartesianVector endMin(slidingFitResult.GetGlobalMinLayerPosition());
        const CartesianVector endMax(slidingFitResult.GetGlobalMaxLayerPosition());
        const CartesianVector dirMin(slidingFitResult.GetGlobalMinLayerDirection().GetUnitVector());
        const CartesianVector dirMax(slidingFitResult.GetGlobalMaxLayerDirection().GetUnitVector());

        const bool startAtMin((endMin - vertexPosition).GetMagnitudeSquared() < (endMax - vertexPosition).GetMagnitudeSquared());
        const CartesianVector start(startAtMin ? endMin : endMax);
        const CartesianVector end(startAtMin ? endMax : endMin);
        const CartesianVector directionAtVertex(startAtMin ? dirMin : dirMax);

        const float directionCorrection(((end - start).GetDotProduct(directionAtVertex) > 0.f) ? 1.f : -1.f);
        weightedDirection += directionAtVertex * directionCorrection * static_cast<float>(clusters3D.front()->GetNCaloHits());
    }

    sliceProperties.m_nuCosWeightedDir = weightedDirection.GetUnitVector().GetZ();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoIdTool::FillCosmicRayProperties(const PfoList *const /*pPfoList*/, SliceProperties &/*sliceProperties*/) const
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool NeutrinoIdTool::GetNeutrinoSliceIndex(const SliceIndexToPropertiesMap &sliceIndexToPropertiesMap, unsigned int &neutrinoSliceIndex) const
{
    float maxLikelihood(-std::numeric_limits<float>::max());

    for (const auto &mapEntry : sliceIndexToPropertiesMap)
    {
        const SliceProperties &sliceProperties(mapEntry.second);

        // TODO This is just a preliminary parameterisation, to be developed and all parameters made configurable
        const float likelihood_nuVtxY((sliceProperties.m_nuVtxY < 93.f) && (sliceProperties.m_nuVtxY > -93.f) ? 0.6f : 0.1f);
        const float likelihood_nuVtxZ((sliceProperties.m_nuVtxY > 10.4f) ? 0.65f : 0.1f);
        const float likelihood_nuCosWeightedDir((sliceProperties.m_nuCosWeightedDir > 0.f) ?
            (0.436f * std::exp(sliceProperties.m_nuCosWeightedDir ) - 0.236f) :
            (0.2f - 0.4f * sliceProperties.m_nuCosWeightedDir));

        const float likelihood(likelihood_nuVtxY * likelihood_nuVtxZ * likelihood_nuCosWeightedDir);

        if (likelihood > maxLikelihood)
        {
            maxLikelihood = likelihood;
            neutrinoSliceIndex = mapEntry.first;
        }
    }

    return (!sliceIndexToPropertiesMap.empty());
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode NeutrinoIdTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
