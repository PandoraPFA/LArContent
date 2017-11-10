/**
 *  @file   larpandoracontent/LArCheating/CheatingNeutrinoIdTool.cc
 *
 *  @brief  Implementation of the cheating neutrino id tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArCheating/CheatingNeutrinoIdTool.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

using namespace pandora;

namespace lar_content
{

CheatingNeutrinoIdTool::CheatingNeutrinoIdTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingNeutrinoIdTool::FillNeutrinoProperties(const PfoList *const pPfoList, SliceProperties &sliceProperties) const
{
    float neutrinoWeight(0.f);

    for (const Pfo *const pNeutrinoPfo : *pPfoList)
    {
        if (!LArPfoHelper::IsNeutrino(pNeutrinoPfo))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        PfoList downstreamPfos;
        LArPfoHelper::GetAllDownstreamPfos(pNeutrinoPfo, downstreamPfos);

        float thisNeutrinoWeight(0.f), thisTotalWeight(0.f);
        LArMCParticleHelper::GetNeutrinoWeight(&downstreamPfos, thisNeutrinoWeight, thisTotalWeight);
        neutrinoWeight += thisNeutrinoWeight;
    }

    sliceProperties.m_weight = neutrinoWeight;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingNeutrinoIdTool::FillCosmicRayProperties(const PfoList *const /*pPfoList*/, SliceProperties &/*sliceProperties*/) const
{
    // Deliberately empty
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CheatingNeutrinoIdTool::GetNeutrinoSliceIndex(const SliceIndexToPropertiesMap &sliceIndexToPropertiesMap, unsigned int &neutrinoSliceIndex) const
{
    float bestWeight(0.f);
    unsigned int bestSliceIndex(std::numeric_limits<unsigned int>::max());

    for (const auto &mapEntry : sliceIndexToPropertiesMap)
    {
        if (mapEntry.second.m_weight > bestWeight)
        {
            bestWeight = mapEntry.second.m_weight;
            bestSliceIndex = mapEntry.first;
        }
    }

    if (bestWeight < std::numeric_limits<float>::epsilon())
        return false;

    neutrinoSliceIndex = bestSliceIndex;
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingNeutrinoIdTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
