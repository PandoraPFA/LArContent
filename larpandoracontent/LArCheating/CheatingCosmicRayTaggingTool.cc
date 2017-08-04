/**
 *  @file   larpandoracontent/LArCheating/CheatingCosmicRayTaggingTool.cc
 *
 *  @brief  Implementation of the cheating cosmic-ray tagging tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

    #include "larpandoracontent/LArHelpers/LArClusterHelper.h" // TODO
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArCheating/CheatingCosmicRayTaggingTool.h"

using namespace pandora;

namespace lar_content
{

CheatingCosmicRayTaggingTool::CheatingCosmicRayTaggingTool() :
    m_minNeutrinoFraction(0.75f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingCosmicRayTaggingTool::FindAmbiguousPfos(const PfoList &parentCosmicRayPfos, PfoList &ambiguousPfos) const
{
    PfoList ambiguousParentPfos;                                                                                                            
                                                                                                                                            
    for (const Pfo *const pParentCosmicRayPfo : parentCosmicRayPfos)                                                                        
    {
        PfoList downstreamPfos;
        LArPfoHelper::GetAllDownstreamPfos(pParentCosmicRayPfo, downstreamPfos);
        const float neutrinoFraction(LArMCParticleHelper::GetNeutrinoFraction(&downstreamPfos));

        if (neutrinoFraction > m_minNeutrinoFraction)
            ambiguousParentPfos.push_back(pParentCosmicRayPfo);
    }

    LArPfoHelper::GetAllConnectedPfos(ambiguousParentPfos, ambiguousPfos);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingCosmicRayTaggingTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinNeutrinoFraction", m_minNeutrinoFraction));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
