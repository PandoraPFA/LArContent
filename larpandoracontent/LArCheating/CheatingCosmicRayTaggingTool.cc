/**
 *  @file   larpandoracontent/LArCheating/CheatingCosmicRayTaggingTool.cc
 *
 *  @brief  Implementation of the cheating cosmic-ray tagging tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArCheating/CheatingCosmicRayTaggingTool.h"

using namespace pandora;

namespace lar_content
{

CheatingCosmicRayTaggingTool::CheatingCosmicRayTaggingTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingCosmicRayTaggingTool::FindAmbiguousPfos(const PfoList &parentCosmicRayPfos, PfoList &ambiguousPfos) const
{
    // TODO
    int counter(0);
    PfoList ambiguousParentPfos;                                                                                                            
                                                                                                                                            
    for (const Pfo *const pParentCosmicRayPfo : parentCosmicRayPfos)                                                                        
    {                                                                                                                                       
        if (++counter % 4 == 0)                                                                                                             
            ambiguousParentPfos.push_back(pParentCosmicRayPfo);                                                                             
    }                                                                                                                                       

    LArPfoHelper::GetAllConnectedPfos(ambiguousParentPfos, ambiguousPfos);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingCosmicRayTaggingTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
