/**
 *  @file   larpandoracontent/LArCheating/CheatingCosmicRayTaggingTool.cc
 *
 *  @brief  Implementation of the cheating cosmic-ray tagging tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArCheating/CheatingCosmicRayTaggingTool.h"
#include "larpandoracontent/LArCheating/CheatingSliceIdBaseTool.h"

using namespace pandora;

namespace lar_content
{

CheatingCosmicRayTaggingTool::CheatingCosmicRayTaggingTool() :
    m_maxCosmicRayFraction(0.25f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingCosmicRayTaggingTool::FindAmbiguousPfos(const PfoList &parentCosmicRayPfos, PfoList &ambiguousPfos, const MasterAlgorithm *const /*pAlgorithm*/)
{
    if (this->GetPandora().GetSettings()->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    PfoList ambiguousParentPfos;

    for (const Pfo *const pParentCosmicRayPfo : parentCosmicRayPfos)
    {
        PfoList downstreamPfos;
        LArPfoHelper::GetAllDownstreamPfos(pParentCosmicRayPfo, downstreamPfos);

        float thisCosmicRayWeight(0.f), thisTotalWeight(0.f);
        CheatingSliceIdBaseTool::GetTargetParticleWeight(&downstreamPfos, thisCosmicRayWeight, thisTotalWeight, LArMCParticleHelper::IsCosmicRay);

        if ((thisTotalWeight > 0.f) && ((thisCosmicRayWeight / thisTotalWeight) < m_maxCosmicRayFraction))
            ambiguousParentPfos.push_back(pParentCosmicRayPfo);
    }

    LArPfoHelper::GetAllConnectedPfos(ambiguousParentPfos, ambiguousPfos);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingCosmicRayTaggingTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxCosmicRayFraction", m_maxCosmicRayFraction));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
