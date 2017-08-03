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

CheatingCosmicRayTaggingTool::CheatingCosmicRayTaggingTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingCosmicRayTaggingTool::FindAmbiguousPfos(const PfoList &parentCosmicRayPfos, PfoList &ambiguousPfos) const
{
    PfoList ambiguousParentPfos;                                                                                                            
                                                                                                                                            
    for (const Pfo *const pParentCosmicRayPfo : parentCosmicRayPfos)                                                                        
    {
        if (LArPfoHelper::IsNeutrino(pParentCosmicRayPfo))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        const float neutrinoScore(LArMCParticleHelper::GetDownstreamNeutrinoScore(pParentCosmicRayPfo));

        // TODO Sort out normalisation of this score
        PfoList downstreamPfos;
        LArPfoHelper::GetAllDownstreamPfos(pParentCosmicRayPfo, downstreamPfos);
        downstreamPfos.sort(LArPfoHelper::SortByNHits);

        float nHits(0.f);

        for (const Pfo *const pDownstreamPfo : downstreamPfos)
        {
            ClusterList twoDClusters;
            LArPfoHelper::GetTwoDClusterList(pDownstreamPfo, twoDClusters);
            twoDClusters.sort(LArClusterHelper::SortByNHits);

            for (const Cluster *const pCluster : twoDClusters)
                nHits += static_cast<float>(pCluster->GetNCaloHits());
        }

        if ((nHits > std::numeric_limits<float>::epsilon()) && ((neutrinoScore / nHits) > 0.75f)) // TODO
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
