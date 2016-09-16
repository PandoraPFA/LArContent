/**
 *  @file   larpandoracontent/LArStitching/StitchingPfoMergingTool.cc
 * 
 *  @brief  Implementation of the stitching pfo merging tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArStitching/StitchingPfoMergingTool.h"

using namespace pandora;

namespace lar_content
{

void StitchingPfoMergingTool::Run(const StitchingAlgorithm *const pAlgorithm, StitchingAlgorithm::StitchingInfo &stitchingInfo)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << this->GetType() << std::endl;

    const PandoraInstanceList &pandoraInstances(MultiPandoraApi::GetDaughterPandoraInstanceList(&(this->GetPandora())));

    if (2 != pandoraInstances.size())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    // ATTN Handy display for initial debug purposes
    const PfoList *pPfoList1(nullptr), *pPfoList2(nullptr), *pRecreatedPfoList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*pandoraInstances.at(0), pPfoList1));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*pandoraInstances.at(1), pPfoList2));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*pAlgorithm, pRecreatedPfoList));

    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
    PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), pPfoList1, "InputPfoList1", BLUE));
    PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), pPfoList2, "InputPfoList2", GREEN));
    PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), pRecreatedPfoList, "RecreatedPfoList", RED));
    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    // ATTN Example, place-holder implementation
    const PfoList *pNewPfoList = nullptr;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*pAlgorithm, pNewPfoList));

    PfoList primaryPfos;

    for (const ParticleFlowObject *const pPfo : *pNewPfoList)
        primaryPfos.push_back(LArPfoHelper::GetParentPfo(pPfo));

    if (2 != primaryPfos.size())
        return;

    // ATTN Can add other information to StitchingInfo (prototype in StitchingAlgorithm.h; values are stored in StitchingObjectCreationTool)
    StitchingAlgorithm::PfoToVolumeIdMap &pfoToVolumeIdMap(stitchingInfo.m_pfoToVolumeIdMap);

    for (const ParticleFlowObject *const pPfo1 : primaryPfos)
    {
        for (const ParticleFlowObject *const pPfo2 : primaryPfos)
        {
            // ATTN Could ultimately check to see whether pfos are in adjacent volumes, not just different volumes
            if ((pPfo1 == pPfo2) || (pfoToVolumeIdMap.at(pPfo1) == pfoToVolumeIdMap.at(pPfo2)))
                continue;

            // ATTN Add key pattern-recognition decision #1 here. Write a function to return this boolean.
            const bool shouldMergePfos(true);

            if (shouldMergePfos)
            {
                // ATTN Add key pattern-recognition decision #2 here. Write a function to return this boolean.
                const ParticleFlowObject *const pPfoToEnlarge(pPfo1);
                const ParticleFlowObject *const pPfoToDelete(pPfo2);

                // ATTN Update stitching information here, just placeholder/example update shown below
                pfoToVolumeIdMap[pPfoToEnlarge] = -1;
                pfoToVolumeIdMap.erase(pPfoToDelete);

                const PfoList enlargePfoList(1, pPfoToEnlarge), deletePfoList(1, pPfoToDelete);
                PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &enlargePfoList, "enlargePfoList", BLUE));
                PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &deletePfoList, "deletePfoList", GREEN));
                PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

                pAlgorithm->MergeAndDeletePfos(pPfoToEnlarge, pPfoToDelete);

                const PfoList mergedPfoList(1, pPfoToEnlarge);
                PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &mergedPfoList, "mergedPfoList", RED));
                PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

                // ATTN Some decisions could be made in very different ways, e.g. which vertices to keep, hierarchy, etc. This is all just a suggestion.

                // Now need to break-out from this example approach as we have invalidated current iterators. Care with subsequent implementation.
                return;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode StitchingPfoMergingTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
