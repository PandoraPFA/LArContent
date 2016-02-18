/**
 *  @file   LArContent/src/LArStitching/StitchingPfoMergingTool.cc
 * 
 *  @brief  Implementation of the stitching pfo merging tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArPfoHelper.h"

#include "LArStitching/StitchingPfoMergingTool.h"

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
        primaryPfos.insert(LArPfoHelper::GetParentPfo(pPfo));

    if (2 != primaryPfos.size())
        return;

    const StitchingAlgorithm::PfoToVolumeIdMap &pfoToVolumeIdMap(stitchingInfo.m_pfoToVolumeIdMap);

    for (const ParticleFlowObject *const pPfo1 : primaryPfos)
    {
        for (const ParticleFlowObject *const pPfo2 : primaryPfos)
        {
            // ATTN Could ultimately check to see whether pfos are in adjacent volumes, not just different volumes
            if ((pPfo1 == pPfo2) || (pfoToVolumeIdMap.at(pPfo1) == pfoToVolumeIdMap.at(pPfo2)))
                continue;

            // ATTN Add key pattern-recognition logic here
            const bool shouldMergePfos(true);

            if (shouldMergePfos)
            {
                // ATTN Need to work out how to handle rebuilding of the two primary pfos into a single coherent pfo
                // My suggestion is to identify merges of individual primary particles, then update hierarchy information
                // Stitching of secondary particles is, of course, possible but needs either a specific treatment or a clever recursive approach

                

PfoList merge1, merge2;
merge1.insert(pPfo1);
merge2.insert(pPfo2);
PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &merge1, "merge1", RED));
PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &merge2, "merge2", GREEN));
PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
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
