/**
 *  @file   LArContent/src/LArStitching/StitchingPfoMergingTool.cc
 * 
 *  @brief  Implementation of the stitching pfo merging tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArStitching/StitchingPfoMergingTool.h"

using namespace pandora;

namespace lar_content
{

void StitchingPfoMergingTool::Run(const StitchingAlgorithm *const pAlgorithm, const MultiPandora &multiPandora,
    StitchingAlgorithm::StitchingInfo &/*stitchingInfo*/)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << this->GetType() << std::endl;

    const MultiPandora::PandoraInstanceList &pandoraInstances(multiPandora.GetDaughterPandoraInstances(&(this->GetPandora())));

    if (2 != pandoraInstances.size())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const PfoList *pPfoList1 = nullptr;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*pandoraInstances.at(0), pPfoList1));

    const PfoList *pPfoList2 = nullptr;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*pandoraInstances.at(1), pPfoList2));

    PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), pPfoList1, "InputPfoList1", BLUE));
    PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), pPfoList2, "InputPfoList2", GREEN));
    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    const PfoList *pNewPfoList = nullptr;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*pAlgorithm, pNewPfoList));
    PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), pNewPfoList, "NewPfoList", RED));
    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode StitchingPfoMergingTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
