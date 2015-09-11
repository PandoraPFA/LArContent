/**
 *  @file   LArContent/src/LArThreeDReco/LArEventBuilding/EventSlicingTool.cc
 * 
 *  @brief  Implementation of the event slicing tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArThreeDReco/LArEventBuilding/EventSlicingTool.h"

using namespace pandora;

namespace lar_content
{

typedef NeutrinoParentAlgorithm::SliceList SliceList;
typedef NeutrinoParentAlgorithm::HitTypeToNameMap HitTypeToNameMap;

void EventSlicingTool::Slice(NeutrinoParentAlgorithm *const pAlgorithm, const HitTypeToNameMap &caloHitListNames,
    const HitTypeToNameMap &/*clusterListNames*/, SliceList &sliceList)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << this->GetType() << std::endl;

    // TODO - This just puts all the input hits in a single slice
    sliceList.push_back(NeutrinoParentAlgorithm::Slice());
    NeutrinoParentAlgorithm::Slice &sliceHack(sliceList.at(0));

    const CaloHitList *pTempListU(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*pAlgorithm, caloHitListNames.at(TPC_VIEW_U), pTempListU));
    sliceHack.m_caloHitListU = *pTempListU;

    const CaloHitList *pTempListV(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*pAlgorithm, caloHitListNames.at(TPC_VIEW_V), pTempListV));
    sliceHack.m_caloHitListV = *pTempListV;

    const CaloHitList *pTempListW(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*pAlgorithm, caloHitListNames.at(TPC_VIEW_W), pTempListW));
    sliceHack.m_caloHitListW = *pTempListW;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventSlicingTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
