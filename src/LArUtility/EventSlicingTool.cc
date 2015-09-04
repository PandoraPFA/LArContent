/**
 *  @file   LArContent/src/LArUtility/EventSlicingTool.cc
 * 
 *  @brief  Implementation of the event slicing tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArUtility/EventSlicingTool.h"

using namespace pandora;

namespace lar_content
{

typedef NeutrinoParentAlgorithm::HitTypeToNameMap HitTypeToNameMap;
typedef NeutrinoParentAlgorithm::SliceList SliceList;

void EventSlicingTool::Slice(NeutrinoParentAlgorithm *const pAlgorithm, const HitTypeToNameMap &/*caloHitListNames*/,
    const HitTypeToNameMap &/*clusterListNames*/, SliceList &/*sliceList*/)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << this->GetType() << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventSlicingTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
