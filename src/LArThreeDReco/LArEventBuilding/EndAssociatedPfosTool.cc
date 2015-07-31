/**
 *  @file   LArContent/src/LArThreeDReco/LArEventBuilding/EndAssociatedPfosTool.cc
 * 
 *  @brief  Implementation of the end associated pfos tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArThreeDReco/LArEventBuilding/EndAssociatedPfosTool.h"

using namespace pandora;

namespace lar_content
{

void EndAssociatedPfosTool::Run(PfoHierarchyAlgorithm *const pAlgorithm, const Vertex *const /*pNeutrinoVertex*/,
    PfoHierarchyAlgorithm::PfoInfoMap &/*pfoInfoMap*/)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << this->GetType() << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EndAssociatedPfosTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
