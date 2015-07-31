/**
 *  @file   LArContent/src/LArThreeDReco/LArEventBuilding/BranchAssociatedPfosTool.cc
 * 
 *  @brief  Implementation of the branch associated pfos tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArThreeDReco/LArEventBuilding/BranchAssociatedPfosTool.h"

using namespace pandora;

namespace lar_content
{

void BranchAssociatedPfosTool::Run(PfoHierarchyAlgorithm *const pAlgorithm, const Vertex *const /*pNeutrinoVertex*/,
    PfoHierarchyAlgorithm::PfoInfoMap &/*pfoInfoMap*/)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << this->GetType() << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode BranchAssociatedPfosTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
