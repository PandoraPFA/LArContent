/**
 *  @file   LArContent/src/LArThreeDReco/LArEventBuilding/VertexAssociatedPfosTool.cc
 * 
 *  @brief  Implementation of the vertex associated pfos tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArThreeDReco/LArEventBuilding/VertexAssociatedPfosTool.h"

using namespace pandora;

namespace lar_content
{

void VertexAssociatedPfosTool::Run(PfoHierarchyAlgorithm *const pAlgorithm, const Vertex *const /*pNeutrinoVertex*/,
    PfoHierarchyAlgorithm::PfoInfoMap &/*pfoInfoMap*/)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << this->GetType() << std::endl;

    
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexAssociatedPfosTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
