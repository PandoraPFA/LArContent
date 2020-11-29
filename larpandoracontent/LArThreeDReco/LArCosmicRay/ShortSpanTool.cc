/**
 *  @file   larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/ShortSpanTool.cc
 *
 *  @brief  Implementation of the clear tracks tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArThreeDReco/LArCosmicRay/ShortSpanTool.h"

using namespace pandora;

namespace lar_content
{

ShortSpanTool::ShortSpanTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ShortSpanTool::Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    bool changesMade(false);

    TensorType::ElementList elementList;
    overlapTensor.GetUnambiguousElements(true, elementList);
    this->InvestigateShortSpans(pAlgorithm, elementList, particlesMade);

    return particlesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShortSpanTool::InvestigateShortSpans(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::ElementList &elementList,
    bool &changesMade) const
{

    TensorType::ElementList elementList;
    overlapTensor.GetUnambiguousElements(true, elementList);

    for (const TensorType::Element &element : elementList)
    {
        const Cluster *pShortCluster(null);
        
        if (!this->GetShortCluster(element, pShortCluster))
            continue;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ShortSpanTool::GetShortCluster(const TensorType::Element &element, const Cluster *&pShortCluster) const
{

}

//------------------------------------------------------------------------------------------------------------------------------------------    
    
StatusCode ShortSpanTool::ReadSettings(const TiXmlHandle xmlHandle)
{

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
