/**
 *  @file   larpandoracontent/LArMonitoring/TwoViewTransverseTracksValidationTool.cc
 *
 *  @brief  Implementation of the two view transverse tracks validation tool.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/TwoViewTransverseTracksValidationTool.h"

using namespace pandora;

namespace lar_content
{

TwoViewTransverseTracksValidationTool::TwoViewTransverseTracksValidationTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewTransverseTracksValidationTool::Run(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2,
    const DiscreteProbabilityVector &discreteProbabilityVector1, const DiscreteProbabilityVector &discreteProbabilityVector2,
    const TwoViewTransverseOverlapResult &overlapResult)
{
    std::cout<<"This runs"<<std::endl;
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoViewTransverseTracksValidationTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
