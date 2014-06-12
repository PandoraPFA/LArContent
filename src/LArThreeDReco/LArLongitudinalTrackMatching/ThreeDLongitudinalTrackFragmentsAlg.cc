/**
 *  @file   LArContent/src/LArThreeDReco/LArLongitudinalTrackMatching/ThreeDLongitudinalTrackFragmentsAlg.cc
 *
 *  @brief  Implementation of the three dimensional longitudinal track fragments algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"
#include "LArHelpers/LArThreeDHelper.h"

#include "LArThreeDReco/LArLongitudinalTrackMatching/ThreeDLongitudinalTrackFragmentsAlg.h"

using namespace pandora;

namespace lar
{

void ThreeDLongitudinalTrackFragmentsAlg::GetProjectedPositions(const TwoDSlidingFitResult&, const TwoDSlidingFitResult&, CartesianPointList&) const
{
  
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDLongitudinalTrackFragmentsAlg::ReadSettings(const TiXmlHandle xmlHandle)
{
    return ThreeDFragmentsBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
