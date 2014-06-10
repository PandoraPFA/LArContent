/**
 *  @file   LArContent/src/LArThreeDReco/LArTransverseTrackMatching/ThreeDTransverseTrackFragmentsAlg.cc
 *
 *  @brief  Implementation of the three dimensional transverse track fragments algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"
#include "LArHelpers/LArThreeDHelper.h"

#include "LArThreeDReco/LArTransverseTrackMatching/ThreeDTransverseTrackFragmentsAlg.h"

using namespace pandora;

namespace lar
{

void ThreeDTransverseTrackFragmentsAlg::GetProjectedPositions(const TwoDSlidingFitResult &fitResult1, const TwoDSlidingFitResult &fitResult2,
    CartesianPointList &projectedPositions) const
{
    const Cluster *pCluster1(fitResult1.GetCluster());
    const Cluster *pCluster2(fitResult2.GetCluster());

    const HitType hitType1(LArThreeDHelper::GetClusterHitType(pCluster1));
    const HitType hitType2(LArThreeDHelper::GetClusterHitType(pCluster2));

    float xMin1(0.f), xMax1(0.f), xMin2(0.f), xMax2(0.f);
    LArClusterHelper::GetClusterSpanX(pCluster1, xMin1, xMax1);
    LArClusterHelper::GetClusterSpanX(pCluster2, xMin2, xMax2);

    const float xMin(std::max(xMin1, xMin2));
    const float xMax(std::min(xMax1, xMax2));
    const float nSamplingPoints(3.f * (xMax - xMin) / m_maxPointDisplacement);

    if (nSamplingPoints < 1.f)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    for (float iSample = 0.f; iSample < nSamplingPoints; iSample += 1.f)
    {
        const float x(xMin + (xMax - xMin) * (0.5f + iSample) / nSamplingPoints);

        try
        {
            float chi2(0.f);
            CartesianVector position1(0.f, 0.f, 0.f), position2(0.f, 0.f, 0.f), position3(0.f, 0.f, 0.f);
            fitResult1.GetGlobalFitPositionAtX(x, position1);
            fitResult2.GetGlobalFitPositionAtX(x, position2);
            LArGeometryHelper::MergeTwoPositions(hitType1, hitType2, position1, position2, position3, chi2);
            projectedPositions.push_back(position3);
        }
        catch (StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDTransverseTrackFragmentsAlg::ReadSettings(const TiXmlHandle xmlHandle)
{
    return ThreeDFragmentsBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
