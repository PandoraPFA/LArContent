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

    // Check hit types
    const HitType hitType1(LArThreeDHelper::GetClusterHitType(pCluster1));
    const HitType hitType2(LArThreeDHelper::GetClusterHitType(pCluster2));
    const HitType hitType3((TPC_VIEW_U != hitType1 && TPC_VIEW_U != hitType2) ? TPC_VIEW_U :
                           (TPC_VIEW_V != hitType1 && TPC_VIEW_V != hitType2) ? TPC_VIEW_V :
                           (TPC_VIEW_W != hitType1 && TPC_VIEW_W != hitType2) ? TPC_VIEW_W : CUSTOM); 

    if (CUSTOM == hitType3)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    // Check absolute and fractional overlap in x coordinate
    float xMin1(0.f), xMax1(0.f), xMin2(0.f), xMax2(0.f);
    LArClusterHelper::GetClusterSpanX(pCluster1, xMin1, xMax1);
    LArClusterHelper::GetClusterSpanX(pCluster2, xMin2, xMax2);
 
    const float xOverlap(std::min(xMax1, xMax2) - std::max(xMin1, xMin2));
    const float xSpan(std::max(xMax1, xMax2) - std::min(xMin1, xMin2));

    if ((xOverlap < m_minXOverlap) || (xSpan < std::numeric_limits<float>::epsilon()) || ((xOverlap / xSpan) < m_minXOverlapFraction))
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
 
    // Identify vertex and end positions (2D)
    const CartesianVector minPosition1(fitResult1.GetGlobalMinLayerPosition());
    const CartesianVector maxPosition1(fitResult1.GetGlobalMaxLayerPosition());
    const CartesianVector minPosition2(fitResult2.GetGlobalMinLayerPosition());
    const CartesianVector maxPosition2(fitResult2.GetGlobalMaxLayerPosition());

    const float dx_A(std::fabs(minPosition2.GetX() - minPosition1.GetX()));
    const float dx_B(std::fabs(maxPosition2.GetX() - maxPosition1.GetX()));
    const float dx_C(std::fabs(maxPosition2.GetX() - minPosition1.GetX()));
    const float dx_D(std::fabs(minPosition2.GetX() - maxPosition1.GetX()));

    if (std::min(dx_C,dx_D) > std::max(dx_A,dx_B) && std::min(dx_A,dx_B) > std::max(dx_C,dx_D))
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    const CartesianVector &vtxPosition1(minPosition1);
    const CartesianVector &endPosition1(maxPosition1);
    const CartesianVector &vtxPosition2((dx_A < dx_C) ? minPosition2 : maxPosition2);
    const CartesianVector &endPosition2((dx_A < dx_C) ? maxPosition2 : minPosition2);

    // Calculate vertex and end positions (3D)
    float vtxChi2(0.f);
    CartesianVector vtxPosition3D(0.f, 0.f, 0.f);
    LArGeometryHelper::MergeTwoPositions3D(hitType1, hitType2, vtxPosition1, vtxPosition2, vtxPosition3D, vtxChi2);

    float endChi2(0.f);
    CartesianVector endPosition3D(0.f, 0.f, 0.f);
    LArGeometryHelper::MergeTwoPositions3D(hitType1, hitType2, endPosition1, endPosition2, endPosition3D, endChi2);

    const CartesianVector vtxProjection1(LArGeometryHelper::ProjectPosition(vtxPosition3D, hitType1));
    const CartesianVector vtxProjection2(LArGeometryHelper::ProjectPosition(vtxPosition3D, hitType2));
    const CartesianVector vtxProjection3(LArGeometryHelper::ProjectPosition(vtxPosition3D, hitType3));

    const CartesianVector endProjection1(LArGeometryHelper::ProjectPosition(endPosition3D, hitType1));
    const CartesianVector endProjection2(LArGeometryHelper::ProjectPosition(endPosition3D, hitType2));
    const CartesianVector endProjection3(LArGeometryHelper::ProjectPosition(endPosition3D, hitType3));
    
    const float nSamplingPoints(3.f * (endProjection3 - vtxProjection3).GetMagnitude() / m_maxPointDisplacement);

    if (nSamplingPoints < 1.f)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    // Loop over trajectory points
    for (float iSample = 0.f; iSample < nSamplingPoints; iSample += 1.f)
    {
        const CartesianVector linearPosition3D(vtxPosition3D + (endPosition3D - vtxPosition3D) * ((0.5f + iSample) / nSamplingPoints));
        const CartesianVector linearPosition1(LArGeometryHelper::ProjectPosition(linearPosition3D, hitType1));
        const CartesianVector linearPosition2(LArGeometryHelper::ProjectPosition(linearPosition3D, hitType2));

        try
        {
            float chi2(0.f);
            CartesianVector fitPosition1(0.f, 0.f, 0.f), fitPosition2(0.f, 0.f, 0.f);
            fitResult1.GetGlobalFitProjection(linearPosition1, fitPosition1);
            fitResult2.GetGlobalFitProjection(linearPosition2, fitPosition2); 

            float rL1(0.f), rL2(0.f), rT1(0.f), rT2(0.f);
            fitResult1.GetLocalPosition(fitPosition1, rL1, rT1);
            fitResult2.GetLocalPosition(fitPosition2, rL2, rT2);

            const TwoDSlidingFitResult::FitSegment& fitSegment1 = fitResult1.GetFitSegment(rL1);
            const TwoDSlidingFitResult::FitSegment& fitSegment2 = fitResult2.GetFitSegment(rL2);

            const float x(0.5 * (fitPosition1.GetX() + fitPosition2.GetX()));
            CartesianVector position1(0.f, 0.f, 0.f), position2(0.f, 0.f, 0.f), position3(0.f, 0.f, 0.f);
            fitResult1.GetTransverseProjection(x, fitSegment1, position1);
            fitResult2.GetTransverseProjection(x, fitSegment2, position2);
            LArGeometryHelper::MergeTwoPositions(hitType1, hitType2, position1, position2, position3, chi2);
            projectedPositions.push_back(position3);
        }
        catch (StatusCodeException &)
        {
        }
    }

    // Bail out if list of projected positions is empty
    if (projectedPositions.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDTransverseTrackFragmentsAlg::ReadSettings(const TiXmlHandle xmlHandle)
{   
    m_minXOverlap = 3.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinXOverlap", m_minXOverlap));

    m_minXOverlapFraction = 0.8f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinXOverlapFraction", m_minXOverlapFraction));

    return ThreeDFragmentsBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
