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

void ThreeDLongitudinalTrackFragmentsAlg::GetProjectedPositions(const TwoDSlidingFitResult &fitResult1, const TwoDSlidingFitResult &fitResult2,
    CartesianPointList &projectedPositions) const
{
    const Cluster *pCluster1(fitResult1.GetCluster());
    const Cluster *pCluster2(fitResult2.GetCluster());

    const HitType hitType1(LArThreeDHelper::GetClusterHitType(pCluster1));
    const HitType hitType2(LArThreeDHelper::GetClusterHitType(pCluster2));
    const HitType hitType3((TPC_VIEW_U == hitType1 && TPC_VIEW_V == hitType2) ? TPC_VIEW_W :
                           (TPC_VIEW_V == hitType1 && TPC_VIEW_W == hitType2) ? TPC_VIEW_U :
                           (TPC_VIEW_W == hitType1 && TPC_VIEW_U == hitType2) ? TPC_VIEW_V : CUSTOM); 

    if (CUSTOM == hitType3)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const CartesianVector minPosition1(fitResult1.GetGlobalMinLayerPosition());
    const CartesianVector maxPosition1(fitResult1.GetGlobalMaxLayerPosition());
    const CartesianVector minPosition2(fitResult2.GetGlobalMinLayerPosition());
    const CartesianVector maxPosition2(fitResult2.GetGlobalMaxLayerPosition());

    const float dx_A(std::fabs(minPosition2.GetX() - minPosition1.GetX()));
    const float dx_B(std::fabs(maxPosition2.GetX() - maxPosition1.GetX()));
    const float dx_C(std::fabs(maxPosition2.GetX() - minPosition1.GetX()));
    const float dx_D(std::fabs(minPosition2.GetX() - maxPosition1.GetX()));

    if (std::min(dx_C,dx_D) > std::max(dx_A,dx_B) && std::min(dx_A,dx_B) > std::max(dx_C,dx_D))
        return;

    const CartesianVector &vtxPosition1(minPosition1);
    const CartesianVector &endPosition1(maxPosition1);
    const CartesianVector &vtxPosition2((dx_A < dx_C) ? minPosition2 : maxPosition2);
    const CartesianVector &endPosition2((dx_A < dx_C) ? maxPosition2 : minPosition2);

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

    for (float iSample = 0.f; iSample < nSamplingPoints; iSample += 1.f)
    {
        const CartesianVector linearPosition3D(vtxPosition3D + (endPosition3D - vtxPosition3D) * ((0.5f + iSample) / nSamplingPoints));
        const CartesianVector linearPosition1(LArGeometryHelper::ProjectPosition(linearPosition3D, hitType1));
        const CartesianVector linearPosition2(LArGeometryHelper::ProjectPosition(linearPosition3D, hitType2));

        try
        {
            float chi2(0.f);
            CartesianVector position1(0.f, 0.f, 0.f), position2(0.f, 0.f, 0.f), position3(0.f, 0.f, 0.f);
            fitResult1.GetGlobalFitProjection(linearPosition1, position1);
            fitResult2.GetGlobalFitProjection(linearPosition2, position2); 
            LArGeometryHelper::MergeTwoPositions(hitType1, hitType2, position1, position2, position3, chi2);
            projectedPositions.push_back(position3); 
        }
        catch (StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDLongitudinalTrackFragmentsAlg::ReadSettings(const TiXmlHandle xmlHandle)
{
    return ThreeDFragmentsBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
