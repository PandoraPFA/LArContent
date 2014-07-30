/**
 *  @file   LArContent/src/LArThreeDReco/LArHitCreation/MultiValuedLongitudinalTrackHitsTool.cc
 *
 *  @brief  Implementation of the longitudinal track hit creation tool.
 *
 *  $Log: $
 */

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"

#include "LArThreeDReco/LArHitCreation/MultiValuedLongitudinalTrackHitsTool.h"

using namespace pandora;

namespace lar
{

void MultiValuedLongitudinalTrackHitsTool::GetThreeDPosition(const CaloHit *const pCaloHit2D, const MatchedSlidingFitMap &matchedSlidingFitMap,
    const CartesianVector &vtx3D, const CartesianVector &end3D, CartesianVector &position3D, float &chiSquared) const
{
    const HitType hitType(pCaloHit2D->GetHitType());
    const HitType hitType1((TPC_VIEW_U == hitType) ? TPC_VIEW_V : (TPC_VIEW_V == hitType) ? TPC_VIEW_W : TPC_VIEW_U);
    const HitType hitType2((TPC_VIEW_U == hitType) ? TPC_VIEW_W : (TPC_VIEW_V == hitType) ? TPC_VIEW_U : TPC_VIEW_V);

    const CartesianVector vtx2D(LArGeometryHelper::ProjectPosition(vtx3D, hitType));
    const CartesianVector end2D(LArGeometryHelper::ProjectPosition(end3D, hitType));

    if ((end2D - vtx2D).GetMagnitudeSquared() < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    const float frac((end2D - vtx2D).GetDotProduct(pCaloHit2D->GetPositionVector() - vtx2D) / (end2D - vtx2D).GetMagnitudeSquared());
    const CartesianVector projection3D(vtx3D + (end3D - vtx3D) * frac);

    bool foundPosition1(false), foundPosition2(false);
    CartesianVector position1(0.f, 0.f, 0.f), position2(0.f, 0.f, 0.f);

    try
    {
        MatchedSlidingFitMap::const_iterator fIter1 = matchedSlidingFitMap.find(hitType1);
        if (matchedSlidingFitMap.end() == fIter1)
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);

        const TwoDSlidingFitResult &fitResult1 = fIter1->second;
        const CartesianVector position2D(LArGeometryHelper::ProjectPosition(projection3D, hitType1));

        fitResult1.GetGlobalFitProjection(position2D, position1);
        foundPosition1 = true;
    }
    catch (StatusCodeException &)
    {
    }

    try
    {
        MatchedSlidingFitMap::const_iterator fIter2 = matchedSlidingFitMap.find(hitType2);
        if (matchedSlidingFitMap.end() == fIter2)
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);

        const TwoDSlidingFitResult &fitResult2 = fIter2->second;
        const CartesianVector position2D(LArGeometryHelper::ProjectPosition(projection3D, hitType2));

        fitResult2.GetGlobalFitProjection(position2D, position2);
        foundPosition2 = true;
    }
    catch (StatusCodeException &)
    {
    }

    if (foundPosition1 && foundPosition2)
    {
        this->GetPosition3D(pCaloHit2D, hitType1, hitType2, position1, position2, position3D, chiSquared);
    }
    else if(foundPosition1)
    {
        this->GetPosition3D(pCaloHit2D, hitType1, position1, position3D, chiSquared);
    }
    else if(foundPosition2)
    {
        this->GetPosition3D(pCaloHit2D, hitType2, position2, position3D, chiSquared);
    }
    else
    {
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MultiValuedLongitudinalTrackHitsTool::ReadSettings(const TiXmlHandle xmlHandle)
{
  return LongitudinalTrackHitsBaseTool::ReadSettings(xmlHandle);
}

} // namespace lar
