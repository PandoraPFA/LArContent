/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/ClearLongitudinalTrackHitsTool.cc
 *
 *  @brief  Implementation of the clear longitudinal track hit creation tool.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArThreeDReco/LArHitCreation/ClearLongitudinalTrackHitsTool.h"

using namespace pandora;

namespace lar_content
{

void ClearLongitudinalTrackHitsTool::GetThreeDPosition(const CaloHit *const pCaloHit2D, const MatchedSlidingFitMap &matchedSlidingFitMap,
    const CartesianVector &vtx3D, const CartesianVector &end3D, CartesianVector &position3D, float &chiSquared) const
{
    const HitType hitType(pCaloHit2D->GetHitType());
    const HitType hitType1((TPC_VIEW_U == hitType) ? TPC_VIEW_V : (TPC_VIEW_V == hitType) ? TPC_VIEW_W : TPC_VIEW_U);
    const HitType hitType2((TPC_VIEW_U == hitType) ? TPC_VIEW_W : (TPC_VIEW_V == hitType) ? TPC_VIEW_U : TPC_VIEW_V);

    const CartesianVector vtx2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), vtx3D, hitType));
    const CartesianVector end2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), end3D, hitType));

    if ((end2D - vtx2D).GetMagnitudeSquared() < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    const float frac((end2D - vtx2D).GetDotProduct(pCaloHit2D->GetPositionVector() - vtx2D) / (end2D - vtx2D).GetMagnitudeSquared());
    const CartesianVector projection3D(vtx3D + (end3D - vtx3D) * frac);

    CartesianPointList fitPositionList1, fitPositionList2;

    MatchedSlidingFitMap::const_iterator fIter1 = matchedSlidingFitMap.find(hitType1);
    if (matchedSlidingFitMap.end() != fIter1)
    {
        const TwoDSlidingFitResult &fitResult1 = fIter1->second;
        const CartesianVector position2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), projection3D, hitType1));

        float rL1(0.f), rT1(0.f);
        CartesianVector position1(0.f, 0.f, 0.f);
        fitResult1.GetLocalPosition(position2D, rL1, rT1);
        const StatusCode statusCode(fitResult1.GetTransverseProjection(pCaloHit2D->GetPositionVector().GetX(), fitResult1.GetFitSegment(rL1), position1));

        if ((STATUS_CODE_SUCCESS != statusCode) && (STATUS_CODE_NOT_FOUND != statusCode))
            throw StatusCodeException(statusCode);

        if (STATUS_CODE_SUCCESS == statusCode)
            fitPositionList1.push_back(position1);
    }

    MatchedSlidingFitMap::const_iterator fIter2 = matchedSlidingFitMap.find(hitType2);
    if (matchedSlidingFitMap.end() != fIter2)
    {
        const TwoDSlidingFitResult &fitResult2 = fIter2->second;
        const CartesianVector position2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), projection3D, hitType2));

        float rL2(0.f), rT2(0.f);
        CartesianVector position2(0.f, 0.f, 0.f);
        fitResult2.GetLocalPosition(position2D, rL2, rT2);
        const StatusCode statusCode(fitResult2.GetTransverseProjection(pCaloHit2D->GetPositionVector().GetX(), fitResult2.GetFitSegment(rL2), position2));

        if ((STATUS_CODE_SUCCESS != statusCode) && (STATUS_CODE_NOT_FOUND != statusCode))
            throw StatusCodeException(statusCode);

        if (STATUS_CODE_SUCCESS == statusCode)
            fitPositionList2.push_back(position2);
    }

    unsigned int nViews(1);
    if (fitPositionList1.size() > 0) ++nViews;
    if (fitPositionList2.size() > 0) ++nViews;

    if (nViews < m_minViews)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    this->GetBestPosition3D(pCaloHit2D, hitType1, hitType2, fitPositionList1, fitPositionList2, position3D, chiSquared);  
}

} // namespace lar_content
