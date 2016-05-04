/**
 *  @file   LArContent/src/LArThreeDReco/LArHitCreation/ClearTransverseTrackHitsTool.cc
 *
 *  @brief  Implementation of the transverse track hit creation tool.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArContent/LArThreeDReco/LArHitCreation/ClearTransverseTrackHitsTool.h"

using namespace pandora;

namespace lar_content
{

void ClearTransverseTrackHitsTool::GetThreeDPosition(const CaloHit *const pCaloHit2D, const MatchedSlidingFitMap &matchedSlidingFitMap,
    CartesianVector &position3D, float &chiSquared) const
{
    const HitType hitType(pCaloHit2D->GetHitType());
    const HitType hitType1((TPC_VIEW_U == hitType) ? TPC_VIEW_V : (TPC_VIEW_V == hitType) ? TPC_VIEW_W : TPC_VIEW_U);
    const HitType hitType2((TPC_VIEW_U == hitType) ? TPC_VIEW_W : (TPC_VIEW_V == hitType) ? TPC_VIEW_U : TPC_VIEW_V);

    CartesianPointList fitPositionList1, fitPositionList2;

    MatchedSlidingFitMap::const_iterator fIter1 = matchedSlidingFitMap.find(hitType1);
    if (matchedSlidingFitMap.end() != fIter1)
    {
        const TwoDSlidingFitResult &fitResult1 = fIter1->second;
        CartesianVector position1(0.f, 0.f, 0.f);
        const StatusCode statusCode(fitResult1.GetExtrapolatedPositionAtX(pCaloHit2D->GetPositionVector().GetX(), position1));

        if ((STATUS_CODE_SUCCESS != statusCode) && (STATUS_CODE_NOT_FOUND != statusCode))
            throw StatusCodeException(statusCode);

        if (STATUS_CODE_SUCCESS == statusCode)
            fitPositionList1.push_back(position1);
    }

    MatchedSlidingFitMap::const_iterator fIter2 = matchedSlidingFitMap.find(hitType2);
    if (matchedSlidingFitMap.end() != fIter2)
    {
        const TwoDSlidingFitResult &fitResult2 = fIter2->second;
        CartesianVector position2(0.f, 0.f, 0.f);
        const StatusCode statusCode(fitResult2.GetExtrapolatedPositionAtX(pCaloHit2D->GetPositionVector().GetX(), position2));

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
