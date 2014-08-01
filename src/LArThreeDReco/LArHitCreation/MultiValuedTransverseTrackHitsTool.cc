/**
 *  @file   LArContent/src/LArThreeDReco/LArHitCreation/MultiValuedTransverseTrackHitsTool.cc
 * 
 *  @brief  Implementation of the multivalued transverse track hit creation tool.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArThreeDReco/LArHitCreation/MultiValuedTransverseTrackHitsTool.h"

using namespace pandora;

namespace lar
{
 
void MultiValuedTransverseTrackHitsTool::GetThreeDPosition(const CaloHit *const pCaloHit2D, const MatchedSlidingFitMap &matchedSlidingFitMap,
    CartesianVector &position3D, float &chiSquared) const
{  
    const HitType hitType(pCaloHit2D->GetHitType());
    const HitType hitType1((TPC_VIEW_U == hitType) ? TPC_VIEW_V : (TPC_VIEW_V == hitType) ? TPC_VIEW_W : TPC_VIEW_U);
    const HitType hitType2((TPC_VIEW_U == hitType) ? TPC_VIEW_W : (TPC_VIEW_V == hitType) ? TPC_VIEW_U : TPC_VIEW_V);

    MatchedSlidingFitMap::const_iterator fIter1 = matchedSlidingFitMap.find(hitType1);
    MatchedSlidingFitMap::const_iterator fIter2 = matchedSlidingFitMap.find(hitType2);

    if (matchedSlidingFitMap.end() == fIter1 || matchedSlidingFitMap.end() == fIter2)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    const TwoDSlidingFitResult &fitResult1 = fIter1->second;
    const TwoDSlidingFitResult &fitResult2 = fIter2->second;

    CartesianPointList fitPositionList1, fitPositionList2;
    fitResult1.GetGlobalFitPositionListAtX(pCaloHit2D->GetPositionVector().GetX(), fitPositionList1);
    fitResult2.GetGlobalFitPositionListAtX(pCaloHit2D->GetPositionVector().GetX(), fitPositionList2);

    if (fitPositionList1.empty() || fitPositionList2.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    this->GetBestPosition3D(pCaloHit2D, hitType1, hitType2, fitPositionList1, fitPositionList2, position3D, chiSquared);
}

} // namespace lar
