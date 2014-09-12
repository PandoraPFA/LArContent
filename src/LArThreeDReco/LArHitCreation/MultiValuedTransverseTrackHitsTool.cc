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

namespace lar_content
{
 
void MultiValuedTransverseTrackHitsTool::GetThreeDPosition(const CaloHit *const pCaloHit2D, const MatchedSlidingFitMap &matchedSlidingFitMap,
    CartesianVector &position3D, float &chiSquared) const
{  
    const HitType hitType(pCaloHit2D->GetHitType());
    const HitType hitType1((TPC_VIEW_U == hitType) ? TPC_VIEW_V : (TPC_VIEW_V == hitType) ? TPC_VIEW_W : TPC_VIEW_U);
    const HitType hitType2((TPC_VIEW_U == hitType) ? TPC_VIEW_W : (TPC_VIEW_V == hitType) ? TPC_VIEW_U : TPC_VIEW_V);

    CartesianPointList fitPositionList1, fitPositionList2;

    try
    {
        MatchedSlidingFitMap::const_iterator iter1 = matchedSlidingFitMap.find(hitType1);
        if (matchedSlidingFitMap.end() == iter1)
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);

        const TwoDSlidingFitResult &fitResult1 = iter1->second;
        fitResult1.GetGlobalFitPositionListAtX(pCaloHit2D->GetPositionVector().GetX(), fitPositionList1);
    }
    catch (StatusCodeException &)
    {
    }

    try
    {
        MatchedSlidingFitMap::const_iterator iter2 = matchedSlidingFitMap.find(hitType2);
        if (matchedSlidingFitMap.end() == iter2)
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);

        const TwoDSlidingFitResult &fitResult2 = iter2->second;
        fitResult2.GetGlobalFitPositionListAtX(pCaloHit2D->GetPositionVector().GetX(), fitPositionList2); 
    }
    catch (StatusCodeException &)
    {
    }

    unsigned int nViews(1);
    if (fitPositionList1.size() > 0) ++nViews;
    if (fitPositionList2.size() > 0) ++nViews;

    if (nViews < m_minViews)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    if (nViews <= 2 && !m_useDeltaXCorrection)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    this->GetBestPosition3D(pCaloHit2D, hitType1, hitType2, fitPositionList1, fitPositionList2, position3D, chiSquared);
}

} // namespace lar_content
