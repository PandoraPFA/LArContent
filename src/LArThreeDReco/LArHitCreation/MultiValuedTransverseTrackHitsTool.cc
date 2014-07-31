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

    CartesianPointList fitPositions1, fitPositions2;
    fitResult1.GetGlobalFitPositionListAtX(pCaloHit2D->GetPositionVector().GetX(), fitPositions1);
    fitResult2.GetGlobalFitPositionListAtX(pCaloHit2D->GetPositionVector().GetX(), fitPositions2);

    bool foundPosition(false);
    CartesianVector bestPosition3D(0.f, 0.f, 0.f);
    float bestChiSquared(std::numeric_limits<float>::max());

    for (CartesianPointList::const_iterator iter1 = fitPositions1.begin(), iter1End = fitPositions1.end(); iter1 != iter1End; ++iter1)
    {
        for (CartesianPointList::const_iterator iter2 = fitPositions2.begin(), iter2End = fitPositions2.end(); iter2 != iter2End; ++iter2)
        {
            CartesianVector thisPosition3D(0.f, 0.f, 0.f);
            float thisChiSquared(std::numeric_limits<float>::max());
            this->GetPosition3D(pCaloHit2D, hitType1, hitType2, *iter1, *iter2, thisPosition3D, thisChiSquared);

            if (thisChiSquared < bestChiSquared)
            {
                foundPosition = true;
                bestPosition3D = thisPosition3D;
                bestChiSquared = thisChiSquared;
            }
        }
    }

    if (!foundPosition)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    position3D = bestPosition3D;
    chiSquared = bestChiSquared;
}

} // namespace lar
