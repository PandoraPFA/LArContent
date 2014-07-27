/**
 *  @file   LArContent/src/LArThreeDReco/LArHitCreation/ClearTransverseTrackHitsTool.cc
 * 
 *  @brief  Implementation of the transverse track hit creation tool.
 * 
 *  $Log: $
 */


#include "LArHelpers/LArClusterHelper.h"

#include "LArThreeDReco/LArHitCreation/ClearTransverseTrackHitsTool.h"

using namespace pandora;

namespace lar
{
 
void ClearTransverseTrackHitsTool::GetThreeDPosition(const CaloHit *const pCaloHit2D, const MatchedSlidingFitMap &matchedSlidingFitMap,
    CartesianVector &position3D, float &chiSquared) const
{  
    const HitType hitType(pCaloHit2D->GetHitType());
    const HitType hitType1((TPC_VIEW_U == hitType) ? TPC_VIEW_V : (TPC_VIEW_V == hitType) ? TPC_VIEW_W : TPC_VIEW_U);
    const HitType hitType2((TPC_VIEW_U == hitType) ? TPC_VIEW_W : (TPC_VIEW_V == hitType) ? TPC_VIEW_U : TPC_VIEW_V);

    bool foundPosition1(false), foundPosition2(false);
    CartesianVector position1(0.f, 0.f, 0.f), position2(0.f, 0.f, 0.f);
   
    try
    {
        MatchedSlidingFitMap::const_iterator fIter1 = matchedSlidingFitMap.find(hitType1);
        if (matchedSlidingFitMap.end() == fIter1)
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    
        const TwoDSlidingFitResult &fitResult1 = fIter1->second;
        this->GetExtrapolatedPositionAtX(fitResult1, pCaloHit2D->GetPositionVector().GetX(), position1);
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
        this->GetExtrapolatedPositionAtX(fitResult2, pCaloHit2D->GetPositionVector().GetX(), position2);
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
        // TODO --- JUST COMBINE TWO VIEWS

        chiSquared = 0.f;
    }
    else if(foundPosition2)
    {
        // TODO --- JUST COMBINE TWO VIEWS

        chiSquared = 0.f;
    }
    else
    {
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClearTransverseTrackHitsTool::GetExtrapolatedPositionAtX(const TwoDSlidingFitResult &fitResult, const float x, CartesianVector &position2D) const
{
    
    // TODO --- EXTRAPOLATE BEYOND RANGE OF SLIDING FIT RESULT

    // fitResult.GetExtrapolatedPositionAtX(x, position2D);

    fitResult.GetGlobalFitPositionAtX(x, position2D);
   
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClearTransverseTrackHitsTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    return TransverseTrackHitsBaseTool::ReadSettings(xmlHandle);
}

} // namespace lar
