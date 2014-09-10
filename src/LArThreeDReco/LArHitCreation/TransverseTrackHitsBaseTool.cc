/**
 *  @file   LArContent/src/LArThreeDReco/LArHitCreation/TransverseTrackHitsBaseTool.cc
 * 
 *  @brief  Implementation of the transverse track hits base tool.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.h"
#include "LArThreeDReco/LArHitCreation/TransverseTrackHitsBaseTool.h"

#include "LArHelpers/LArGeometryHelper.h"

using namespace pandora;

namespace lar_content
{

void TransverseTrackHitsBaseTool::CreateThreeDHits(ThreeDHitCreationAlgorithm *pAlgorithm, const CaloHitList &inputTwoDHits, 
    const MatchedSlidingFitMap &matchedSlidingFitMap, CaloHitList &newThreeDHits) const
{   
    for (CaloHitList::const_iterator iter = inputTwoDHits.begin(), iterEnd = inputTwoDHits.end(); iter != iterEnd; ++iter)
    {
        try
        {
            CaloHit *pCaloHit2D(*iter);

            CartesianVector position3D(0.f, 0.f, 0.f);
            float chiSquared1(std::numeric_limits<float>::max()), chiSquared2(0.f);
            this->GetThreeDPosition(pCaloHit2D, matchedSlidingFitMap, position3D, chiSquared1);
            this->GetTransverseChi2(pCaloHit2D, matchedSlidingFitMap, position3D, chiSquared2);

            if (chiSquared1 + chiSquared2 > m_chiSquaredCut)
                throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);

            CaloHit *pCaloHit3D(NULL);
            pAlgorithm->CreateThreeDHit(pCaloHit2D, position3D, pCaloHit3D);
            newThreeDHits.insert(pCaloHit3D);
        }
        catch (StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseTrackHitsBaseTool::GetTransverseChi2(const CaloHit *const pCaloHit2D, const MatchedSlidingFitMap &matchedSlidingFitMap,
    const CartesianVector &position3D, float &chiSquared) const
{  
    // TODO Develop a proper treatment of the |dz/dx| * sigmaX uncertainty
    chiSquared = 0.f;

    for (MatchedSlidingFitMap::const_iterator iter = matchedSlidingFitMap.begin(), iterEnd = matchedSlidingFitMap.end();
        iter != iterEnd; ++iter)
    {
        const HitType hitType(iter->first);

        if (pCaloHit2D->GetHitType() == hitType)
            continue;

        const CartesianVector position2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), position3D, hitType));
        const TwoDSlidingFitResult &fitResult(iter->second);

        chiSquared += this->GetTransverseChi2(position2D, fitResult);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float TransverseTrackHitsBaseTool::GetTransverseChi2(const CartesianVector &position2D, const TwoDSlidingFitResult &fitResult) const
{
    const float minLayerX(fitResult.GetGlobalMinLayerPosition().GetX());
    const float maxLayerX(fitResult.GetGlobalMaxLayerPosition().GetX());

    const float minX(std::min(minLayerX, maxLayerX));
    const float maxX(std::max(minLayerX, maxLayerX));

    if (((position2D.GetX() - minX) > -std::numeric_limits<float>::epsilon()) && 
        ((position2D.GetX() - maxX) < +std::numeric_limits<float>::epsilon()))
        return 0.f;

    const float minLayerDistanceSquared((position2D - fitResult.GetGlobalMinLayerPosition()).GetMagnitudeSquared());
    const float maxLayerDistanceSquared((position2D - fitResult.GetGlobalMaxLayerPosition()).GetMagnitudeSquared());

    float px(0.f), pz(0.f);

    if (minLayerDistanceSquared < maxLayerDistanceSquared)
    {
        px = fitResult.GetGlobalMinLayerDirection().GetX();
        pz = fitResult.GetGlobalMinLayerDirection().GetZ();
    }
    else if(maxLayerDistanceSquared < minLayerDistanceSquared)
    {
        px = fitResult.GetGlobalMaxLayerDirection().GetX();
        pz = fitResult.GetGlobalMaxLayerDirection().GetZ();
    }

    if (std::fabs(px) > std::numeric_limits<float>::epsilon())
    {
        return (0.5f * (pz * pz) / (px * px)); // TODO: Figure out the appropriate scale factor
    }

    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

} // namespace lar_content
