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

using namespace pandora;

namespace lar
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
            float chiSquared(std::numeric_limits<float>::max());
            this->GetThreeDPosition(pCaloHit2D, matchedSlidingFitMap, position3D, chiSquared);

            if (chiSquared > m_chiSquaredCut)
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

} // namespace lar
