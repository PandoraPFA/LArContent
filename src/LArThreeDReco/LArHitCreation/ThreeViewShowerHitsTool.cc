/**
 *  @file   LArContent/src/LArThreeDReco/LArHitCreation/ThreeViewShowerHitsTool.cc
 * 
 *  @brief  Implementation of the three view shower hits tool.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArThreeDReco/LArHitCreation/ThreeViewShowerHitsTool.h"

using namespace pandora;

namespace lar
{

void ThreeViewShowerHitsTool::GetThreeDPosition(const CaloHit *const pCaloHit2D, const CaloHitList &caloHitList1, const CaloHitList &caloHitList2,
    CartesianVector &position3D, float &chiSquared) const
{
    if (caloHitList1.empty() || caloHitList2.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    const HitType hitType1((*caloHitList1.begin())->GetHitType());
    const HitType hitType2((*caloHitList2.begin())->GetHitType());

    CartesianPointList hitPositionList1, hitPositionList2;

    for (CaloHitList::const_iterator iter1 = caloHitList1.begin(), iterEnd1 = caloHitList1.end(); iter1 != iterEnd1; ++iter1)
    {
        hitPositionList1.push_back((*iter1)->GetPositionVector());
    }

    for (CaloHitList::const_iterator iter2 = caloHitList2.begin(), iterEnd2 = caloHitList2.end(); iter2 != iterEnd2; ++iter2)
    {
        hitPositionList2.push_back((*iter2)->GetPositionVector());
    }

    this->GetBestPosition3D(pCaloHit2D, hitType1, hitType2, hitPositionList1, hitPositionList2, position3D, chiSquared); 
}

} // namespace lar
