/**
 *  @file   larpandoracontent/LArHelpers/LArDiscreteCumulativeDistributionHelper.cc
 *
 *  @brief  Implementation of the discrete cumulative distribution helper class.
 *
 *  $Log: $
 */

#include <algorithm>
#include <iostream>
#include "larpandoracontent/LArObjects/LArDiscreteCumulativeDistribution.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

namespace lar_content
{

float LArDiscreteCumulativeDistributionHelper::CalculateKSTestStatistic(const DiscreteCumulativeDistribution &distributionA, const DiscreteCumulativeDistribution &distributionB)
{
    float ks = std::numeric_limits<float>::min();
    for (size_t iElement = 0; iElement < distributionA.GetSize(); ++iElement)
    {
        float xA(0), yA(0);
        distributionA.GetXandY(iElement,xA,yA);
        float yB = FindY(distributionB, xA);
        ks = std::max(ks, std::abs(yA-yB));
    }
    for (size_t iElement = 0; iElement < distributionB.GetSize(); ++iElement)
    {
        float xB(0), yB(0);
        distributionB.GetXandY(iElement,xB,yB);
        float yA = FindY(distributionA, xB);
        ks = std::max(ks, std::abs(yB-yA));
    }
    /*
    //float dp = std::numeric_limits<float>::min();
    //float dm = std::numeric_limits<float>::min();
    for (size_t iDist = 0; iDist < distributionA.GetSize(); iDist++)
    {
        float xA(0), xB(0);
        float yA(0), yB(0);
        distributionA.GetXandY(iDist,xA,yA);
        distributionB.GetXandY(iDist,xB,yB);
        //std::cout<<"xA: " << xA << "  xB: " << xB << "  yAIn: " << distributionA.GetYFromInputVector(iDist) << "  yA: " << yA << "  yB: " << yB << "  ks: " << std::abs(yA-yB) << std::endl;
        ks = std::max(ks, std::abs(yA-yB));
        //dp = std::max(dp, yA-yB);
        //dm = std::max(dm, yB-yA);
    }
    */
    return ks;
}

float LArDiscreteCumulativeDistributionHelper::FindY(const DiscreteCumulativeDistribution &distribution, const float &x)
{
    float y = 0;
    for (size_t iElement = 0; iElement < distribution.GetSize(); ++iElement)
    {
        float xElement(0), yElement(0);
        distribution.GetXandY(iElement, xElement, yElement);
        if (xElement > x) 
            break;
        else
            y = yElement;
    }
    return y;
}

void LArDiscreteCumulativeDistributionHelper::CreateDistributionFromCaloHits(const pandora::CaloHitList &caloHitList, const DiscreteCumulativeDistribution &distribution)
{
    caloHitList.sort(LArClusterHelper::SortHitsByPositionInX);
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        distribution.CollectInputData(pCaloHit->GetPositionVector().GetX(), pCaloHit->GetInputEnergy());
    }
    distribution.CreateCumulativeDistribution();
    return;
}

} // namespace lar_content
