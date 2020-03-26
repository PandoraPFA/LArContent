/**
 *  @file   larpandoracontent/LArHelpers/LArDiscreteCumulativeDistributionHelper.cc
 *
 *  @brief  Implementation of the discrete cumulative distribution helper class.
 *
 *  $Log: $
 */

#include <algorithm>
#include <iostream>

#include "larpandoracontent/LArHelpers/LArDiscreteCumulativeDistributionHelper.h"

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
    //TODO throw status code when distribution is empty
    float y(0.f);
    float xMin(0.f), xMax(0.f);
    distribution.GetXandY(0, xMin, y);
    distribution.GetXandY(distribution.GetSize()-1, xMax, y);

    if (x < xMin)
        return 0.f;

    if (x > xMax)
        return 1.;

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

void LArDiscreteCumulativeDistributionHelper::CreateDistributionFromCaloHits(pandora::CaloHitList caloHitList, DiscreteCumulativeDistribution &distribution)
{
    caloHitList.sort(LArClusterHelper::SortHitsByPositionInX);
    for (const pandora::CaloHit *const pCaloHit : caloHitList)
    {
        distribution.CollectInputData(pCaloHit->GetPositionVector().GetX(), pCaloHit->GetInputEnergy());
    }
    distribution.CreateCumulativeDistribution();
    return;
}

float LArDiscreteCumulativeDistributionHelper::CalculatePValueWithKSTestStatistic(const DiscreteCumulativeDistribution &distributionA, const DiscreteCumulativeDistribution &distributionB)
{   
    float ks(CalculateKSTestStatistic(distributionA, distributionB));

    float sum(CalculatePValueSumTerm(1, ks, distributionA, distributionB));
    float prev_sumTerm(std::numeric_limits<float>::epsilon());
    float sumTerm(std::numeric_limits<float>::max());

    int i(1);
    while ((abs(prev_sumTerm+sumTerm)/sum)>std::numeric_limits<float>::epsilon())
    {
        i++;
        prev_sumTerm = sumTerm;
	sumTerm = CalculatePValueSumTerm(i, ks, distributionA, distributionB);	
	sum += sumTerm;
    }

    float PValue(1-2*sum);

    return PValue;
}

float LArDiscreteCumulativeDistributionHelper::CalculatePValueSumTerm(const int i, const float &ks, const DiscreteCumulativeDistribution &distributionA, const DiscreteCumulativeDistribution &distributionB)
{
    float sumTerm(pow(-1,i-1)*exp(-2*pow(i*ks*sqrt(distributionA.GetSize()*distributionB.GetSize()/(distributionA.GetSize()+distributionB.GetSize())),2)));

    return sumTerm;
}
} // namespace lar_content
