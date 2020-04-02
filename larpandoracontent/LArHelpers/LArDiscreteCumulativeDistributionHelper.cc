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

  void LArDiscreteCumulativeDistributionHelper::SplitCaloHitList(const int &nSegments, const float &xOverlap, const pandora::CaloHitList &overlapHits, std::vector<pandora::CaloHitList> &segmentedOverlapHits)
  {
    float splitXOverlap(xOverlap/nSegments);

    for (int k = 0; k<nSegments; k++)
      {
	pandora::CaloHitList caloHitList;
	for (const pandora::CaloHit *const pCaloHit : overlapHits)
	  {
	    if (pCaloHit->GetPositionVector().GetX()>k*splitXOverlap && pCaloHit->GetPositionVector().GetX()<(k+1)*splitXOverlap)
	      {
	      caloHitList.push_back(pCaloHit);
	      }

	    //if (segmentedOverlapHits[k].size()<25 && k>0)                                                                                                
	    //segmentedOverlapHits[k-1].insert(segmentedOverlapHits[k-1].end(), segmentedOverlapHits[k].begin(), segmentedOverlapHits[k].end());           
	  }
	segmentedOverlapHits.push_back(caloHitList);
	std::cout << segmentedOverlapHits.size() << std::endl;
      }
    return;
  }

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

    std::cout << "ks = " << ks << std::endl;
    return ks;
}
float LArDiscreteCumulativeDistributionHelper::CalculateKuiperTestStatistic(const DiscreteCumulativeDistribution &distributionA, const DiscreteCumulativeDistribution &distributionB)
{

    float Dp = std::numeric_limits<float>::min();
    float Dm = std::numeric_limits<float>::min();
    
    for (size_t iElement = 0; iElement < distributionA.GetSize(); ++iElement)
      {
        float xA(0), yA(0);
        distributionA.GetXandY(iElement,xA,yA);
        float yB = FindY(distributionB, xA);
        Dp = std::max(Dp, yA-yB);
        Dm = std::max(Dm, yB-yA);
      }
    for (size_t iElement = 0; iElement < distributionB.GetSize(); ++iElement)
      {
        float xB(0), yB(0);
        distributionB.GetXandY(iElement,xB,yB);
        float yA = FindY(distributionA, xB);
        Dp = std::max(Dp, yA-yB);
        Dm = std::max(Dm, yB-yA);
      }


    float kuip (Dp+Dm);

    return kuip;
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

    float sum(CalculatePValueSumKSTerm(1, ks, distributionA, distributionB));
    float prev_sumTerm(std::numeric_limits<float>::epsilon());
    float sumTerm(999);

    int i(1);
    while (fabs((prev_sumTerm+sumTerm)/sum)>std::numeric_limits<float>::epsilon())
    {
        i++;
        prev_sumTerm = sumTerm;
	sumTerm = CalculatePValueSumKSTerm(i, ks, distributionA, distributionB);	
	sum += sumTerm;
    }
    
    float PValue(2*sum);
    std::cout << "PValue = " << PValue << std::endl;
    return PValue;
}

float LArDiscreteCumulativeDistributionHelper::CalculatePValueSumKSTerm(const int i, const float &ks, const DiscreteCumulativeDistribution &distributionA, const DiscreteCumulativeDistribution &distributionB)
{
    float sumTerm(pow(-1,i-1)*exp(-2*pow(i*ks*sqrt(distributionA.GetSize()*distributionB.GetSize()/(distributionA.GetSize()+distributionB.GetSize())),2)));

    return sumTerm;
}
  float LArDiscreteCumulativeDistributionHelper::CalculatePValueWithKuiperTestStatistic(const DiscreteCumulativeDistribution &distributionA, const DiscreteCumulativeDistribution &distributionB)
  {
    float kuip(CalculateKuiperTestStatistic(distributionA, distributionB));

    float sum(CalculatePValueSumKuiperTerm(1, kuip, distributionA, distributionB));
    float prev_sumTerm(std::numeric_limits<float>::epsilon());
    float sumTerm(999);

    int i(1);
    while (fabs((prev_sumTerm+sumTerm)/sum)>std::numeric_limits<float>::epsilon())
      {
        i++;
        prev_sumTerm = sumTerm;
        sumTerm = CalculatePValueSumKuiperTerm(i, kuip, distributionA, distributionB);
        sum += sumTerm;
      }

    float PValue(2*sum);
    std::cout << "PValue = " << PValue << std::endl;
    return PValue;
  }

  float LArDiscreteCumulativeDistributionHelper::CalculatePValueSumKuiperTerm(const int i, const float &kuip, const DiscreteCumulativeDistribution &distributionA, const DiscreteCumulativeDistribution &distributionB)
  {
    float sumTerm((4*pow(i*kuip*sqrt(distributionA.GetSize()*distributionB.GetSize()/(distributionA.GetSize()+distributionB.GetSize())),2)-1)*exp(-2*pow(i*kuip*sqrt(distributionA.GetSize()*distributionB.GetSize()/(distributionA.GetSize()+distributionB.GetSize())),2)));

    return sumTerm;
  }

} // namespace lar_content
