/**
 *  @file   larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewTransverseTracksAlgorithm.cc
 *
 *  @brief  Implementation of the two view transverse tracks algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"


#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArDiscreteCumulativeDistributionHelper.h"


#include "larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewTransverseTracksAlgorithm.h"

#include "larpandoracontent/LArControlFlow/MultiPandoraApi.h"

using namespace pandora;

namespace lar_content
{

TwoViewTransverseTracksAlgorithm::TwoViewTransverseTracksAlgorithm() :
    m_nMaxMatrixToolRepeats(1000)
{

}

TwoViewTransverseTracksAlgorithm::~TwoViewTransverseTracksAlgorithm(){

}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewTransverseTracksAlgorithm::CalculateOverlapResult(const Cluster *const pCluster1, const Cluster *const pCluster2, const Cluster *const)
{
    std::cout << "=======================NEXTCOMPARISON======================" << std::endl;
    float xMin1(0.f), xMax1(0.f), xMin2(0.f), xMax2(0.f);
    LArClusterHelper::GetClusterSpanX(pCluster1, xMin1, xMax1);
    LArClusterHelper::GetClusterSpanX(pCluster2, xMin2, xMax2);

    const float xSpan1(xMax1 - xMin1), xSpan2(xMax2 - xMin2);

    if ((xSpan1 < std::numeric_limits<float>::epsilon()) || (xSpan2 < std::numeric_limits<float>::epsilon()))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const float xOverlapMin(std::max(xMin1, xMin2));
    const float xOverlapMax(std::min(xMax1, xMax2));
    const float xOverlap(xOverlapMax - xOverlapMin);
    TwoViewXOverlap twoViewXOverlap(xMin1, xMax1, xMin2, xMax2, xOverlap);

    float zMin1(0.f), zMax1(0.f);
    float zMin2(0.f), zMax2(0.f);
    LArClusterHelper::GetClusterSpanZ(pCluster1, xMin1, xMax1, zMin1, zMax1);
    LArClusterHelper::GetClusterSpanZ(pCluster2, xMin2, xMax2, zMin2, zMax2);
    CartesianVector boundingBoxMin1(xOverlapMin, 0.f, zMin1), boundingBoxMax1(xOverlapMax, 0.f, zMax1);
    CartesianVector boundingBoxMin2(xOverlapMin, 0.f, zMin2), boundingBoxMax2(xOverlapMax, 0.f, zMax2);
    pandora::CaloHitList overlapHits1, overlapHits2;
    LArClusterHelper::GetCaloHitListInBoundingBox(pCluster1, boundingBoxMin1, boundingBoxMax1, overlapHits1);
    LArClusterHelper::GetCaloHitListInBoundingBox(pCluster2, boundingBoxMin2, boundingBoxMax2, overlapHits2);

    DiscreteCumulativeDistribution disCumulDist1, disCumulDist2;
    pandora::CaloHitList hitList1, hitList2;
    LArClusterHelper::GetCaloHitList(pCluster1, hitList1);
    LArClusterHelper::GetCaloHitList(pCluster2, hitList2);

    LArDiscreteCumulativeDistributionHelper::CreateDistributionFromCaloHits(overlapHits1, disCumulDist1);
    LArDiscreteCumulativeDistributionHelper::CreateDistributionFromCaloHits(overlapHits2, disCumulDist2);

    if (0 == disCumulDist1.GetSize() || 0 == disCumulDist2.GetSize())
        return;

    Spline1f spline1(CreateSplineFromCumulativeDistribution(disCumulDist1));
    Spline1f spline2(CreateSplineFromCumulativeDistribution(disCumulDist2));

    DiscreteCumulativeDistribution resampledDisCumulDist1;
    DiscreteCumulativeDistribution resampledDisCumulDist2;

    int nSamplingPoints(std::min(disCumulDist1.GetSize(), disCumulDist2.GetSize())/5);
    nSamplingPoints = std::max(nSamplingPoints,10);
    for (float xPos = xOverlapMin; xPos < xOverlapMax; xPos += (xOverlapMax-xOverlapMin)/nSamplingPoints)
    {
        float q1 = spline1(xPos).coeff(0);
        resampledDisCumulDist1.CollectCumulativeData(xPos, q1);
        float q2 = spline2(xPos).coeff(0);
        resampledDisCumulDist2.CollectCumulativeData(xPos, q2);
    }

    ChargeProfile profile1(CreateProfileFromCumulativeDistribution(resampledDisCumulDist1));
    ChargeProfile profile2(CreateProfileFromCumulativeDistribution(resampledDisCumulDist2));

    size_t sizeWindowInBins(std::max(nSamplingPoints/3,10)); //ATTN - No real reason to have nSamplingPoints/3
    ScoreProfile xMatchingScore;
    float fracGoodScore(0.);

    if (profile1.size()>sizeWindowInBins)
    {
        xMatchingScore = SlidingWindowMatchingScore(sizeWindowInBins, profile1, profile2, fracGoodScore);
	       fracGoodScore /= (nSamplingPoints-sizeWindowInBins+1);
    }
    //float matchingScore(1-CalculateCorrelationCoefficient(profile1, profile2));
    float matchingScore(fracGoodScore);
    
    /*
    int NHits1(overlapHits1.size()), NHits2(overlapHits2.size());
    int k(2);
    int nSegments(1),nMaxSegments(7);
    while ( (std::min(NHits1,NHits2)/k)>25 && k<nMaxSegments+1 )
      {
        nSegments = k;
        k++;
      }

    std::vector<pandora::CaloHitList> segmentedOverlapHits1, segmentedOverlapHits2;

    LArDiscreteCumulativeDistributionHelper::SplitCaloHitList(nSegments, xOverlap, overlapHits1, segmentedOverlapHits1);
    LArDiscreteCumulativeDistributionHelper::SplitCaloHitList(nSegments, xOverlap, overlapHits2, segmentedOverlapHits2);

    for (int i = 0; i<nSegments; i++)
      {
        DiscreteCumulativeDistribution disCumulDist1, disCumulDist2;
        //LArDiscreteCumulativeDistributionHelper::CreateDistributionFromCaloHits(overlapHits1, disCumulDist1);                                                                                             
        //LArDiscreteCumulativeDistributionHelper::CreateDistributionFromCaloHits(overlapHits2, disCumulDist2);                                                                                             

	LArDiscreteCumulativeDistributionHelper::CreateDistributionFromCaloHits(segmentedOverlapHits1[i], disCumulDist1);
	LArDiscreteCumulativeDistributionHelper::CreateDistributionFromCaloHits(segmentedOverlapHits2[i], disCumulDist2);

        if (0 == disCumulDist1.GetSize() || 0 == disCumulDist2.GetSize())
          return;

        float matchingScoreKuip(LArDiscreteCumulativeDistributionHelper::CalculatePValueWithKuiperTestStatistic(disCumulDist1, disCumulDist2));
        float matchingScoreKS(LArDiscreteCumulativeDistributionHelper::CalculatePValueWithKSTestStatistic(disCumulDist1, disCumulDist2));
	std::cout << "MatchingScoreKUIPER: " << matchingScoreKuip << std::endl;
	std::cout << "MatchingScoreKS: " << matchingScoreKS << std::endl;
      }
    */
    TwoViewTransverseOverlapResult twoViewTransverseOverlapResult(matchingScore, twoViewXOverlap);

    if (xOverlap > std::numeric_limits<float>::epsilon())
        this->GetMatchingControl().GetOverlapMatrix().SetOverlapResult(pCluster1, pCluster2, twoViewTransverseOverlapResult);


    std::vector<float> x_U;
    std::vector<float> y_U;
    std::vector<float> x_V;
    std::vector<float> y_V;
    std::vector<float> y_prof_U;
    std::vector<float> y_prof_V;
    std::vector<float> score_prof;

//    for (size_t i = 0; i < disCumulDist1.GetSize(); i++){
//        float x,y;
//        disCumulDist1.GetXandY(i,x,y);
//        x_U.push_back(x);
//        y_U.push_back(y);
//    }
//    for (size_t i = 0; i < disCumulDist2.GetSize(); i++){
//        float x,y;
//        disCumulDist2.GetXandY(i,x,y);
//        x_V.push_back(x);
//        y_V.push_back(y);
//    }
    for (size_t i = 0; i < resampledDisCumulDist1.GetSize(); i++){
        float x,y;
        resampledDisCumulDist1.GetXandY(i,x,y);
        x_U.push_back(x);
        y_U.push_back(y);
        y_prof_U.push_back(profile1[i].second);
    }
    for (size_t i = 0; i < resampledDisCumulDist2.GetSize(); i++){
        float x,y;
        resampledDisCumulDist2.GetXandY(i,x,y);
        x_V.push_back(x);
        y_V.push_back(y);
        y_prof_V.push_back(profile2[i].second);
	score_prof.push_back(xMatchingScore[i].second);
    }

    const MCParticle* particle1(MCParticleHelper::GetMainMCParticle(pCluster1)); 
    const MCParticle* particle2(MCParticleHelper::GetMainMCParticle(pCluster2)); 
    int sameParticle(particle1->GetUid() == particle2->GetUid());
    int pdgU(particle1->GetParticleId());
    int pdgV(particle2->GetParticleId());
    int isPrimaryU(particle1->IsRootParticle());
    int isPrimaryV(particle2->IsRootParticle());

    int clusterSizeU = pCluster1->GetOrderedCaloHitList().size();
    int clusterSizeV = pCluster2->GetOrderedCaloHitList().size();
    int overlapNHitsU = disCumulDist1.GetSize();
    int overlapNHitsV = disCumulDist2.GetSize();
    int resampleOverlapNBinsU = profile1.size();
    int resampleOverlapNBinsV = profile2.size();

    float correlation(CalculateCorrelationCoefficient(profile1,profile2));
    float pvalue(CalculateTTestPValue(profile1,profile2));
    float xOverlapFractionU(twoViewXOverlap.GetXOverlapFractionU());
    float xOverlapFractionV(twoViewXOverlap.GetXOverlapFractionV());



    const pandora::Pandora * primary_pandora = MultiPandoraApi::GetPrimaryPandoraInstance(&(this->GetPandora()));

    PANDORA_MONITORING_API(SetTreeVariable(*primary_pandora, "matchtree", "sameparticle", sameParticle));
    PANDORA_MONITORING_API(SetTreeVariable(*primary_pandora, "matchtree", "pdgu", pdgU));
    PANDORA_MONITORING_API(SetTreeVariable(*primary_pandora, "matchtree", "pdgv", pdgV));
    PANDORA_MONITORING_API(SetTreeVariable(*primary_pandora, "matchtree", "isprimaryu", isPrimaryU));
    PANDORA_MONITORING_API(SetTreeVariable(*primary_pandora, "matchtree", "isprimaryv", isPrimaryV));

    PANDORA_MONITORING_API(SetTreeVariable(*primary_pandora, "matchtree", "clustersizeu", clusterSizeU));
    PANDORA_MONITORING_API(SetTreeVariable(*primary_pandora, "matchtree", "clustersizev", clusterSizeV));
    PANDORA_MONITORING_API(SetTreeVariable(*primary_pandora, "matchtree", "overlapnhitsu", overlapNHitsU));
    PANDORA_MONITORING_API(SetTreeVariable(*primary_pandora, "matchtree", "overlapnhitsv", overlapNHitsV));
    PANDORA_MONITORING_API(SetTreeVariable(*primary_pandora, "matchtree", "resampleoverlapnbinsu", resampleOverlapNBinsU));
    PANDORA_MONITORING_API(SetTreeVariable(*primary_pandora, "matchtree", "resampleoverlapnbinsv", resampleOverlapNBinsV));
    PANDORA_MONITORING_API(SetTreeVariable(*primary_pandora, "matchtree", "correlation", correlation));
    PANDORA_MONITORING_API(SetTreeVariable(*primary_pandora, "matchtree", "pvalue", pvalue));

    PANDORA_MONITORING_API(SetTreeVariable(*primary_pandora, "matchtree", "matchscore", matchingScore));
    PANDORA_MONITORING_API(SetTreeVariable(*primary_pandora, "matchtree", "xoverlap", xOverlap));
    PANDORA_MONITORING_API(SetTreeVariable(*primary_pandora, "matchtree", "xoverlapfracu", xOverlapFractionU));
    PANDORA_MONITORING_API(SetTreeVariable(*primary_pandora, "matchtree", "xoverlapfracv", xOverlapFractionV));


    //PANDORA_MONITORING_API(SetTreeVariable(*primary_pandora, "matchtree", "xu", &x_U));
    //PANDORA_MONITORING_API(SetTreeVariable(*primary_pandora, "matchtree", "yu", &y_U));
    //PANDORA_MONITORING_API(SetTreeVariable(*primary_pandora, "matchtree", "xv", &x_V));
    //PANDORA_MONITORING_API(SetTreeVariable(*primary_pandora, "matchtree", "yv", &y_V));
    //PANDORA_MONITORING_API(SetTreeVariable(*primary_pandora, "matchtree", "yprofu", &y_prof_U));
    //PANDORA_MONITORING_API(SetTreeVariable(*primary_pandora, "matchtree", "yprofv", &y_prof_V));
    //PANDORA_MONITORING_API(SetTreeVariable(*primary_pandora, "matchtree", "xmatchingscore", &score_prof));




    PANDORA_MONITORING_API(FillTree(*primary_pandora, "matchtree"));


    std::cout<<"Cluster 1 NHits: " << pCluster1->GetOrderedCaloHitList().size() << "  Cluster 2 NHits: " << pCluster2->GetOrderedCaloHitList().size() << std::endl;
    std::cout<<"KS PValue: " << LArDiscreteCumulativeDistributionHelper::CalculatePValueWithKSTestStatistic(resampledDisCumulDist1, resampledDisCumulDist2) << std::endl;
    std::cout<<"Kuiper PValue: " << LArDiscreteCumulativeDistributionHelper::CalculatePValueWithKuiperTestStatistic(resampledDisCumulDist1, resampledDisCumulDist2) << std::endl;
    std::cout<<"XOverlap: " << xOverlap << std::endl;
    std::cout<<"Correlation: " << correlation << std::endl;    
    std::cout<<"p-value: " << CalculateTTestPValue(profile1,profile2) << std::endl;
    std::cout<<"fracGoodScore = " << matchingScore << std::endl;
    std::cout<<"thanks"<<std::endl;



}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewTransverseTracksAlgorithm::ExamineOverlapContainer()
{
    const pandora::Pandora * primary_pandora = MultiPandoraApi::GetPrimaryPandoraInstance(&(this->GetPandora()));
    PANDORA_MONITORING_API(Create(*primary_pandora));
    PANDORA_MONITORING_API(SaveTree(*primary_pandora, "matchtree", "output.root", "RECREATE"));
    //PANDORA_MONITORING_API(Delete(*primary_pandora));


    unsigned int repeatCounter(0);

    for (MatrixToolVector::const_iterator iter = m_algorithmToolVector.begin(), iterEnd = m_algorithmToolVector.end(); iter != iterEnd; )
    {
        if ((*iter)->Run(this, this->GetMatchingControl().GetOverlapMatrix()))
        {
            iter = m_algorithmToolVector.begin();

            if (++repeatCounter > m_nMaxMatrixToolRepeats)
                break;
        }
        else
        {
            ++iter;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoViewTransverseTracksAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle,
        "TrackTools", algorithmToolVector));

    for (AlgorithmToolVector::const_iterator iter = algorithmToolVector.begin(), iterEnd = algorithmToolVector.end(); iter != iterEnd; ++iter)
    {
        TransverseMatrixTool *const pTransverseMatrixTool(dynamic_cast<TransverseMatrixTool*>(*iter));

        if (!pTransverseMatrixTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolVector.push_back(pTransverseMatrixTool);
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NMaxMatrixToolRepeats", m_nMaxMatrixToolRepeats));

    return BaseAlgorithm::ReadSettings(xmlHandle);
}


TwoViewTransverseTracksAlgorithm::Spline1f TwoViewTransverseTracksAlgorithm::CreateSplineFromCumulativeDistribution(const DiscreteCumulativeDistribution &cumulativeDistribution)
{
    Eigen::RowVectorXf xVals(cumulativeDistribution.GetSize());
    Eigen::RowVectorXf yVals(cumulativeDistribution.GetSize());

    float x(0.f), y(0.f);
    for (size_t iValue = 0; iValue < cumulativeDistribution.GetSize(); iValue++)
    {
        cumulativeDistribution.GetXandY(iValue, x, y);
        xVals(iValue) = x;
        yVals(iValue) = y;
    }
    Spline1f spline(Eigen::SplineFitting<Spline1f>::Interpolate(yVals, 1, xVals));
    return spline;
}

TwoViewTransverseTracksAlgorithm::ChargeProfile TwoViewTransverseTracksAlgorithm::CreateProfileFromCumulativeDistribution(const DiscreteCumulativeDistribution &cumulativeDistribution)
{
    ChargeProfile profile;
    for (size_t iValue = 0; iValue < cumulativeDistribution.GetSize(); iValue++)
    {
        float x(0.f), y(0.f);
        float prevY(0.f);
        if (iValue > 0) 
            cumulativeDistribution.GetXandY(iValue-1,x,prevY);
        cumulativeDistribution.GetXandY(iValue,x,y);
        profile.emplace_back(x,y-prevY);
    }
    return profile;
}

float TwoViewTransverseTracksAlgorithm::CalculateCorrelationCoefficient(const ChargeProfile &profile1, const ChargeProfile &profile2)
{
    float covariance(0.f);
    float variance1(0.f);
    float variance2(0.f);

    float yAv1(0.f), yAv2(0.f);
    for (size_t iElement = 0; iElement < profile1.size(); iElement++)
    {
        yAv1+=profile1[iElement].second;
        yAv2+=profile2[iElement].second;
    }
    if (profile1.size() > 0)
    {
        yAv1 /= profile1.size();
        yAv2 /= profile2.size();
    }

    for (size_t iElement = 0; iElement < profile1.size(); iElement++)
    {
        //float x(profile1[iElement].first);
        float y1(profile1[iElement].second);
        float y2(profile2[iElement].second);

        covariance += (y1-yAv1)*(y2-yAv2);
        variance1 += (y1-yAv1)*(y1-yAv1);
        variance2 += (y2-yAv2)*(y2-yAv2);
    }
    if (profile1.size() > 1)
    {
        covariance /= (profile1.size()-1);
        variance1 /= (profile1.size()-1);
        variance2 /= (profile2.size()-1);
    }

    float correlation(0.f);
    if (variance1*variance2 > 0)
        correlation = covariance/=(sqrt(variance1*variance2));

    return correlation;
}

float TwoViewTransverseTracksAlgorithm::CalculateTTestValue(const float x, const float coefficient, const float dof)
{
    return coefficient*std::pow( 1.0 + x*x/dof, -0.5 *(dof + 1.0));
}

float TwoViewTransverseTracksAlgorithm::CalculateTTestPValue(const ChargeProfile &profile1, const ChargeProfile &profile2)
{
    if (profile1.size() != profile2.size())
        return -999.;
    float correlation(CalculateCorrelationCoefficient(profile1,profile2));
    float dof(profile1.size() - 2);
    float tTestStatistic(correlation*sqrt(dof)/(sqrt(1.-correlation*correlation)));
    float tDistCoeff(std::tgamma(0.5 * (dof+1.)) / std::tgamma(0.5*dof) / (std::sqrt(dof*3.14159265359)));

    int nSteps(10000);
    float upperLimit(15.f);
    float dx((upperLimit-tTestStatistic)/nSteps);
    float integral(CalculateTTestValue(tTestStatistic, tDistCoeff, dof) + CalculateTTestValue(upperLimit, tDistCoeff, dof));
    for (int iStep = 1; iStep < nSteps; iStep++)
        integral+=2. * CalculateTTestValue(tTestStatistic + iStep*dx, tDistCoeff, dof);
    integral *= dx/2.0;


    return integral;
}

  TwoViewTransverseTracksAlgorithm::ScoreProfile TwoViewTransverseTracksAlgorithm::SlidingWindowMatchingScore(const size_t &sizeWindowInBins, const ChargeProfile &profile1, const ChargeProfile &profile2, float &fracGoodScore)
{
  ScoreProfile xMatchingScore;
  for (size_t k = 0; k<profile1.size(); k++)
    {
      if ( (k<(profile1.size()-sizeWindowInBins/2)) && (k>(sizeWindowInBins/2-1)) )
	{
      ChargeProfile windowProfile1(GetWindow(k, sizeWindowInBins, profile1));
      ChargeProfile windowProfile2(GetWindow(k, sizeWindowInBins, profile2));
      xMatchingScore.emplace_back(profile1[k].first, 1-CalculateTTestPValue(windowProfile1, windowProfile2));
      if (1-CalculateTTestPValue(windowProfile1, windowProfile2)>0.99999) fracGoodScore++; //ATTN - Arbitrary 0.99999 value of matching score
	}
      else
	{
	  xMatchingScore.emplace_back(profile1[k].first, 0.5);
	}
    }
  return xMatchingScore;
}

TwoViewTransverseTracksAlgorithm::ChargeProfile   TwoViewTransverseTracksAlgorithm::GetWindow(size_t &i, const size_t &sizeWindowInBins, const ChargeProfile &profile)
{
  ChargeProfile windowedProfile;
  for (size_t iElement = (i-sizeWindowInBins/2); iElement < (i+sizeWindowInBins/2 + 1); iElement++)
    {
      windowedProfile.emplace_back(profile[iElement].first, profile[iElement].second);
    }
  return windowedProfile;

}
  TwoViewTransverseTracksAlgorithm::ScoreProfile TwoViewTransverseTracksAlgorithm::SlidingWindowMatchingScore(const size_t &sizeWindowInBins, const DiscreteCumulativeDistribution &disCumulDist1, const DiscreteCumulativeDistribution &disCumulDist2, float &fracGoodScore)
{
  float x(0), y(0);
  ScoreProfile xMatchingScore;
  for (size_t iElement = 0; iElement<disCumulDist1.GetSize(); iElement++)
    {
      disCumulDist1.GetXandY(iElement, x, y);
      if ( (iElement<(disCumulDist1.GetSize()-sizeWindowInBins/2)) && (iElement>(sizeWindowInBins/2-1)) )
	{
      ChargeProfile windowProfile1(GetWindow(iElement, sizeWindowInBins, disCumulDist1));
      ChargeProfile windowProfile2(GetWindow(iElement, sizeWindowInBins, disCumulDist2));
      xMatchingScore.emplace_back(x, LArDiscreteCumulativeDistributionHelper::CalculatePValueWithKSTestStatistic(disCumulDist1, disCumulDist2));
      //std::cout<<LArDiscreteCumulativeDistributionHelper::CalculatePValueWithKSTestStatistic(disCumulDist1, disCumulDist2)<<std::endl;
      if (LArDiscreteCumulativeDistributionHelper::CalculatePValueWithKSTestStatistic(disCumulDist1, disCumulDist2)>0.9) fracGoodScore++; //ATTN - Arbitrary 0.6 value of matching score
	}
      else
	{
	  xMatchingScore.emplace_back(x, 0.5);
	}
    }
  return xMatchingScore;
}

  TwoViewTransverseTracksAlgorithm::ChargeProfile  TwoViewTransverseTracksAlgorithm::GetWindow(size_t &i, const size_t &sizeWindowInBins, const DiscreteCumulativeDistribution &disCumulDist)
{
  float x(0), y(0);
  ChargeProfile windowedProfile; //ATTN - Should use another type?
  for (size_t iElement = (i-sizeWindowInBins/2); iElement < (i+sizeWindowInBins/2 + 1); iElement++)
    {
      disCumulDist.GetXandY(iElement, x, y);
      windowedProfile.emplace_back(x, y);
    }
  return windowedProfile;

}

} // namespace lar_content
