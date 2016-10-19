/**
 *  @file   larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/TracksCrossingGapsTool.cc
 * 
 *  @brief  Implementation of the long tracks tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/LongTracksTool.h"
#include "larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/TracksCrossingGapsTool.h"


using namespace pandora;

namespace lar_content
{

TracksCrossingGapsTool::TracksCrossingGapsTool() ://notice the change in matched points/fraction as in CrossGapsAssociationAlgorithm
    m_minMatchedFraction(0.5f),
    m_minMatchedSamplingPoints(10),
    m_minXOverlapFraction(0.9f),
    m_minMatchedSamplingPointRatio(2),
    m_maxGapTolerance(2.f),
    m_sampleStepSize(0.5f),
    m_maxAngleRatio(2)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TracksCrossingGapsTool::Run(ThreeDTransverseTracksAlgorithm *const pAlgorithm, TensorType &overlapTensor)
{
  if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
    std::cout << "----> Running Algorithm Tool: " << this << ", " << this->GetType() << std::endl;
  
  ProtoParticleVector protoParticleVector;
  this->FindTracks(pAlgorithm, overlapTensor, protoParticleVector);
  
  const bool particlesMade(pAlgorithm->CreateThreeDParticles(protoParticleVector));
  return particlesMade;
}
  
//------------------------------------------------------------------------------------------------------------------------------------------

void TracksCrossingGapsTool::FindTracks(ThreeDTransverseTracksAlgorithm *const pAlgorithm, const TensorType &overlapTensor, ProtoParticleVector &protoParticleVector) const
{

  ClusterSet usedClusters;
  ClusterVector sortedKeyClusters;
  overlapTensor.GetSortedKeyClusters(sortedKeyClusters);

    for (const Cluster *const pKeyCluster : sortedKeyClusters)
    {
        if (!pKeyCluster->IsAvailable())
            continue;

        unsigned int nU(0), nV(0), nW(0);
        TensorType::ElementList elementList;
        overlapTensor.GetConnectedElements(pKeyCluster, true, elementList, nU, nV, nW);

        IteratorList iteratorList;
        this->SelectElements(pAlgorithm, elementList, usedClusters, iteratorList);

        // Check that elements are not directly connected and are significantly longer than any other directly connected elements
        for (IteratorList::const_iterator iIter = iteratorList.begin(), iIterEnd = iteratorList.end(); iIter != iIterEnd; ++iIter)
        {
            if (LongTracksTool::HasLongDirectConnections(iIter, iteratorList))
                continue;

            if (!LongTracksTool::IsLongerThanDirectConnections(iIter, elementList, m_minMatchedSamplingPointRatio, usedClusters))
                continue;

            ProtoParticle protoParticle;
            protoParticle.m_clusterListU.push_back((*iIter)->GetClusterU());
            protoParticle.m_clusterListV.push_back((*iIter)->GetClusterV());
            protoParticle.m_clusterListW.push_back((*iIter)->GetClusterW());
            protoParticleVector.push_back(protoParticle);

            usedClusters.insert((*iIter)->GetClusterU());
            usedClusters.insert((*iIter)->GetClusterV());
            usedClusters.insert((*iIter)->GetClusterW());
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TracksCrossingGapsTool::SelectElements(ThreeDTransverseTracksAlgorithm *const pAlgorithm, const TensorType::ElementList &elementList, 
					    const pandora::ClusterSet &usedClusters, IteratorList &iteratorList) const
{
  for (TensorType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
    {
      if (usedClusters.count(eIter->GetClusterU()) || usedClusters.count(eIter->GetClusterV()) || usedClusters.count(eIter->GetClusterW()))
	continue;
      
      if (eIter->GetOverlapResult().GetMatchedFraction() < m_minMatchedFraction)
	continue;
      
      if (eIter->GetOverlapResult().GetNMatchedSamplingPoints() < m_minMatchedSamplingPoints)
	continue;
      
      const XOverlap &xOverlap(eIter->GetOverlapResult().GetXOverlap());

      //LORENA
      /*std::cout << " --- SelectElements --- " << std::endl;
      std::cout << "Overlap Span according to tensor: " << xOverlap.GetXOverlapSpan() << std::endl;
      std::cout << "  -fraction in U from tensor: " <<  xOverlap.GetXOverlapSpan()/xOverlap.GetXSpanU() << std::endl;
      std::cout << "  -fraction in V from tensor: " <<  xOverlap.GetXOverlapSpan()/xOverlap.GetXSpanV() << std::endl;
      std::cout << "  -fraction in W from tensor: " <<  xOverlap.GetXOverlapSpan()/xOverlap.GetXSpanW() << std::endl;*/

      if (xOverlap.GetXOverlapSpan() < std::numeric_limits<float>::epsilon())
	continue;
      
      //Calculate effective overlap fraction, including the information about gaps
      float xOverlapFractionU(0.f), xOverlapFractionV(0.f), xOverlapFractionW(0.f);
      PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,this->CalculateEffectiveOverlapFractions(pAlgorithm, *eIter,xOverlapFractionU,xOverlapFractionV,xOverlapFractionW));
      //LORENA
      /*std::cout << "Effective overlap Span: " << std::endl;
      std::cout << "  - in U: " << xOverlapFractionU << std::endl;
      std::cout << "  - in V: " << xOverlapFractionV << std::endl;
      std::cout << "  - in W: " << xOverlapFractionW << std::endl;*/

      if ((xOverlap.GetXSpanU() > std::numeric_limits<float>::epsilon()) && (xOverlapFractionU > m_minXOverlapFraction) &&
	  (xOverlap.GetXSpanV() > std::numeric_limits<float>::epsilon()) && (xOverlapFractionV > m_minXOverlapFraction) &&
	  (xOverlap.GetXSpanW() > std::numeric_limits<float>::epsilon()) && (xOverlapFractionW > m_minXOverlapFraction))
        {
	  iteratorList.push_back(eIter);
        }
      
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------     
StatusCode TracksCrossingGapsTool::CalculateEffectiveOverlapFractions(ThreeDTransverseTracksAlgorithm *const pAlgorithm, const TensorType::Element &element, float &xOverlapFractionU, float &xOverlapFractionV, float &xOverlapFractionW) const
{
  
  //get the three clusters for the element                                                                                                               
  const Cluster *const clusterU = element.GetClusterU();                                                                                          
  const Cluster *const clusterV = element.GetClusterV();                                                                                               
  const Cluster *const clusterW = element.GetClusterW(); 

  //TODO - another way to catch exception here?
  if(!clusterU || !clusterV || !clusterW)
    return STATUS_CODE_NOT_FOUND;
  
  float xMinEffU(-std::numeric_limits<float>::max()), xMinEffV(-std::numeric_limits<float>::max()), xMinEffW(-std::numeric_limits<float>::max());  
  float xMaxEffU(+std::numeric_limits<float>::max()), xMaxEffV(+std::numeric_limits<float>::max()), xMaxEffW(+std::numeric_limits<float>::max());
  
  LArClusterHelper::GetClusterSpanX(clusterU, xMinEffU, xMaxEffU);
  LArClusterHelper::GetClusterSpanX(clusterV, xMinEffV, xMaxEffV);
  LArClusterHelper::GetClusterSpanX(clusterW, xMinEffW, xMaxEffW);


  PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, this->CalculateEffectiveOverlapSpan(pAlgorithm,clusterU,xMinEffU,xMaxEffU,clusterV,xMinEffV,xMaxEffV,clusterW,xMinEffW,xMaxEffW));

  const float effXSpanU(xMaxEffU-xMinEffU), effXSpanV(xMaxEffV-xMinEffV), effXSpanW(xMaxEffW-xMinEffW);
  const float xMinEff(std::max(xMinEffU, std::max(xMinEffV,xMinEffW)));
  const float xMaxEff(std::min(xMaxEffU, std::min(xMaxEffV,xMaxEffW)));
  float effXOverlapSpan(xMaxEff-xMinEff); //start from original overlap span from tensor
  
  //  std::cout << "effXOverlapSpan=" << effXOverlapSpan << std::endl;
  //LORENA
  const CartesianVector startPointEff(xMinEff,0,-1000);
  const CartesianVector endPointEff(xMinEff,0,1000);
  const CartesianVector startPointEffmax(xMaxEff,0,-1000);
  const CartesianVector endPointEffmax(xMaxEff,0,1000);

  PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &startPointEff, &endPointEff, "MinEff", BLACK, 1,2));
  PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &startPointEffmax, &endPointEffmax, "MaxEff", BLACK, 1,1));
  PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

  //Note: the XSpan in each view is now updated with the fraction in gaps, in principle should be less than 1
  //but xOverlap.GetXSpanView() sometimes gives a span smaller than the difference xMax-xMin
  xOverlapFractionU = std::min(1.f, (effXOverlapSpan / effXSpanU));
  xOverlapFractionV = std::min(1.f, (effXOverlapSpan / effXSpanV)); 
  xOverlapFractionW = std::min(1.f, (effXOverlapSpan / effXSpanW));
  
  return STATUS_CODE_SUCCESS;
 
}

//------------------------------------------------------------------------------------------------------------------------------------------     
StatusCode TracksCrossingGapsTool::CalculateEffectiveOverlapSpan(ThreeDTransverseTracksAlgorithm *const pAlgorithm,const pandora::Cluster *const pClusterU, float &xMinEffU, float &xMaxEffU, const pandora::Cluster *const pClusterV, float &xMinEffV, float &xMaxEffV, const pandora::Cluster *const pClusterW, float &xMinEffW, float &xMaxEffW) const
{

  const float xMin(std::min(xMinEffU, std::min(xMinEffV,xMinEffW)));
  const float xMax(std::max(xMaxEffU, std::max(xMaxEffV,xMaxEffW)));
  float xMinEff(std::max(xMinEffU, std::max(xMinEffV,xMinEffW)));
  float xMaxEff(std::min(xMaxEffU, std::min(xMaxEffV,xMaxEffW)));

  float dxUmin(0.f), dxVmin(0.f), dxWmin(0.f);
  float dxUmax(0.f), dxVmax(0.f), dxWmax(0.f);

  //LORENA
  /*std::cout << "xMinU = " << xMinU << ", xMinV= " << xMinV << ", xMinW = " << xMinW << std::endl;
  std::cout << "xMaxU = " << xMaxU << ", xMaxV= " << xMaxV << ", xMaxW = " << xMaxW << std::endl;

  std::cout << "xMin = " << xMin << std::endl;
  std::cout << "xMax = " << xMax << std::endl;

  std::cout << "xMinEff = " << xMinEff << std::endl;
  std::cout << "xMaxEff = " << xMaxEff << std::endl;


    const CartesianVector startPointU(xMinU,0,-1000);
  const CartesianVector endPointU(xMinU,0,1000);
  const CartesianVector startPointUmax(xMaxU,0,-1000);
  const CartesianVector endPointUmax(xMaxU,0,1000);

  PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &startPointU, &endPointU, "MinU", RED, 1,2));
  PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &startPointUmax, &endPointUmax, "MaxU", RED, 1,1));

 const CartesianVector startPointV(xMinV,0,-1000);
  const CartesianVector endPointV(xMinV,0,1000);
  const CartesianVector startPointVmax(xMaxV,0,-1000);
  const CartesianVector endPointVmax(xMaxV,0,1000);

  PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &startPointV, &endPointV, "MinV", GREEN, 1,2));
  PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &startPointVmax, &endPointVmax, "MaxV", GREEN, 1,1));

 const CartesianVector startPointW(xMinW,0,-1000);
  const CartesianVector endPointW(xMinW,0,1000);
  const CartesianVector startPointWmax(xMaxW,0,-1000);
  const CartesianVector endPointWmax(xMaxW,0,1000);

  PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &startPointW, &endPointW, "MinW", BLUE, 1,2));
  PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &startPointWmax, &endPointWmax, "MaxW", BLUE, 1,1));*/
  

  //find layers and sampling points
  const unsigned int nSamplingPointsLeft(1+std::round((xMinEff-xMin)/ m_sampleStepSize));//TODO alternative for round?
  const unsigned int nSamplingPointsRight(1+std::round((xMax-xMaxEff)/ m_sampleStepSize));

  //LORENA  
  //std::cout << "nSamplingPointsLeft=" << nSamplingPointsLeft << std::endl;
  //std::cout << "nSamplingPointsRight=" << nSamplingPointsRight << std::endl;

  //visual debugging
  //LORENA
  /*  ClusterList clusterListU, clusterListV, clusterListW;                                                                                                                 
  clusterListU.push_back(pClusterU);                                                                                                                                          
  clusterListV.push_back(pClusterV);                                                                                                                                          
  clusterListW.push_back(pClusterW);                                                                                                                                          
  PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));                                                            
  PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clusterListU, "UCluster", RED));                                                      
  PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clusterListV, "VCluster", GREEN));                                                    
  PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clusterListW, "WCluster", BLUE));   
  PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));*/

  for (unsigned int iSample = 1; iSample <= nSamplingPointsLeft; ++iSample)
    {
      bool gapInU(false), gapInV(false), gapInW(false); //to update each cluster span
      const float xSample(std::max(xMin,xMinEff - iSample*m_sampleStepSize));//protection not to go further than limits
      //      std::cout << "xSample = " << xSample << std::endl;//LORENA
      if(!this->PassesGapsChecks(pAlgorithm,xSample,pClusterU,pClusterV,pClusterW,gapInU,gapInV,gapInW))
	break;
      //Note: break instead of continue here to avoid finding a non-related gap far from the cluster itself if the x range is long
      //if a point is close to gap but not yet there should not be a problem with a generous gap tolerance 
      if(gapInU) dxUmin = xMinEffU-xSample;
      if(gapInV) dxVmin = xMinEffV-xSample;
      if(gapInW) dxWmin = xMinEffW-xSample;
    }

  for (unsigned int iSample = 1; iSample <= nSamplingPointsRight; ++iSample)
    {
      bool gapInU(false), gapInV(false), gapInW(false); //to update each cluster span
      const float xSample(std::min(xMax,xMaxEff + iSample*m_sampleStepSize));
      //      std::cout << "xSample = " << xSample << std::endl;//LORENA
      if(!this->PassesGapsChecks(pAlgorithm,xSample,pClusterU,pClusterV,pClusterW,gapInU,gapInV,gapInW))
	break;
      if(gapInU) dxUmax = xSample-xMaxEffU;
      if(gapInV) dxVmax = xSample-xMaxEffV;
      if(gapInW) dxWmax = xSample-xMaxEffW;
    }

  //update span in each view
  xMinEffU -=dxUmin;
  xMaxEffU +=dxUmax;
  xMinEffV -=dxVmin;
  xMaxEffV +=dxVmax;
  xMinEffW -=dxWmin;
  xMaxEffW +=dxWmax;

  /*  std::cout << "dxMin = " << dxMin << " dxMax = " << dxMax << std::endl;
  std::cout << "dxUMin = " << dxUmin << " dxUMax = " << dxUmax << std::endl;
  std::cout << "dxVMin = " << dxVmin << " dxVMax = " << dxVmax << std::endl;
  std::cout << "dxWMin = " << dxWmin << " dxWMax = " << dxWmax << std::endl;*/

  //LORENA
  /*  std::cout << "Effective overlap span between: " << xMinEff << " and " << xMaxEff << std::endl;
  std::cout << "Effective  span in U between: " << (xMaxU+dxUmax) << " and " << (xMinU-dxUmin) << " and " << effXSpanU+dxUmax+dxUmin << std::endl;
  std::cout << "Effective  span in V between: " << (xMaxV+dxVmax) << " and " << (xMinV-dxVmin) << std::endl;
  std::cout << "Effective  span in W between: " << (xMaxW+dxWmax) << " and " << (xMinW-dxWmin) << std::endl;*/

  /*    std::cout << "Effective  span ratio in U  " << effXOverlapSpan/effXSpanU <<std::endl;
  std::cout << "Effective  span ratio in V  " << effXOverlapSpan/effXSpanV <<std::endl;
  std::cout << "Effective  span ratio in W  " << effXOverlapSpan/effXSpanW <<std::endl;*/



  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------     
bool TracksCrossingGapsTool::PassesGapsChecks(ThreeDTransverseTracksAlgorithm *const pAlgorithm, const float &xSample, const pandora::Cluster *const pClusterU,const pandora::Cluster *const pClusterV, const pandora::Cluster *const pClusterW, bool &gapInU, bool &gapInV, bool &gapInW) const 
{

  const TwoDSlidingFitResult &slidingFitResultU(pAlgorithm->GetCachedSlidingFitResult(pClusterU));                                                 
  const TwoDSlidingFitResult &slidingFitResultV(pAlgorithm->GetCachedSlidingFitResult(pClusterV));                                                 
  const TwoDSlidingFitResult &slidingFitResultW(pAlgorithm->GetCachedSlidingFitResult(pClusterW));                                                 
  
  CartesianVector fitUPosition(0.f, 0.f, 0.f), fitVPosition(0.f, 0.f, 0.f), fitWPosition(0.f, 0.f, 0.f);
  //Check first we don't have the global x position in the three clusters, if we do, then there are no gaps involved 
  if((STATUS_CODE_SUCCESS == slidingFitResultU.GetGlobalFitPositionAtX(xSample,fitUPosition)) &&
     (STATUS_CODE_SUCCESS == slidingFitResultV.GetGlobalFitPositionAtX(xSample,fitVPosition)) &&
     (STATUS_CODE_SUCCESS == slidingFitResultW.GetGlobalFitPositionAtX(xSample,fitWPosition)))
    return false;

  //Note: order is important in function below - first one is assumed to be the one with gap initially
  //however, once inside CheckXPositionInGap we check for gaps also in the other two views
  //this is a way to avoid repeating this three times
  else 
    {
      if((STATUS_CODE_SUCCESS != slidingFitResultU.GetGlobalFitPositionAtX(xSample,fitUPosition)) && (!this->IsEndOfCluster(xSample, pClusterU, slidingFitResultU)))
	return(this->CheckXPositionInGap(xSample,pClusterU,slidingFitResultU,pClusterV,slidingFitResultV,pClusterW,slidingFitResultW, gapInU,gapInV,gapInW));
      else if((STATUS_CODE_SUCCESS != slidingFitResultV.GetGlobalFitPositionAtX(xSample,fitVPosition)) && (!this->IsEndOfCluster(xSample, pClusterV, slidingFitResultV)))
	return(this->CheckXPositionInGap(xSample,pClusterV,slidingFitResultV,pClusterU,slidingFitResultU,pClusterW,slidingFitResultW, gapInV,gapInU,gapInW));
      else 
	return(this->CheckXPositionInGap(xSample,pClusterW,slidingFitResultW,pClusterU,slidingFitResultU,pClusterV,slidingFitResultV, gapInW,gapInU,gapInV));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------     
bool TracksCrossingGapsTool::CheckXPositionInGap(const float &xSample, const pandora::Cluster *const pClusterGap, const TwoDSlidingFitResult &slidingFitResultGap, const pandora::Cluster *const pSecondCluster, const TwoDSlidingFitResult &secondSlidingFitResult, const pandora::Cluster *const pThirdCluster, const TwoDSlidingFitResult &thirdSlidingFitResult,bool &gapInFirst, bool &gapInSecond, bool &gapInThird)const
{
  CartesianVector secondFitPosition(0.f, 0.f, 0.f), thirdFitPosition(0.f, 0.f, 0.f);
  CartesianVector secondFitDirection(0.f, 0.f, 0.f), thirdFitDirection(0.f, 0.f, 0.f);

  //If we have the global position at X from the two other clusters, calculate projection in the first view and check for gaps
  if((STATUS_CODE_SUCCESS == secondSlidingFitResult.GetGlobalFitPositionAtX(xSample,secondFitPosition)) &&
     (STATUS_CODE_SUCCESS == thirdSlidingFitResult.GetGlobalFitPositionAtX(xSample,thirdFitPosition)))
    {
      const float second(secondFitPosition.GetZ()), third(thirdFitPosition.GetZ());
      const float first(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), LArClusterHelper::GetClusterHitType(pSecondCluster),  LArClusterHelper::GetClusterHitType(pThirdCluster), second, third));
      
      CartesianVector innerCoordinate(0.f, 0.f, 0.f), outerCoordinate(0.f, 0.f, 0.f), startCoordinate(0.f, 0.f, 0.f);
      LArClusterHelper::GetExtremalCoordinates(pClusterGap, innerCoordinate, outerCoordinate);
      startCoordinate = ((std::fabs(innerCoordinate.GetX()-xSample) < std::fabs(outerCoordinate.GetX()-xSample))? innerCoordinate : outerCoordinate);
      const CartesianVector samplingPoint(xSample,startCoordinate.GetY(),first); //probably unnecesary, just set Y = 0?

      PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &samplingPoint, "SamplingPoint", ORANGE, 1));//LORENA
      gapInFirst = LArGeometryHelper::IsInGap(this->GetPandora(), samplingPoint, LArClusterHelper::GetClusterHitType(pClusterGap), m_maxGapTolerance);
      return(gapInFirst);
    }

  //if we dont have a projection at x in the other two clusters, check if they are in gaps or at the end of the cluster
  else if((STATUS_CODE_SUCCESS != secondSlidingFitResult.GetGlobalFitPositionAtX(xSample,secondFitPosition)) &&
	  (STATUS_CODE_SUCCESS != thirdSlidingFitResult.GetGlobalFitPositionAtX(xSample,thirdFitPosition)))
    {
      const bool endSecond(this->IsEndOfCluster(xSample, pSecondCluster, secondSlidingFitResult));
      const bool endThird(this->IsEndOfCluster(xSample, pThirdCluster, thirdSlidingFitResult));      
      
      //now we must use a different method, using the direction at the extreme coordinates of the cluster, as there is no uv2w, etc
      gapInFirst = this->CheckGaps(xSample, pClusterGap,slidingFitResultGap);
      if(!endSecond)
	gapInSecond = this->CheckGaps(xSample, pSecondCluster,secondSlidingFitResult);  
      if(!endThird)
	gapInThird = this->CheckGaps(xSample, pThirdCluster,thirdSlidingFitResult);
      
      //now to pass the gaps check it has to be true in the first cluster and one of the other two clusters
      //      std::cout << "gapInFirst = " << gapInFirst << " , gapInSecond = " << gapInSecond << " , gapInThird = " << gapInThird << ", endSecond = " << endSecond << " , endThird + " << endThird << std::endl;
      return((gapInFirst && gapInSecond && endThird) ||
	     (gapInFirst && gapInThird && endSecond) ||
	     (gapInFirst && endSecond && endThird));
    }
  
  //Now check whether there is a second gap involved
  else if (STATUS_CODE_SUCCESS != secondSlidingFitResult.GetGlobalFitPositionAtX(xSample,secondFitPosition))
    {
      const bool endSecond(this->IsEndOfCluster(xSample, pSecondCluster, secondSlidingFitResult));
      gapInFirst = this->CheckGaps(xSample, pClusterGap,slidingFitResultGap);
      gapInSecond = this->CheckGaps(xSample, pSecondCluster,secondSlidingFitResult);
      return ((gapInFirst && gapInSecond) || (gapInFirst && endSecond));
    }
  else
    {
      const bool endThird(this->IsEndOfCluster(xSample, pThirdCluster, thirdSlidingFitResult));       
      gapInFirst = this->CheckGaps(xSample, pClusterGap,slidingFitResultGap);
      gapInThird = this->CheckGaps(xSample, pThirdCluster,thirdSlidingFitResult);
      return ((gapInFirst && gapInThird) || (gapInFirst && endThird));
    }
}
 //------------------------------------------------------------------------------------------------------------------------------------------     
bool TracksCrossingGapsTool::IsEndOfCluster(const float &xSample,  const pandora::Cluster *const pCluster, const TwoDSlidingFitResult &slidingFitResult) const
{
  CartesianVector innerCoordinate(0.f, 0.f, 0.f), outerCoordinate(0.f, 0.f, 0.f);
  LArClusterHelper::GetExtremalCoordinates(pCluster, innerCoordinate, outerCoordinate);
  //std::cout << "IsEndOfCluster" << std::fabs(innerCoordinate.GetX()-xSample) << " and " << std::fabs(outerCoordinate.GetX()-xSample) << std::endl;

  const float pitch(slidingFitResult.GetLayerPitch());
  //std::cout << "pitch = " << pitch << std::endl;
  
  return( ((std::fabs(innerCoordinate.GetX()-xSample)) < pitch ) ||
	  ((std::fabs(outerCoordinate.GetX()-xSample)) < pitch ) );
}

 //------------------------------------------------------------------------------------------------------------------------------------------     
bool TracksCrossingGapsTool::CheckGaps(const float &xSample, const pandora::Cluster *const pCluster, const TwoDSlidingFitResult &slidingFitResult) const
{
  const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
  CartesianVector innerCoordinate(0.f, 0.f, 0.f), outerCoordinate(0.f, 0.f, 0.f);
  LArClusterHelper::GetExtremalCoordinates(pCluster, innerCoordinate, outerCoordinate);

  CartesianVector startPosition(0.f, 0.f, 0.f), startDirection(0.f, 0.f, 0.f);
  
  if((xSample>innerCoordinate.GetX()) && (xSample<outerCoordinate.GetX())) //the xSample point is inside the cluster, should not happen but
    {

    //LORENA: this should only happen when there is a gap within the cluster, and then sliding fit seems to fail to give the position for X 
      startPosition = innerCoordinate;
      startDirection = slidingFitResult.GetGlobalMinLayerDirection();
    } 
  else
    {
      startPosition = ((std::fabs(innerCoordinate.GetX()-xSample) < std::fabs(outerCoordinate.GetX()-xSample)) ?  innerCoordinate :  outerCoordinate);
      startDirection = ((std::fabs(innerCoordinate.GetX()-xSample) < std::fabs(outerCoordinate.GetX()-xSample)) ? (slidingFitResult.GetGlobalMinLayerDirection() * -1.f) :  slidingFitResult.GetGlobalMaxLayerDirection());
    }
    
  const float ext(std::fabs((startPosition.GetX()-xSample)/startDirection.GetX()));//find the factor we need to multiply the direction to reach xSample
  const CartesianVector samplingPoint(startPosition + startDirection * ext);

  //visual debugging
  //LORENA
  //  PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
  PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &samplingPoint, "SamplingPointCheckGaps", BLACK, 1));
     //PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
  //std::cout << "HitType:" << hitType << std::endl;

  return(LArGeometryHelper::IsInGap(this->GetPandora(), samplingPoint, hitType, m_sampleStepSize));
  
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TracksCrossingGapsTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedFraction", m_minMatchedFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedSamplingPoints", m_minMatchedSamplingPoints));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinXOverlapFraction", m_minXOverlapFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
	"MinMatchedSamplingPointRatio", m_minMatchedSamplingPointRatio));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,                                                     
	"MaxGapTolerance", m_maxGapTolerance)); 

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
	"SampleStepSize", m_sampleStepSize));                                                                                                  
    if (m_sampleStepSize < std::numeric_limits<float>::epsilon())                                                                                    
      {                                                                                                                                                         
        std::cout << "TracksCrossingGapsTool: Invalid value for SampleStepSize " << m_sampleStepSize << std::endl;                                   
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);                                                                                 
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
	"MaxAngleRatio", m_maxAngleRatio));
		
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
