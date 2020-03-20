/**
 *  @file   larpandoracontent/LArTwoDReco/HitWidthClusterMergingAlgorithm.cc
 *
 *  @brief  Implementation of the hit width cluster merging algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArTwoDReco/HitWidthClusterMergingAlgorithm.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

using namespace pandora;

namespace lar_content
{


HitWidthClusterMergingAlgorithm::HitWidthClusterMergingAlgorithm() :
  m_clusterListName(),
  m_maxConstituentHitWidth(0.5),
  m_hitWidthScalingFactor(1.0),
  m_useSlidingLinearFit(false),
  m_layerFitHalfWindow(20),
  m_fittingWeight(10),
  m_minClusterWeight(0.5),      
  m_maxXMergeDistance(5),     
  m_maxZMergeDistance(2),     
  m_maxMergeCosOpeningAngle(0.97),
  m_maxDirectionDeviationCosAngle(0.90)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitWidthClusterMergingAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
 
    LArHitWidthHelper::ClusterToParametersMapStore* pClusterToParametersMapStore = LArHitWidthHelper::ClusterToParametersMapStore::Instance();
    LArHitWidthHelper::ClusterToParametersMap* pClusterToParametersMap = pClusterToParametersMapStore->GetMap();

    // Clear map if already full i.e. from other view clustering
    if(!pClusterToParametersMap->empty())
        pClusterToParametersMap->clear();

    for (const Cluster *const pCluster : *pClusterList)
    {
        // The original cluster weight, with no hit scaling or hit padding
        if(LArHitWidthHelper::GetOriginalTotalClusterWeight(pCluster) < m_minClusterWeight)
            continue;

        pClusterToParametersMap->insert(std::pair<const Cluster*, LArHitWidthHelper::ClusterParameters>(pCluster, LArHitWidthHelper::ClusterParameters(pCluster, m_maxConstituentHitWidth, m_useSlidingLinearFit, m_hitWidthScalingFactor)));

        clusterVector.push_back(pCluster);
    }

    //ORDER BY MAX EXTREMAL X COORDINATE
    std::sort(clusterVector.begin(), clusterVector.end(), LArHitWidthHelper::SortByHigherXExtrema);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitWidthClusterMergingAlgorithm::PopulateClusterAssociationMap(const ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const
{

    LArHitWidthHelper::ClusterToParametersMapStore* pClusterToParametersMapStore = LArHitWidthHelper::ClusterToParametersMapStore::Instance();
    LArHitWidthHelper::ClusterToParametersMap* pClusterToParametersMap = pClusterToParametersMapStore->GetMap();

    if(pClusterToParametersMap->empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    // ATTN this method assumes that clusters have been sorted by extremal x position (low x -> high x)
    for (ClusterVector::const_iterator iterCurrentCluster = clusterVector.begin(); iterCurrentCluster != clusterVector.end(); ++iterCurrentCluster)
    {
        const Cluster *const pCurrentCluster = *iterCurrentCluster;
        const LArHitWidthHelper::ClusterToParametersMap::const_iterator currentParametersIter = pClusterToParametersMap->find(pCurrentCluster);

        if(currentParametersIter == pClusterToParametersMap->end())
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);

        for (ClusterVector::const_iterator iterTestCluster = iterCurrentCluster; iterTestCluster != clusterVector.end(); ++iterTestCluster)
        {
            if (iterCurrentCluster == iterTestCluster)
                continue;

	    const Cluster *const pTestCluster = *iterTestCluster;
            const LArHitWidthHelper::ClusterToParametersMap::const_iterator testParametersIter = pClusterToParametersMap->find(pTestCluster);

	    if(testParametersIter == pClusterToParametersMap->end())
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);

            if (!this->AreClustersAssociated(currentParametersIter->second, testParametersIter->second))
	        continue;

            clusterAssociationMap[pCurrentCluster].m_forwardAssociations.insert(pTestCluster);
            clusterAssociationMap[pTestCluster].m_backwardAssociations.insert(pCurrentCluster);
        }
    }

    // remove all 'shortcut' routes
    this->CleanupClusterAssociations(clusterVector, clusterAssociationMap);

}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitWidthClusterMergingAlgorithm::CleanupClusterAssociations(const ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const
{

  // Create temporary map so can delete elements whilst still iterating over them
  ClusterAssociationMap tempMap(clusterAssociationMap);

  for(const Cluster *const pCluster : clusterVector) 
  {
      const ClusterAssociationMap::const_iterator primaryMapIter = clusterAssociationMap.find(pCluster);

      if(primaryMapIter == clusterAssociationMap.end())
          continue;

      const ClusterSet &primaryForwardAssociations(primaryMapIter->second.m_forwardAssociations);

      // loop through the primary associations of each cluster
      // remove clusters that are present in secondary associations of other primary clusters 
      for(const Cluster *const pConsideredCluster : primaryForwardAssociations)
      {
          for(const Cluster *const pPrimaryCluster : primaryForwardAssociations)
          {
              if(pConsideredCluster == pPrimaryCluster)
                  continue;

              const ClusterAssociationMap::const_iterator secondaryMapIter = clusterAssociationMap.find(pPrimaryCluster);

	          // if primary cluster has no associations (this shouldn't ever be the case)
              if(secondaryMapIter == clusterAssociationMap.end()) 
                  continue;

              const ClusterSet &secondaryForwardAssociations(secondaryMapIter->second.m_forwardAssociations);

              // if cluster is present in the forward associations of any cluster at the same level remove
              if(secondaryForwardAssociations.find(pConsideredCluster) != secondaryForwardAssociations.end()) 
              {
                  ClusterSet &tempPrimaryForwardAssociations(tempMap.find(pCluster)->second.m_forwardAssociations);
                  const ClusterSet::const_iterator forwardAssociationToRemove(tempPrimaryForwardAssociations.find(pConsideredCluster));

                  // if association has already been removed
                  if(forwardAssociationToRemove == tempPrimaryForwardAssociations.end())
                      continue;

                  ClusterSet &tempPrimaryBackwardAssociations(tempMap.find(pConsideredCluster)->second.m_backwardAssociations);
                  const ClusterSet::const_iterator backwardAssociationToRemove(tempPrimaryBackwardAssociations.find(pCluster));

                  // if association has already been removed
                  if(backwardAssociationToRemove == tempPrimaryBackwardAssociations.end())
                      continue;

                  tempPrimaryForwardAssociations.erase(forwardAssociationToRemove);
                  tempPrimaryBackwardAssociations.erase(backwardAssociationToRemove);
              }
          }
      }
  }
  clusterAssociationMap = tempMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool HitWidthClusterMergingAlgorithm::AreClustersAssociated(const LArHitWidthHelper::ClusterParameters &currentFitParameters, const LArHitWidthHelper::ClusterParameters &testFitParameters) const
{

    CartesianVector currentClusterDirection(0,0,0);
    CartesianVector testClusterDirection(0,0,0);

    if(m_useSlidingLinearFit)
    {
        const CartesianPointVector currentConstituentHitPositionVector(LArHitWidthHelper::GetConstituentHitPositionVector(currentFitParameters.GetConstituentHitVector()));
        const CartesianPointVector testConstituentHitPositionVector(LArHitWidthHelper::GetConstituentHitPositionVector(testFitParameters.GetConstituentHitVector()));

        const float layerPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));

        TwoDSlidingFitResult currentFitResult(&currentConstituentHitPositionVector, m_layerFitHalfWindow, layerPitch);
        TwoDSlidingFitResult testFitResult(&testConstituentHitPositionVector, m_layerFitHalfWindow, layerPitch);

        StatusCode currentStatus(currentFitResult.GetGlobalFitDirectionAtX(currentFitParameters.GetHigherXExtrema().GetX(), currentClusterDirection));
        StatusCode testStatus(testFitResult.GetGlobalFitDirectionAtX(testFitParameters.GetLowerXExtrema().GetX(), testClusterDirection));
      
        if(currentStatus != STATUS_CODE_SUCCESS || testStatus != STATUS_CODE_SUCCESS)
            return false;
    }
    else
    {
        currentClusterDirection = GetClusterDirection(currentFitParameters, currentFitParameters.GetHigherXExtrema());
        testClusterDirection = GetClusterDirection(testFitParameters, testFitParameters.GetLowerXExtrema());
    }

    if(testFitParameters.GetLowerXExtrema().GetX() > (currentFitParameters.GetHigherXExtrema().GetX() + m_maxXMergeDistance))
      return false;

    if(testFitParameters.GetLowerXExtrema().GetZ() > (currentFitParameters.GetHigherXExtrema().GetZ() + m_maxZMergeDistance) || testFitParameters.GetLowerXExtrema().GetZ() < (currentFitParameters.GetHigherXExtrema().GetZ() - m_maxZMergeDistance))
      return false;
    
    if(fabs(currentClusterDirection.GetCosOpeningAngle(testClusterDirection)) < m_maxMergeCosOpeningAngle)
      return false;

    
    // check that the new direction is consistent with the old clusters
    if(m_useSlidingLinearFit)
    {
        const CartesianPointVector currentConstituentHitPositionVector(LArHitWidthHelper::GetConstituentHitPositionVector(currentFitParameters.GetConstituentHitVector()));
        const CartesianPointVector testConstituentHitPositionVector(LArHitWidthHelper::GetConstituentHitPositionVector(testFitParameters.GetConstituentHitVector()));

        CartesianPointVector newConstituentHitPositionVector(currentConstituentHitPositionVector);
        newConstituentHitPositionVector.insert(newConstituentHitPositionVector.end(), testConstituentHitPositionVector.begin(), testConstituentHitPositionVector.end());
      
        const float layerPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
	
        TwoDSlidingFitResult newFitResult(&newConstituentHitPositionVector, m_layerFitHalfWindow, layerPitch);
	
        CartesianVector newClusterDirectionCurrent(0,0,0);
        CartesianVector newClusterDirectionTest(0,0,0);
				       
        StatusCode currentStatus(newFitResult.GetGlobalFitDirectionAtX(currentFitParameters.GetHigherXExtrema().GetX(), newClusterDirectionCurrent));
        StatusCode testStatus(newFitResult.GetGlobalFitDirectionAtX(testFitParameters.GetLowerXExtrema().GetX(), newClusterDirectionTest));
      
        if(currentStatus != STATUS_CODE_SUCCESS || testStatus != STATUS_CODE_SUCCESS)
            return false;

        if(std::fabs(newClusterDirectionCurrent.GetCosOpeningAngle(currentClusterDirection)) < m_maxDirectionDeviationCosAngle || std::fabs(newClusterDirectionTest.GetCosOpeningAngle(testClusterDirection)) < m_maxDirectionDeviationCosAngle)
            return false;
        
    }
    else
    {
        const LArHitWidthHelper::ConstituentHitVector currentConstituentHitVector(currentFitParameters.GetConstituentHitVector());
        const LArHitWidthHelper::ConstituentHitVector testConstituentHitVector(testFitParameters.GetConstituentHitVector());

        LArHitWidthHelper::ConstituentHitVector newConstituentHitVector(currentConstituentHitVector);
        newConstituentHitVector.insert(newConstituentHitVector.end(), testConstituentHitVector.begin(), testConstituentHitVector.end());
      
        LArHitWidthHelper::ClusterParameters newParameters(nullptr, currentFitParameters.GetNumCaloHits() + testFitParameters.GetNumCaloHits(), currentFitParameters.GetTotalWeight() + testFitParameters.GetTotalWeight(), newConstituentHitVector, CartesianVector(0,0,0), CartesianVector(0,0,0));

        const CartesianVector midpoint((currentFitParameters.GetHigherXExtrema() + testFitParameters.GetLowerXExtrema())*0.5);
        
        CartesianVector newClusterDirection(GetClusterDirection(newParameters, midpoint));
        
        if(std::fabs(newClusterDirection.GetCosOpeningAngle(currentClusterDirection)) < m_maxDirectionDeviationCosAngle || std::fabs(newClusterDirection.GetCosOpeningAngle(testClusterDirection)) < m_maxDirectionDeviationCosAngle)
            return false;
    }
    
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool HitWidthClusterMergingAlgorithm::IsExtremalCluster(const bool isForward, const Cluster *const pCurrentCluster,  const Cluster *const pTestCluster) const
{

    LArHitWidthHelper::ClusterToParametersMapStore* pClusterToParametersMapStore = LArHitWidthHelper::ClusterToParametersMapStore::Instance();
    LArHitWidthHelper::ClusterToParametersMap* pClusterToParametersMap = pClusterToParametersMapStore->GetMap();

    if(pClusterToParametersMap->empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    LArHitWidthHelper::ClusterToParametersMap::const_iterator currentIter(pClusterToParametersMap->find(pCurrentCluster));
    LArHitWidthHelper::ClusterToParametersMap::const_iterator testIter(pClusterToParametersMap->find(pTestCluster));

    if(currentIter == pClusterToParametersMap->end() || testIter == pClusterToParametersMap->end())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    CartesianVector currentHigherXExtrema(currentIter->second.GetHigherXExtrema()), testHigherXExtrema(testIter->second.GetHigherXExtrema());

    float currentMaxX(currentHigherXExtrema.GetX()), testMaxX(testHigherXExtrema.GetX());

    if (isForward)
    {
        if (std::fabs(testMaxX - currentMaxX) > std::numeric_limits<float>::epsilon())
            return (testMaxX > currentMaxX);
    }
    else
    {
        if (std::fabs(testMaxX - currentMaxX) > std::numeric_limits<float>::epsilon())
            return (testMaxX < currentMaxX);
    }
  
    return LArClusterHelper::SortByNHits(pTestCluster, pCurrentCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector HitWidthClusterMergingAlgorithm::GetClusterDirection(const LArHitWidthHelper::ClusterParameters &clusterFitParameters, const CartesianVector &fitReferencePoint) const 
{

    // If cluster composed of one hit, return a transverse line as fit
    if(clusterFitParameters.GetNumCaloHits() == 1)
        return CartesianVector(1, 0, 0);

    // minimise the longitudinal distance in fit (works best for transverse fits)
    CartesianVector LSTransverseClusterFitDirection(0,0,0);
    CartesianVector LSTransverseIntercept(0,0,0);
    float LSTransverseChiSquared(0);

    // minimise the transverse coordinate in fit (works best for longitudinal fits)
    CartesianVector LSLongitudinalClusterFitDirection(0,0,0);
    CartesianVector LSLongitudinalIntercept(0,0,0);
    float LSLongitudinalChiSquared(0);

    GetWeightedGradient(clusterFitParameters, true, LSTransverseClusterFitDirection, LSTransverseIntercept, LSTransverseChiSquared, fitReferencePoint);
    GetWeightedGradient(clusterFitParameters, false, LSLongitudinalClusterFitDirection, LSLongitudinalIntercept, LSLongitudinalChiSquared, fitReferencePoint);
    
    // return fit with the lowest chi-squared
    if(LSTransverseChiSquared < LSLongitudinalChiSquared)
        return LSTransverseClusterFitDirection;

    return LSLongitudinalClusterFitDirection;
}
  
//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector HitWidthClusterMergingAlgorithm::GetClusterZIntercept(const LArHitWidthHelper::ClusterParameters &clusterFitParameters, const CartesianVector &fitReferencePoint) const
{

    // If cluster composed of one hit, return a transverse line as fit
    if(clusterFitParameters.GetNumCaloHits() == 1)
        return CartesianVector(0, 0, clusterFitParameters.GetConstituentHitVector().begin()->GetPositionVector().GetZ());

    // minimise the longitudinal distance in fit (works best for transverse fits)
    CartesianVector LSTransverseClusterFitDirection(0,0,0);
    CartesianVector LSTransverseIntercept(0,0,0);
    float LSTransverseChiSquared(0);

    // minimise the transverse coordinate in fit (works best for longitudinal fits)
    CartesianVector LSLongitudinalClusterFitDirection(0,0,0);
    CartesianVector LSLongitudinalIntercept(0,0,0);
    float LSLongitudinalChiSquared(0);

    GetWeightedGradient(clusterFitParameters, true, LSTransverseClusterFitDirection, LSTransverseIntercept, LSTransverseChiSquared, fitReferencePoint);
    GetWeightedGradient(clusterFitParameters, false, LSLongitudinalClusterFitDirection, LSLongitudinalIntercept, LSLongitudinalChiSquared, fitReferencePoint);

    // return fit with the lowest chi-squared
    if(LSTransverseChiSquared < LSLongitudinalChiSquared)
        return LSTransverseIntercept;

    return LSLongitudinalIntercept;
}

//------------------------------------------------------------------------------------------------------------------------------------------


void HitWidthClusterMergingAlgorithm::GetWeightedGradient(const LArHitWidthHelper::ClusterParameters &clusterFitParameters, bool isTransverse, CartesianVector &direction, CartesianVector &intercept, float &chiSquared, const CartesianVector &fitReferencePoint) const
{

    // cannot make a longitudinal fit to a single hit cluster
    if(!isTransverse && clusterFitParameters.GetNumCaloHits() == 1) 
    {
        std::cout << "WARNING - CANNOT MAKE LONGITUDINAL FIT TO SINGLE HIT CLUSTER, RETURNED" << std::endl;
        throw;
    }

    float weightSum(0);
    float weightedXSum(0);
    float weightedZSum(0);

    bool isXConstant(true);
    bool isZConstant(true);

    LArHitWidthHelper::ConstituentHitVector constituentHitVector(clusterFitParameters.GetConstituentHitVector());

    // sort hits with respect to their distance to the fitReferencePoint (closest -> furthest)
    std::sort(constituentHitVector.begin(), constituentHitVector.end(), LArHitWidthHelper::ConstituentHit::SortByDistanceToPoint(fitReferencePoint));

    for(const LArHitWidthHelper::ConstituentHit &constituentHit : constituentHitVector) 
    {
        const CartesianVector hitPosition = constituentHit.GetPositionVector();
        const float hitWeight = constituentHit.GetHitWidth();

        if(std::fabs(clusterFitParameters.GetConstituentHitVector().begin()->GetPositionVector().GetX() - hitPosition.GetX()) > std::numeric_limits<float>::epsilon())
            isXConstant = false;

        if(std::fabs(clusterFitParameters.GetConstituentHitVector().begin()->GetPositionVector().GetZ() - hitPosition.GetZ()) > std::numeric_limits<float>::epsilon())
            isZConstant = false;

        weightedXSum += hitPosition.GetX() * hitWeight;
        weightedZSum += hitPosition.GetZ() * hitWeight;
        weightSum += hitWeight;
        
        if(weightSum > m_fittingWeight)
            break;
    }

    // TO FIT A STRAIGHT LINE TO CLUSTERS WITH CONSTANT X OR Z
    if(isXConstant) 
    {
        direction = CartesianVector(0, 0, 1);
        intercept = CartesianVector(0, 0, 0);
        chiSquared = 0;
        return;
    }

    if(isZConstant) 
    {
        direction = CartesianVector(1, 0, 0);
        intercept = CartesianVector(0, 0, clusterFitParameters.GetConstituentHitVector().begin()->GetPositionVector().GetZ());
        chiSquared = 0;
        return;
    }

    float weightedXMean(weightedXSum/weightSum);
    float weightedZMean(weightedZSum/weightSum); 

    float numerator(0);
    float denominator(0);
    float chi(0);

    float weightCount(0);
    for(const LArHitWidthHelper::ConstituentHit &constituentHit : constituentHitVector) 
    {
        const CartesianVector hitPosition = constituentHit.GetPositionVector();
        const float hitWeight = constituentHit.GetHitWidth();
        weightCount += hitWeight;
        
        numerator += hitWeight * (hitPosition.GetX() - weightedXMean) * (hitPosition.GetZ() - weightedZMean);
        isTransverse ? denominator += hitWeight * pow(hitPosition.GetX() - weightedXMean, 2) : denominator += hitWeight * pow(hitPosition.GetZ() - weightedZMean, 2);

        if(weightCount > m_fittingWeight)
            break;
    }

    float gradient = numerator/denominator;
    isTransverse ? intercept.SetValues(0, 0, weightedZMean - gradient * weightedXMean) : intercept.SetValues(weightedXMean - gradient * weightedZMean, 0, 0);

    weightCount = 0;
    for(const LArHitWidthHelper::ConstituentHit &constituentHit : constituentHitVector) 
    {
        const CartesianVector hitPosition = constituentHit.GetPositionVector();
        const float hitWeight = constituentHit.GetHitWidth();
        weightCount += hitWeight;
        
        isTransverse ? chi += hitWeight*pow(hitPosition.GetZ() - intercept.GetZ() - gradient*hitPosition.GetX(), 2) : chi += hitWeight * pow(hitPosition.GetX() - intercept.GetX() - gradient * hitPosition.GetZ(), 2);

        if(weightCount > m_fittingWeight)
            break;
    }

    // change coordinates to y=mx+c fit and normalise
    isTransverse? direction = CartesianVector(1.0, 0, gradient).GetUnitVector() : direction = CartesianVector(gradient, 0, 1.0).GetUnitVector();

    if(!isTransverse)
        intercept.SetValues(0, 0, -intercept.GetX()/gradient);

    chiSquared = chi;

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------


StatusCode HitWidthClusterMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, 
        "ClusterListName", m_clusterListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "UseSlidingLinearFit", m_useSlidingLinearFit));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FittingWeight", m_fittingWeight));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "LayerFitHalfWindow", m_layerFitHalfWindow));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterWeight", m_minClusterWeight));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxXMergeDistance", m_maxXMergeDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxZMergeDistance", m_maxZMergeDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxMergeCosOpeningAngle", m_maxMergeCosOpeningAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxDirectionDeviationCosAngle", m_maxDirectionDeviationCosAngle));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxConstituentHitWidth", m_maxConstituentHitWidth));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "HitWidthScalingFactor", m_hitWidthScalingFactor));


    return ClusterAssociationAlgorithm::ReadSettings(xmlHandle);
}
  
} // namespace lar_content



