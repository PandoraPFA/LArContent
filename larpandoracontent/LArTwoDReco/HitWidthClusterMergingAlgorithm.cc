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



using namespace pandora;

namespace lar_content
{


HitWidthClusterMergingAlgorithm::HitWidthClusterMergingAlgorithm() :
  m_clusterListName(),
  m_minClusterWeight(0.5),      
  m_maxXMergeDistance(5),     
  m_maxZMergeDistance(2),     
  m_maxMergeCosOpeningAngle(0.97), 
  m_maxConstituentHitWidth(0.5)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------


void HitWidthClusterMergingAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    for (const Cluster *const pCluster : *pClusterList)
    {

        if(LArHitWidthHelper::GetTotalClusterWeight(pCluster) < m_minClusterWeight)
            continue;

        clusterVector.push_back(pCluster);
    }

    //ORDER BY MAX EXTREMAL X COORDINATE
    std::sort(clusterVector.begin(), clusterVector.end(), SortByHigherXExtrema(m_maxConstituentHitWidth));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitWidthClusterMergingAlgorithm::PopulateClusterAssociationMap(const ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const
{

    LArHitWidthHelper::ClusterToParametersMap clusterToFitParametersMap;

    for (const Cluster *const pCluster : clusterVector)
    {
        clusterToFitParametersMap.insert(std::pair(pCluster, LArHitWidthHelper::ClusterParameters(pCluster, m_maxConstituentHitWidth)));
    }


    // ATTN This method assumes that clusters have been sorted by x position (low x -> high x)
    for (ClusterVector::const_iterator iterCurrentCluster = clusterVector.begin(); iterCurrentCluster != clusterVector.end(); ++iterCurrentCluster)
    {

        const Cluster *const pCurrentCluster = *iterCurrentCluster;
	const LArHitWidthHelper::ClusterParameters currentFitParameters = clusterToFitParametersMap.at(pCurrentCluster);

        for (ClusterVector::const_iterator iterTestCluster = iterCurrentCluster; iterTestCluster != clusterVector.end(); ++iterTestCluster)
        {
            if (iterCurrentCluster == iterTestCluster)
                continue;

	    const Cluster *const pTestCluster = *iterTestCluster;
            const LArHitWidthHelper::ClusterParameters testFitParameters = clusterToFitParametersMap.at(pTestCluster);

            if (!this->AreClustersAssociated(currentFitParameters, testFitParameters))
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
    CartesianVector currentClusterDirection(GetClusterDirection(currentFitParameters));
    CartesianVector testClusterDirection(GetClusterDirection(testFitParameters));

    if(testFitParameters.m_lowerXExtrema.GetX() > (currentFitParameters.m_higherXExtrema.GetX() + m_maxXMergeDistance)) {  
      return false;
    }

    if(testFitParameters.m_lowerXExtrema.GetZ() > (currentFitParameters.m_higherXExtrema.GetZ() + m_maxZMergeDistance) || testFitParameters.m_lowerXExtrema.GetZ() < (currentFitParameters.m_higherXExtrema.GetZ() - m_maxZMergeDistance)) {
      return false;
    }
    
    if(fabs(currentClusterDirection.GetCosOpeningAngle(testClusterDirection)) < m_maxMergeCosOpeningAngle) {
      return false;
    }
    
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool HitWidthClusterMergingAlgorithm::IsExtremalCluster(const bool isForward, const Cluster *const pCurrentCluster,  const Cluster *const pTestCluster) const
{

    LArHitWidthHelper::ConstituentHitVector currentConstituentHitVector(LArHitWidthHelper::GetConstituentHits(pCurrentCluster, m_maxConstituentHitWidth));
    LArHitWidthHelper::ConstituentHitVector testConstituentHitVector(LArHitWidthHelper::GetConstituentHits(pTestCluster, m_maxConstituentHitWidth));

    CartesianVector currentHigherXExtrema(LArHitWidthHelper::GetExtremalCoordinatesHigherX(currentConstituentHitVector));
    CartesianVector testHigherXExtrema(LArHitWidthHelper::GetExtremalCoordinatesHigherX(testConstituentHitVector));

    float testMaxX(testHigherXExtrema.GetX()), currentMaxX(currentHigherXExtrema.GetX());

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

  CartesianVector HitWidthClusterMergingAlgorithm::GetClusterDirection(const LArHitWidthHelper::ClusterParameters &clusterFitParameters) const 
{

    // If cluster composed of one hit, return a transverse line as fit
    if(clusterFitParameters.m_numCaloHits == 1) {
      return CartesianVector(1, 0, 0);
    }
   

    // minimise the longitudinal distance in fit (works best for transverse fits)
    CartesianVector LSTransverseClusterFitDirection(0,0,0);
    CartesianVector LSTransverseIntercept(0,0,0);
    float LSTransverseChiSquared(0);

    // minimise the transverse coordinate in fit (works best for longitudinal fits)
    CartesianVector LSLongitudinalClusterFitDirection(0,0,0);
    CartesianVector LSLongitudinalIntercept(0,0,0);
    float LSLongitudinalChiSquared(0);
    
    GetWeightedGradient(clusterFitParameters, true, LSTransverseClusterFitDirection, LSTransverseIntercept, LSTransverseChiSquared);    
    GetWeightedGradient(clusterFitParameters, false, LSLongitudinalClusterFitDirection, LSLongitudinalIntercept, LSLongitudinalChiSquared);

    // return fit with the lowest chi-squared
    if(LSTransverseChiSquared < LSLongitudinalChiSquared)
      return LSTransverseClusterFitDirection;

    return LSLongitudinalClusterFitDirection;

}
  
//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector HitWidthClusterMergingAlgorithm::GetClusterZIntercept(const LArHitWidthHelper::ClusterParameters &clusterFitParameters) const
{

    // If cluster composed of one hit, return a transverse line as fit
    if(clusterFitParameters.m_numCaloHits == 1) {
      return CartesianVector(0, 0, clusterFitParameters.m_constituentHitVector.begin()->m_positionVector.GetZ());
    }      

    // minimise the longitudinal distance in fit (works best for transverse fits)
    CartesianVector LSTransverseClusterFitDirection(0,0,0);
    CartesianVector LSTransverseIntercept(0,0,0);
    float LSTransverseChiSquared(0);

    // minimise the transverse coordinate in fit (works best for longitudinal fits)
    CartesianVector LSLongitudinalClusterFitDirection(0,0,0);
    CartesianVector LSLongitudinalIntercept(0,0,0);
    float LSLongitudinalChiSquared(0);

    GetWeightedGradient(clusterFitParameters, true, LSTransverseClusterFitDirection, LSTransverseIntercept, LSTransverseChiSquared);
    GetWeightedGradient(clusterFitParameters, false, LSLongitudinalClusterFitDirection, LSLongitudinalIntercept, LSLongitudinalChiSquared);

    // return fit with the lowest chi-squared
    if(LSTransverseChiSquared < LSLongitudinalChiSquared)
      return LSTransverseIntercept;

    return LSLongitudinalIntercept;

}

//------------------------------------------------------------------------------------------------------------------------------------------


  void HitWidthClusterMergingAlgorithm::GetWeightedGradient(const LArHitWidthHelper::ClusterParameters &clusterFitParameters, bool isTransverse, CartesianVector &direction, CartesianVector &intercept, float &chiSquared) const
{

    // cannot make a longitudinal fit to a single hit cluster
    if(!isTransverse && clusterFitParameters.m_numCaloHits == 1) 
    {
      std::cout << "WARNING - CANNOT MAKE LONGITUDINAL FIT TO SINGLE HIT CLUSTER, RETURNED" << std::endl;
      throw;
    }

    float weightSum(0);
    float weightedXSum(0);
    float weightedZSum(0);

    bool isXConstant(true);
    bool isZConstant(true);

    LArHitWidthHelper::ConstituentHitVector constituentHitVector(clusterFitParameters.m_constituentHitVector);

    // Include in fit hits nearest the relevant x extrema
    // Stop including hits when cumulative weight exceeds m_fittingSampleWeight 
    // or have already included all cluster constituent hits     
    for(const LArHitWidthHelper::ConstituentHit &constituentHit : constituentHitVector) 
    {
        const CartesianVector hitPosition = constituentHit.m_positionVector;
	const float hitWeight = constituentHit.m_hitWidth;

	if(std::fabs(clusterFitParameters.m_constituentHitVector.begin()->m_positionVector.GetX() - hitPosition.GetX()) > std::numeric_limits<float>::epsilon())
	  isXConstant = false;

	if(std::fabs(clusterFitParameters.m_constituentHitVector.begin()->m_positionVector.GetZ() - hitPosition.GetZ()) > std::numeric_limits<float>::epsilon())
	  isZConstant = false;

        weightedXSum += hitPosition.GetX() * hitWeight;
	weightedZSum += hitPosition.GetZ() * hitWeight;
    }

    // TO FIT A STRAIGHT LINE TO CLUSTERS WITH CONSTANT X OR Z
    if(isXConstant) 
    {
      direction = CartesianVector(0, 0, 1);
      intercept = CartesianVector(0, 0, -std::numeric_limits<float>::max());
      chiSquared = 0;
      return;
    }

    if(isZConstant) 
    {
      direction = CartesianVector(1, 0, 0);
      intercept = CartesianVector(0, 0, clusterFitParameters.m_constituentHitVector.begin()->m_positionVector.GetZ());
      chiSquared = 0;
      return;
    }

    float weightedXMean(weightedXSum/weightSum);
    float weightedZMean(weightedZSum/weightSum); 

    float numerator(0);
    float denominator(0);
    float chi(0);


    for(const LArHitWidthHelper::ConstituentHit &constituentHit : constituentHitVector) 
    {
        const CartesianVector hitPosition = constituentHit.m_positionVector;
	const float hitWeight = constituentHit.m_hitWidth;

	numerator += hitWeight * (hitPosition.GetX() - weightedXMean) * (hitPosition.GetZ() - weightedZMean);
	isTransverse ? denominator += hitWeight * pow(hitPosition.GetX() - weightedXMean, 2) : denominator += hitWeight * pow(hitPosition.GetZ() - weightedZMean, 2);
    }

    float gradient = numerator/denominator;
    isTransverse ? intercept.SetValues(0, 0, weightedZMean - gradient * weightedXMean) : intercept.SetValues(weightedXMean - gradient * weightedZMean, 0, 0);


    for(const LArHitWidthHelper::ConstituentHit &constituentHit : constituentHitVector) 
    {
        const CartesianVector hitPosition = constituentHit.m_positionVector;
	const float hitWeight = constituentHit.m_hitWidth;

	isTransverse ? chi += hitWeight*pow(hitPosition.GetZ() - intercept.GetZ() - gradient*hitPosition.GetX(), 2) : chi += hitWeight * pow(hitPosition.GetX() - intercept.GetX() - gradient * hitPosition.GetZ(), 2);
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

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListName", m_clusterListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterWeight", m_minClusterWeight));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxXMergeDistance", m_maxXMergeDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxZMergeDistance", m_maxZMergeDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxMergeCosOpeningAngle", m_maxMergeCosOpeningAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxConstituentHitWidth", m_maxConstituentHitWidth));

    return ClusterAssociationAlgorithm::ReadSettings(xmlHandle);
}
  
} // namespace lar_content



