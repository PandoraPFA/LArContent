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
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"

using namespace pandora;

namespace lar_content
{

HitWidthClusterMergingAlgorithm::HitWidthClusterMergingAlgorithm() :
  m_maxConstituentHitWidth(0.5f),
  m_hitWidthScalingFactor(1.f),
  m_fittingWeight(10.f),
  m_minClusterWeight(0.5f),      
  m_maxXMergeDistance(5.f),     
  m_maxZMergeDistance(2.f),     
  m_minMergeCosOpeningAngle(0.97f),
  m_minDirectionDeviationCosAngle(0.9f),
  m_minClusterSparseness(0.3f),
  m_useOldDirectionMethod(true),
  m_useClosestMergePoint(false),
  m_doubleFittingWeight(false)
{
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

void HitWidthClusterMergingAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    // clear map if already full i.e. from other view clustering
    if (!m_clusterToParametersMap.empty())
        m_clusterToParametersMap.clear();
    
    for (const Cluster *const pCluster : *pClusterList)
    {
        // the original cluster weight, with no hit scaling or hit padding
        if (LArHitWidthHelper::GetOriginalTotalClusterWeight(pCluster) < m_minClusterWeight)
            continue;
        
        m_clusterToParametersMap.insert(std::pair<const Cluster*, LArHitWidthHelper::ClusterParameters>(pCluster,
            LArHitWidthHelper::ClusterParameters(pCluster, m_maxConstituentHitWidth, false, m_hitWidthScalingFactor)));
        

        LArHitWidthHelper::ClusterParameters clusterParameters(pCluster, m_maxConstituentHitWidth, false, m_hitWidthScalingFactor);

        float clusterSparseness(1.f - (static_cast<float>(pCluster->GetNCaloHits()) / static_cast<float>(clusterParameters.GetConstituentHitVector().size())));

        if (clusterSparseness < m_minClusterSparseness && pCluster->GetNCaloHits() != 1)
            continue;
        
        /*
        ClusterList theCluster;
        theCluster.push_back(pCluster);
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &theCluster, "CLUSTER", VIOLET);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        */                                                
        
        clusterVector.push_back(pCluster);
    }

    std::sort(clusterVector.begin(), clusterVector.end(), LArHitWidthHelper::SortByHigherXExtrema(m_clusterToParametersMap));

    /*
    ClusterList consideredClusters(clusterVector.begin(), clusterVector.end());
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &consideredClusters, "CLUSTER", VIOLET);
    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    */
    
    
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitWidthClusterMergingAlgorithm::PopulateClusterAssociationMap(const ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const
{
    // ATTN this method assumes that clusters have been sorted by extremal x position (low higherXExtrema -> high higherXExtrema)
    for (ClusterVector::const_iterator iterCurrentCluster = clusterVector.begin(); iterCurrentCluster != clusterVector.end(); ++iterCurrentCluster)
    {
        const Cluster *const pCurrentCluster = *iterCurrentCluster;
        const LArHitWidthHelper::ClusterParameters &currentClusterParameters(LArHitWidthHelper::GetClusterParameters(pCurrentCluster, m_clusterToParametersMap));

        for (ClusterVector::const_iterator iterTestCluster = iterCurrentCluster; iterTestCluster != clusterVector.end(); ++iterTestCluster)
        {
            if (iterCurrentCluster == iterTestCluster)
                continue;

            const Cluster *const pTestCluster = *iterTestCluster;
            const LArHitWidthHelper::ClusterParameters &testClusterParameters(LArHitWidthHelper::GetClusterParameters(pTestCluster, m_clusterToParametersMap));

            if (!this->AreClustersAssociated(currentClusterParameters, testClusterParameters))
                continue;

            clusterAssociationMap[pCurrentCluster].m_forwardAssociations.insert(pTestCluster);
            clusterAssociationMap[pTestCluster].m_backwardAssociations.insert(pCurrentCluster);
        }
    }
    
    this->RemoveShortcutAssociations(clusterVector, clusterAssociationMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
bool HitWidthClusterMergingAlgorithm::IsExtremalCluster(const bool isForward, const Cluster *const pCurrentCluster,  const Cluster *const pTestCluster) const
{
    //ATTN - cannot use map since higherXExtrema may have changed during merging
    const LArHitWidthHelper::ConstituentHitVector currentConstituentHitVector(LArHitWidthHelper::GetConstituentHits(pCurrentCluster, m_maxConstituentHitWidth, m_hitWidthScalingFactor, false));
    const LArHitWidthHelper::ConstituentHitVector testConstituentHitVector(LArHitWidthHelper::GetConstituentHits(pTestCluster, m_maxConstituentHitWidth, m_hitWidthScalingFactor, false));
    const CartesianVector currentHigherXExtrema(LArHitWidthHelper::GetExtremalCoordinatesHigherX(currentConstituentHitVector));
    const CartesianVector testHigherXExtrema(LArHitWidthHelper::GetExtremalCoordinatesHigherX(testConstituentHitVector));
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

void HitWidthClusterMergingAlgorithm::FindClosestPointToPosition(const CartesianVector &position, const LArHitWidthHelper::ConstituentHitVector &constituentHitVector,
    CartesianVector &closestPoint) const
{
    float minDistanceSquared(std::numeric_limits<float>::max());
    for (const LArHitWidthHelper::ConstituentHit &constituentHit : constituentHitVector)
    {
        const CartesianVector &hitPosition(constituentHit.GetPositionVector());

        const float separationDistanceSquared(hitPosition.GetDistanceSquared(position));
        if (separationDistanceSquared < minDistanceSquared)
        {
            minDistanceSquared = separationDistanceSquared;
            closestPoint = hitPosition;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool HitWidthClusterMergingAlgorithm::AreClustersAssociated(const LArHitWidthHelper::ClusterParameters &currentFitParameters,
    const LArHitWidthHelper::ClusterParameters &testFitParameters) const
{    
    // check merging points not too far away in x
    if (testFitParameters.GetLowerXExtrema().GetX() > (currentFitParameters.GetHigherXExtrema().GetX() + m_maxXMergeDistance))
        return false;

    // check merging points not too far away in z
    if (testFitParameters.GetLowerXExtrema().GetZ() > (currentFitParameters.GetHigherXExtrema().GetZ() + m_maxZMergeDistance) ||
        testFitParameters.GetLowerXExtrema().GetZ() < (currentFitParameters.GetHigherXExtrema().GetZ() - m_maxZMergeDistance))
    {
        return false;
    }

    // find appropriate test merging point and check closeness in z
    CartesianVector testMergePoint(0.f, 0.f, 0.f);
    if (testFitParameters.GetLowerXExtrema().GetX() < currentFitParameters.GetHigherXExtrema().GetX())
    {
        if (m_useClosestMergePoint)
        {
            this->FindClosestPointToPosition(currentFitParameters.GetHigherXExtrema(), testFitParameters.GetConstituentHitVector(), testMergePoint);
        }
        else
        {
            testMergePoint = testFitParameters.GetLowerXExtrema();
        }
    }
    else
    {
        testMergePoint = testFitParameters.GetLowerXExtrema();
    }

    // check merging points not too far away in z
    if (testMergePoint.GetZ() > (currentFitParameters.GetHigherXExtrema().GetZ() + m_maxZMergeDistance) ||
        testMergePoint.GetZ() < (currentFitParameters.GetHigherXExtrema().GetZ() - m_maxZMergeDistance))
    {
        return false;
    }


    
    CartesianVector currentClusterDirection(0.f, 0.f, 0.f), testClusterDirection(0.f, 0.f, 0.f);
    if (m_useOldDirectionMethod)
    {
        currentClusterDirection = this->GetClusterDirection(currentFitParameters.GetConstituentHitVector(), currentFitParameters.GetHigherXExtrema());
        testClusterDirection = this->GetClusterDirection(testFitParameters.GetConstituentHitVector(), testMergePoint);
    }
    else
    {
        this->GetWeightedGradient(currentFitParameters.GetConstituentHitVector(), currentClusterDirection, currentFitParameters.GetHigherXExtrema(), m_fittingWeight);
        this->GetWeightedGradient(testFitParameters.GetConstituentHitVector(), testClusterDirection, testMergePoint, m_fittingWeight);
    }
    
    // check clusters have a similar direction
    if (currentClusterDirection.GetCosOpeningAngle(testClusterDirection) < m_minMergeCosOpeningAngle)
        return false;

    // check that the new direction is consistent with the old clusters
    LArHitWidthHelper::ConstituentHitVector newConstituentHitVector(currentFitParameters.GetConstituentHitVector());
    newConstituentHitVector.insert(newConstituentHitVector.end(), testFitParameters.GetConstituentHitVector().begin(), testFitParameters.GetConstituentHitVector().end());
    
    const CartesianVector midpoint((currentFitParameters.GetHigherXExtrema() + testMergePoint) * 0.5);
    CartesianVector newClusterDirection(0.f, 0.f, 0.f);
    if (m_useOldDirectionMethod)
    {
        newClusterDirection = this->GetClusterDirection(newConstituentHitVector, midpoint);
    }
    else
    {
        if (m_doubleFittingWeight)
        {
            this->GetWeightedGradient(newConstituentHitVector, newClusterDirection, midpoint, m_fittingWeight * 2.0);
        }
        else
        {
            this->GetWeightedGradient(newConstituentHitVector, newClusterDirection, midpoint, m_fittingWeight);
        }
    }

    if (newClusterDirection.GetCosOpeningAngle(currentClusterDirection) < m_minDirectionDeviationCosAngle ||
        newClusterDirection.GetCosOpeningAngle(testClusterDirection) < m_minDirectionDeviationCosAngle)
    {
        return false;
    }

    ///////////////////////
    /*
    ClusterList theTestCluster;
    theTestCluster.push_back(testFitParameters.GetClusterAddress());
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &theTestCluster, "TEST CLUSTER", RED);
    ClusterList theCurrentCluster;
    theCurrentCluster.push_back(currentFitParameters.GetClusterAddress());
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &theCurrentCluster, "CURRENT CLUSTER", BLUE);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testMergePoint, "CLOSEST MERGE POINT", RED, 2);
    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    */
    //////////////////////
    
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector HitWidthClusterMergingAlgorithm::GetClusterDirection(const LArHitWidthHelper::ConstituentHitVector &constituentHitVector,
    const CartesianVector &fitReferencePoint) const 
{
    // minimise the transverse distance in fit
    CartesianVector lsTransverseClusterFitDirection(0.f, 0.f, 0.f), lsTransverseIntercept(0.f, 0.f, 0.f);
    float lsTransverseChiSquared(0.f);

    // minimise the longitudinal distance in fit
    CartesianVector lsLongitudinalClusterFitDirection(0.f, 0.f, 0.f), lsLongitudinalIntercept(0.f, 0.f, 0.f);
    float lsLongitudinalChiSquared(0.f);

    this->GetWeightedGradient(constituentHitVector, false, lsTransverseClusterFitDirection, lsTransverseIntercept, lsTransverseChiSquared, fitReferencePoint);
    this->GetWeightedGradient(constituentHitVector, true, lsLongitudinalClusterFitDirection, lsLongitudinalIntercept, lsLongitudinalChiSquared, fitReferencePoint);
    
    // return fit with the lowest chi-squared
    if (lsTransverseChiSquared < lsLongitudinalChiSquared)
        return lsTransverseClusterFitDirection;

    return lsLongitudinalClusterFitDirection;
}
  
//------------------------------------------------------------------------------------------------------------------------------------------

void HitWidthClusterMergingAlgorithm::GetWeightedGradient(const LArHitWidthHelper::ConstituentHitVector &unsortedConstituentHitVector, const bool isLongitudinal,
        CartesianVector &direction, CartesianVector &zIntercept, float &chiSquared, const CartesianVector &fitReferencePoint) const
{
    float weightSum(0.f), weightedXSum(0.f), weightedZSum(0.f);
    bool isXConstant(true), isZConstant(true);
    LArHitWidthHelper::ConstituentHitVector constituentHitVector(unsortedConstituentHitVector);
    
    // sort hits with respect to their distance to the fitReferencePoint (closest -> furthest)
    std::sort(constituentHitVector.begin(), constituentHitVector.end(), LArHitWidthHelper::ConstituentHit::SortByDistanceToPoint(fitReferencePoint));
    
    // calculate weightedXMean and weightedZMean for fit
    for (const LArHitWidthHelper::ConstituentHit &constituentHit : constituentHitVector) 
    {
        const CartesianVector &hitPosition(constituentHit.GetPositionVector());
        const float hitWeight(constituentHit.GetHitWidth());

        if (std::fabs(constituentHitVector.begin()->GetPositionVector().GetX() - hitPosition.GetX()) > std::numeric_limits<float>::epsilon())
            isXConstant = false;

        if (std::fabs(constituentHitVector.begin()->GetPositionVector().GetZ() - hitPosition.GetZ()) > std::numeric_limits<float>::epsilon())
            isZConstant = false;

        weightedXSum += hitPosition.GetX() * hitWeight;
        weightedZSum += hitPosition.GetZ() * hitWeight;
        weightSum += hitWeight;

        // do not exceed the specified hit weight in the fit
        if(weightSum > m_fittingWeight)
            break;
    }

    if (weightSum < std::numeric_limits<float>::epsilon())
    {
        std::cout << "HitWidthClusterMergingAlgorithm::GetWeightedGradient - hit weight in fit is negative or equivalent to zero" << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);
    }

    // return vertical fit for a cluster with constant x
    if (isXConstant) 
    {
        direction = CartesianVector(0.f, 0.f, 1.f);
        zIntercept = CartesianVector(0.f, 0.f, 0.f);
        chiSquared = 0.f;
        return;
    }

    // return horizontal fit for a cluster with constant z
    if (isZConstant) 
    {
        direction = CartesianVector(1.f, 0.f, 0.f);
        zIntercept = CartesianVector(0.f, 0.f, constituentHitVector.begin()->GetPositionVector().GetZ());
        chiSquared = 0.f;
        return;
    }

    float weightedXMean(weightedXSum/weightSum), weightedZMean(weightedZSum/weightSum), numerator(0.f), denominator(0.f), weightCount(0.f);
    for (const LArHitWidthHelper::ConstituentHit &constituentHit : constituentHitVector) 
    {
        const CartesianVector &hitPosition(constituentHit.GetPositionVector());
        const float hitWeight(constituentHit.GetHitWidth());
        weightCount += hitWeight;
        
        numerator += hitWeight * (hitPosition.GetX() - weightedXMean) * (hitPosition.GetZ() - weightedZMean);
        denominator += isLongitudinal ? hitWeight * pow(hitPosition.GetX() - weightedXMean, 2) : hitWeight * pow(hitPosition.GetZ() - weightedZMean, 2);

        // do not exceed the specified hit weight in the fit
        if (weightCount > m_fittingWeight)
            break;
    }

    float gradient(numerator/denominator), chi(0.f);
    float intercept(isLongitudinal ? weightedZMean - gradient * weightedXMean : weightedXMean - gradient * weightedZMean);

    weightCount = 0.f;
    for (const LArHitWidthHelper::ConstituentHit &constituentHit : constituentHitVector) 
    {
        const CartesianVector hitPosition(constituentHit.GetPositionVector());
        const float hitWeight(constituentHit.GetHitWidth());
        weightCount += hitWeight;
        
        chi += isLongitudinal ? hitWeight * pow(hitPosition.GetZ() - intercept - gradient * hitPosition.GetX(), 2) : hitWeight * pow(hitPosition.GetX() - intercept - gradient * hitPosition.GetZ(), 2);

        // do not exceed the specified hit weight in the fit
        if (weightCount > m_fittingWeight)
            break;
    }

    // change coordinates to z=mx+c fit and normalise
    direction = isLongitudinal ? CartesianVector(1.f, 0.f, gradient).GetUnitVector() : CartesianVector(gradient, 0.f, 1.f).GetUnitVector();

    if (direction.GetX() < std::numeric_limits<float>::epsilon())
        direction = direction * (-1.f);
    
    zIntercept = isLongitudinal ? CartesianVector(0.f, 0.f, intercept) : CartesianVector(0.f, 0.f, -intercept/gradient);
    chiSquared = chi;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitWidthClusterMergingAlgorithm::GetConstituentHitSubsetVector(const LArHitWidthHelper::ConstituentHitVector &constituentHitVector, const CartesianVector &fitReferencePoint,
    const float fittingWeight, LArHitWidthHelper::ConstituentHitVector &constituentHitSubsetVector) const
{
    LArHitWidthHelper::ConstituentHitVector sortedConstituentHitVector(constituentHitVector);

    // sort hits with respect to their distance to the fitReferencePoint (closest -> furthest)
    std::sort(sortedConstituentHitVector.begin(), sortedConstituentHitVector.end(), LArHitWidthHelper::ConstituentHit::SortByDistanceToPoint(fitReferencePoint));

    float weightCount(0.f);
    for (const LArHitWidthHelper::ConstituentHit &constituentHit : sortedConstituentHitVector)
    {
        constituentHitSubsetVector.push_back(constituentHit);

        weightCount += constituentHit.GetHitWidth();
        
        if (weightCount > fittingWeight)
            break;
    }    
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitWidthClusterMergingAlgorithm::GetFittingAxes(const LArHitWidthHelper::ConstituentHitVector &constituentHitSubsetVector, CartesianVector &axisDirection,
    CartesianVector &orthoDirection) const
{
    CartesianPointVector constituentHitSubsetPositionVector(LArHitWidthHelper::GetConstituentHitPositionVector(constituentHitSubsetVector));
 
    if (constituentHitSubsetPositionVector.size() < 2)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    CartesianVector centroid(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVecs;
    LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
    LArPcaHelper::RunPca(constituentHitSubsetPositionVector, centroid, eigenValues, eigenVecs);
    
    axisDirection = eigenVecs.at(0);

    // Use y-axis to generate an orthogonal axis (assuming that cluster occupies X-Z plane)
    const CartesianVector yAxis(0.f, 1.f, 0.f);
    orthoDirection = yAxis.GetCrossProduct(axisDirection).GetUnitVector();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitWidthClusterMergingAlgorithm::GetFittingCoordinates(const CartesianVector &axisDirection, const CartesianVector &constituentHitPosition, float &rL, float &rT) const
{
    const CartesianVector xAxis(1.f, 0.f, 0.f);
    const float openingAngle(axisDirection.GetOpeningAngle(xAxis)), c(std::cos(openingAngle)), s(std::sin(openingAngle));

    rL = (c * constituentHitPosition.GetX()) + (s * constituentHitPosition.GetZ());
    rT = (c * constituentHitPosition.GetZ()) - (s * constituentHitPosition.GetX()); 
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitWidthClusterMergingAlgorithm::GetGlobalDirection(const CartesianVector &axisDirection, const float gradient, CartesianVector &globalDirection) const
{
    const CartesianVector xAxis(1.f, 0.f, 0.f);
    const float openingAngle(axisDirection.GetOpeningAngle(xAxis)), c(std::cos(openingAngle)), s(std::sin(openingAngle));
    const float deltaL(1.f), deltaT(gradient);

    const float x = (c * deltaL) - (s * deltaT);
    const float z = (c * deltaT) + (s * deltaL);
    
    globalDirection.SetValues(x, 0.f, z);
    globalDirection = globalDirection.GetUnitVector();
}
//------------------------------------------------------------------------------------------------------------------------------------------

void HitWidthClusterMergingAlgorithm::GetWeightedGradient(const LArHitWidthHelper::ConstituentHitVector &constituentHitVector, CartesianVector &direction,
    const CartesianVector &fitReferencePoint, const float fittingWeight) const
{
    float weightSum(0.f), weightedLSum(0.f), weightedTSum(0.f);
    bool isLConstant(true), isTConstant(true);

    // get fitting subset vector
    LArHitWidthHelper::ConstituentHitVector constituentHitSubsetVector;
    this->GetConstituentHitSubsetVector(constituentHitVector, fitReferencePoint, fittingWeight, constituentHitSubsetVector);

    // determine the fitting axes
    CartesianVector axisDirection(0.f, 0.f, 0.f), orthoDirection(0.f, 0.f, 0.f);
    this->GetFittingAxes(constituentHitSubsetVector, axisDirection, orthoDirection);

    // to check if cosntituent hit positions are constant in rL or constant in rT (would lead to a division by zero)
    float firstHitL(0.f), firstHitT(0.f);
    this->GetFittingCoordinates(axisDirection, constituentHitSubsetVector.begin()->GetPositionVector(), firstHitL, firstHitT);
        
    for (const LArHitWidthHelper::ConstituentHit &constituentHit : constituentHitSubsetVector) 
    {
        const float hitWeight(constituentHit.GetHitWidth());
        float rL(0.f), rT(0.f);   
        this->GetFittingCoordinates(axisDirection, constituentHit.GetPositionVector(), rL, rT);

        if (std::fabs(firstHitL - rL) > std::numeric_limits<float>::epsilon())
            isLConstant = false;

        if (std::fabs(firstHitT - rT) > std::numeric_limits<float>::epsilon())
            isTConstant = false;

        weightedLSum += rL * hitWeight;
        weightedTSum += rT * hitWeight;
        weightSum += hitWeight;
    }

    if (weightSum < std::numeric_limits<float>::epsilon())
    {
        std::cout << "HitWidthClusterMergingAlgorithm::GetWeightedGradient - hit weight in fit is negative or equivalent to zero" << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);
    }

    // return ortho direction for a cluster with constant rL
    if (isLConstant) 
    {
        direction = orthoDirection;
        return;
    }

    // return axis direction for a cluster with constant rT
    if (isTConstant) 
    {
        direction = axisDirection;
        return;
    }

    const float weightedLMean(weightedLSum/weightSum), weightedTMean(weightedTSum/weightSum);
    float numerator(0.f), denominator(0.f);
    
    for (const LArHitWidthHelper::ConstituentHit &constituentHit : constituentHitSubsetVector) 
    {
        const float hitWeight(constituentHit.GetHitWidth());
        float rL(0.f), rT(0.f);   
        this->GetFittingCoordinates(axisDirection, constituentHit.GetPositionVector(), rL, rT);
        
        numerator += hitWeight * (rL - weightedLMean) * (rT - weightedTMean);
        denominator += hitWeight * pow(rL - weightedLMean, 2);
    }

    const float gradient(numerator/denominator);

    // change coordinates to z=mx+c fit and normalise
    this->GetGlobalDirection(axisDirection, gradient, direction);

    // ATTN fitting convention is to point in direction of increasing x
    if (direction.GetX() < 0.f)
        direction *= -1.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitWidthClusterMergingAlgorithm::RemoveShortcutAssociations(const ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const
{
    // Create temporary map so can delete elements whilst still iterating over them
    ClusterAssociationMap tempMap(clusterAssociationMap);

    for (const Cluster *const pCluster : clusterVector)
    {
        const ClusterAssociationMap::const_iterator primaryMapIter = clusterAssociationMap.find(pCluster);

        if (primaryMapIter == clusterAssociationMap.end())
            continue;

        // put ClusterSet into ClusterVector
        ClusterVector primaryForwardAssociations(primaryMapIter->second.m_forwardAssociations.begin(), primaryMapIter->second.m_forwardAssociations.end());
        std::sort(primaryForwardAssociations.begin(), primaryForwardAssociations.end(), LArClusterHelper::SortByNHits);
        
        // remove primary clusters that are present in secondary associations of other primary clusters
        for (const Cluster *const pConsideredCluster : primaryForwardAssociations)
        {
            for (const Cluster *const pPrimaryCluster : primaryForwardAssociations)
            {
                if (pConsideredCluster == pPrimaryCluster)
                    continue;

                const ClusterAssociationMap::const_iterator secondaryMapIter = clusterAssociationMap.find(pPrimaryCluster);

                // if primary cluster has no associations (this shouldn't ever be the case)
                if (secondaryMapIter == clusterAssociationMap.end())
                    continue;

                const ClusterSet &secondaryForwardAssociations(secondaryMapIter->second.m_forwardAssociations);

                if (secondaryForwardAssociations.find(pConsideredCluster) != secondaryForwardAssociations.end())
                {
                    ClusterSet &tempPrimaryForwardAssociations(tempMap.find(pCluster)->second.m_forwardAssociations);
                    const ClusterSet::const_iterator forwardAssociationToRemove(tempPrimaryForwardAssociations.find(pConsideredCluster));

                    // if association has already been removed
                    if(forwardAssociationToRemove == tempPrimaryForwardAssociations.end())
                        continue;

                    ClusterSet &tempPrimaryBackwardAssociations(tempMap.find(pConsideredCluster)->second.m_backwardAssociations);
                    const ClusterSet::const_iterator backwardAssociationToRemove(tempPrimaryBackwardAssociations.find(pCluster));

                    // if association has already been removed
                    if (backwardAssociationToRemove == tempPrimaryBackwardAssociations.end())
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

StatusCode HitWidthClusterMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FittingWeight", m_fittingWeight));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterWeight", m_minClusterWeight));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxXMergeDistance", m_maxXMergeDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxZMergeDistance", m_maxZMergeDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMergeCosOpeningAngle", m_minMergeCosOpeningAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinDirectionDeviationCosAngle", m_minDirectionDeviationCosAngle));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxConstituentHitWidth", m_maxConstituentHitWidth));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "HitWidthScalingFactor", m_hitWidthScalingFactor));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "UseOldDirectionMethod", m_useOldDirectionMethod));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "MinClusterSparseness", m_minClusterSparseness));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "UseClosestMergePoint", m_useClosestMergePoint));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "DoubleFittingWeight", m_doubleFittingWeight));

    return ClusterAssociationAlgorithm::ReadSettings(xmlHandle);
}
  
} // namespace lar_content
