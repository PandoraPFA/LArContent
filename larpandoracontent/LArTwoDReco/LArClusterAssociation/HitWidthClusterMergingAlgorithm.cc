/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/HitWidthClusterMergingAlgorithm.cc
 *
 *  @brief  Implementation of the hit width cluster merging algorithm class.
 *
 *  $Log: $
 */
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/HitWidthClusterMergingAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"

using namespace pandora;

namespace lar_content
{

HitWidthClusterMergingAlgorithm::HitWidthClusterMergingAlgorithm() :
    m_maxConstituentHitWidth(0.5f),
    m_hitWidthScalingFactor(1.f),
    m_fittingWeight(20.f),
    m_minClusterWeight(0.5f),
    m_maxXMergeDistance(5.f),
    m_maxZMergeDistance(2.f),
    m_minMergeCosOpeningAngle(0.97f),
    m_minDirectionDeviationCosAngle(0.9f),
    m_minClusterSparseness(0.3f)
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

        const unsigned int numberOfProposedConstituentHits(
            LArHitWidthHelper::GetNProposedConstituentHits(pCluster, m_maxConstituentHitWidth, m_hitWidthScalingFactor));

        if (numberOfProposedConstituentHits == 0)
            continue;

        // clusterSparseness [0 -> 1] where a higher value indicates sparseness
        const float clusterSparseness(1.f - (static_cast<float>(pCluster->GetNCaloHits()) / static_cast<float>(numberOfProposedConstituentHits)));

        if (clusterSparseness < m_minClusterSparseness)
            continue;

        m_clusterToParametersMap.insert(std::pair<const Cluster *, LArHitWidthHelper::ClusterParameters>(
            pCluster, LArHitWidthHelper::ClusterParameters(pCluster, m_maxConstituentHitWidth, false, m_hitWidthScalingFactor)));

        clusterVector.push_back(pCluster);
    }

    std::sort(clusterVector.begin(), clusterVector.end(), LArHitWidthHelper::SortByHigherXExtrema(m_clusterToParametersMap));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitWidthClusterMergingAlgorithm::PopulateClusterAssociationMap(const ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const
{
    // ATTN this method assumes that clusters have been sorted by extremal x position (low higherXExtrema -> high higherXExtrema)
    for (ClusterVector::const_iterator iterCurrentCluster = clusterVector.begin(); iterCurrentCluster != clusterVector.end(); ++iterCurrentCluster)
    {
        const Cluster *const pCurrentCluster = *iterCurrentCluster;
        const LArHitWidthHelper::ClusterParameters &currentClusterParameters(
            LArHitWidthHelper::GetClusterParameters(pCurrentCluster, m_clusterToParametersMap));

        for (ClusterVector::const_iterator iterTestCluster = iterCurrentCluster; iterTestCluster != clusterVector.end(); ++iterTestCluster)
        {
            if (iterCurrentCluster == iterTestCluster)
                continue;

            const Cluster *const pTestCluster = *iterTestCluster;
            const LArHitWidthHelper::ClusterParameters &testClusterParameters(
                LArHitWidthHelper::GetClusterParameters(pTestCluster, m_clusterToParametersMap));

            if (!this->AreClustersAssociated(currentClusterParameters, testClusterParameters))
                continue;

            clusterAssociationMap[pCurrentCluster].m_forwardAssociations.insert(pTestCluster);
            clusterAssociationMap[pTestCluster].m_backwardAssociations.insert(pCurrentCluster);
        }
    }

    this->RemoveShortcutAssociations(clusterVector, clusterAssociationMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool HitWidthClusterMergingAlgorithm::IsExtremalCluster(const bool isForward, const Cluster *const pCurrentCluster, const Cluster *const pTestCluster) const
{
    //ATTN - cannot use map since higherXExtrema may have changed during merging
    const LArHitWidthHelper::ConstituentHitVector currentConstituentHitVector(
        LArHitWidthHelper::GetConstituentHits(pCurrentCluster, m_maxConstituentHitWidth, m_hitWidthScalingFactor, false));
    const LArHitWidthHelper::ConstituentHitVector testConstituentHitVector(
        LArHitWidthHelper::GetConstituentHits(pTestCluster, m_maxConstituentHitWidth, m_hitWidthScalingFactor, false));
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

bool HitWidthClusterMergingAlgorithm::AreClustersAssociated(
    const LArHitWidthHelper::ClusterParameters &currentFitParameters, const LArHitWidthHelper::ClusterParameters &testFitParameters) const
{
    // check cluster extrema are not too far away in x
    if (testFitParameters.GetLowerXExtrema().GetX() > (currentFitParameters.GetHigherXExtrema().GetX() + m_maxXMergeDistance))
        return false;

    // check cluster extrema are not too far away in z
    if (testFitParameters.GetLowerXExtrema().GetZ() > (currentFitParameters.GetHigherXExtrema().GetZ() + m_maxZMergeDistance) ||
        testFitParameters.GetLowerXExtrema().GetZ() < (currentFitParameters.GetHigherXExtrema().GetZ() - m_maxZMergeDistance))
    {
        return false;
    }

    // find appropriate test merging point (may be different to lowerXExtrema when hits overlap)
    CartesianVector testMergePoint(0.f, 0.f, 0.f);
    if (testFitParameters.GetLowerXExtrema().GetX() < currentFitParameters.GetHigherXExtrema().GetX())
    {
        this->FindClosestPointToPosition(currentFitParameters.GetHigherXExtrema(), testFitParameters.GetConstituentHitVector(), testMergePoint);

        // check closeness in z is maintained
        if (testMergePoint.GetZ() > (currentFitParameters.GetHigherXExtrema().GetZ() + m_maxZMergeDistance) ||
            testMergePoint.GetZ() < (currentFitParameters.GetHigherXExtrema().GetZ() - m_maxZMergeDistance))
        {
            return false;
        }
    }
    else
    {
        testMergePoint = testFitParameters.GetLowerXExtrema();
    }

    CartesianVector currentClusterDirection(0.f, 0.f, 0.f), testClusterDirection(0.f, 0.f, 0.f);

    try
    {
        this->GetClusterDirection(currentFitParameters.GetConstituentHitVector(), currentClusterDirection,
            currentFitParameters.GetHigherXExtrema(), m_fittingWeight);
        this->GetClusterDirection(testFitParameters.GetConstituentHitVector(), testClusterDirection, testMergePoint, m_fittingWeight);
    }
    catch (const StatusCodeException &)
    {
        return false;
    }

    // check clusters have a similar direction
    if (currentClusterDirection.GetCosOpeningAngle(testClusterDirection) < m_minMergeCosOpeningAngle)
        return false;

    // check that the new direction is consistent with the old clusters
    LArHitWidthHelper::ConstituentHitVector newConstituentHitVector(currentFitParameters.GetConstituentHitVector());
    newConstituentHitVector.insert(newConstituentHitVector.end(), testFitParameters.GetConstituentHitVector().begin(),
        testFitParameters.GetConstituentHitVector().end());

    const CartesianVector midpoint((currentFitParameters.GetHigherXExtrema() + testMergePoint) * 0.5);
    CartesianVector newClusterDirection(0.f, 0.f, 0.f);

    try
    {
        this->GetClusterDirection(newConstituentHitVector, newClusterDirection, midpoint, m_fittingWeight);
    }
    catch (const StatusCodeException &)
    {
        return false;
    }

    if (newClusterDirection.GetCosOpeningAngle(currentClusterDirection) < m_minDirectionDeviationCosAngle ||
        newClusterDirection.GetCosOpeningAngle(testClusterDirection) < m_minDirectionDeviationCosAngle)
    {
        return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitWidthClusterMergingAlgorithm::FindClosestPointToPosition(const CartesianVector &position,
    const LArHitWidthHelper::ConstituentHitVector &constituentHitVector, CartesianVector &closestPoint) const
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

void HitWidthClusterMergingAlgorithm::GetClusterDirection(const LArHitWidthHelper::ConstituentHitVector &constituentHitVector,
    CartesianVector &direction, const CartesianVector &fitReferencePoint, const float fittingWeight) const
{
    float weightSum(0.f), weightedLSum(0.f), weightedTSum(0.f);
    bool isLConstant(true), isTConstant(true);

    // get fitting subset vector
    LArHitWidthHelper::ConstituentHitVector constituentHitSubsetVector;
    this->GetConstituentHitSubsetVector(constituentHitVector, fitReferencePoint, fittingWeight, constituentHitSubsetVector);

    // determine the fitting axes
    CartesianVector axisDirection(0.f, 0.f, 0.f), orthoDirection(0.f, 0.f, 0.f);
    this->GetFittingAxes(constituentHitSubsetVector, axisDirection, orthoDirection);

    // the comparison to check if constituent hit positions are constant in rL or rT (would lead to a division by zero)
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
        // ATTN direction convention is to point in direction of increasing x
        direction = (orthoDirection.GetX() < 0.f) ? orthoDirection * (-1.f) : orthoDirection;
        return;
    }

    // return axis direction for a cluster with constant rT
    if (isTConstant)
    {
        // ATTN direction convention is to point in direction of increasing x
        direction = (axisDirection.GetX() < 0.f) ? axisDirection * (-1.f) : axisDirection;
        return;
    }

    const float weightedLMean(weightedLSum / weightSum), weightedTMean(weightedTSum / weightSum);
    float numerator(0.f), denominator(0.f);

    for (const LArHitWidthHelper::ConstituentHit &constituentHit : constituentHitSubsetVector)
    {
        const float hitWeight(constituentHit.GetHitWidth());
        float rL(0.f), rT(0.f);
        this->GetFittingCoordinates(axisDirection, constituentHit.GetPositionVector(), rL, rT);

        numerator += hitWeight * (rL - weightedLMean) * (rT - weightedTMean);
        denominator += hitWeight * pow(rL - weightedLMean, 2);
    }

    if (denominator < std::numeric_limits<float>::epsilon())
    {
        std::cout << "HitWidthClusterMergingAlgorithm::GetWeightedGradient - denominator is zero" << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);
    }
    const float gradient(numerator / denominator);

    // change coordinates to z=mx+c fit and normalise
    this->GetGlobalDirection(axisDirection, gradient, direction);

    // ATTN direction convention is to point in direction of increasing x
    if (direction.GetX() < 0.f)
        direction *= -1.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitWidthClusterMergingAlgorithm::GetConstituentHitSubsetVector(const LArHitWidthHelper::ConstituentHitVector &constituentHitVector,
    const CartesianVector &fitReferencePoint, const float fittingWeight, LArHitWidthHelper::ConstituentHitVector &constituentHitSubsetVector) const
{
    LArHitWidthHelper::ConstituentHitVector sortedConstituentHitVector(constituentHitVector);

    // sort hits with respect to their distance to the fitReferencePoint (closest -> furthest)
    std::sort(sortedConstituentHitVector.begin(), sortedConstituentHitVector.end(),
        LArHitWidthHelper::ConstituentHit::SortByDistanceToPoint(fitReferencePoint));

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

void HitWidthClusterMergingAlgorithm::GetFittingAxes(const LArHitWidthHelper::ConstituentHitVector &constituentHitSubsetVector,
    CartesianVector &axisDirection, CartesianVector &orthoDirection) const
{
    CartesianPointVector constituentHitSubsetPositionVector(LArHitWidthHelper::GetConstituentHitPositionVector(constituentHitSubsetVector));

    if (constituentHitSubsetPositionVector.size() < 2)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    CartesianVector centroid(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVecs;
    LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
    LArPcaHelper::RunPca(constituentHitSubsetPositionVector, centroid, eigenValues, eigenVecs);

    axisDirection = eigenVecs.at(0);

    // ATTN fitting convention is to point in direction of increasing z
    if (axisDirection.GetZ() < 0.f)
        axisDirection *= -1.f;

    // Use y-axis to generate an orthogonal axis (assuming that cluster occupies x-z plane)
    const CartesianVector yAxis(0.f, 1.f, 0.f);
    orthoDirection = yAxis.GetCrossProduct(axisDirection).GetUnitVector();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitWidthClusterMergingAlgorithm::GetFittingCoordinates(
    const CartesianVector &axisDirection, const CartesianVector &constituentHitPosition, float &rL, float &rT) const
{
    // axisDirection is in positive z by convention so use opening angle to obtain coordinates in rotated frame
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
        ClusterVector primaryForwardAssociations(
            primaryMapIter->second.m_forwardAssociations.begin(), primaryMapIter->second.m_forwardAssociations.end());
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
                    if (tempMap.find(pCluster) == tempMap.end())
                        continue;
                    ClusterSet &tempPrimaryForwardAssociations(tempMap.at(pCluster).m_forwardAssociations);
                    const ClusterSet::const_iterator forwardAssociationToRemove(tempPrimaryForwardAssociations.find(pConsideredCluster));

                    // if association has already been removed
                    if (forwardAssociationToRemove == tempPrimaryForwardAssociations.end())
                        continue;

                    if (tempMap.find(pConsideredCluster) == tempMap.end())
                        continue;
                    ClusterSet &tempPrimaryBackwardAssociations(tempMap.at(pConsideredCluster).m_backwardAssociations);
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

    clusterAssociationMap = std::move(tempMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HitWidthClusterMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FittingWeight", m_fittingWeight));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterWeight", m_minClusterWeight));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxXMergeDistance", m_maxXMergeDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxZMergeDistance", m_maxZMergeDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinMergeCosOpeningAngle", m_minMergeCosOpeningAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinDirectionDeviationCosAngle", m_minDirectionDeviationCosAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxConstituentHitWidth", m_maxConstituentHitWidth));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "HitWidthScalingFactor", m_hitWidthScalingFactor));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterSparseness", m_minClusterSparseness));

    return ClusterAssociationAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
