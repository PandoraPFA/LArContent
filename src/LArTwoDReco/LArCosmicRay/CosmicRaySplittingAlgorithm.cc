/**
 *  @file   LArContent/src/LArTwoDReco/LArCosmicRay/CosmicRaySplittingAlgorithm.cc
 *
 *  @brief  Implementation of the cosmic ray splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArTwoDReco/LArCosmicRay/CosmicRaySplittingAlgorithm.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"
#include "LArHelpers/LArPointingClusterHelper.h"

#include "LArPlugins/LArTransformationPlugin.h"

using namespace pandora;

namespace lar_content
{

CosmicRaySplittingAlgorithm::CosmicRaySplittingAlgorithm() :
    m_clusterMinLength(10.f),
    m_halfWindowLayers(30),
    m_samplingPitch(1.f),
    m_maxCosSplittingAngle(0.9925f),
    m_minCosMergingAngle(0.94f),
    m_maxTransverseDisplacement(5.f),
    m_maxLongitudinalDisplacement(30.f),
    m_maxLongitudinalDisplacementSquared(m_maxLongitudinalDisplacement * m_maxLongitudinalDisplacement)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRaySplittingAlgorithm::Run()
{
    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

    // Get ordered list of clean clusters
    ClusterVector clusterVector;
    this->GetListOfCleanClusters(pClusterList, clusterVector);

    // Calculate sliding fit results for clean clusters
    TwoDSlidingFitResultMap slidingFitResultMap;
    this->BuildSlidingFitResultMap(clusterVector, slidingFitResultMap);

    // Loop over clusters, identify and perform splits
    ClusterList splitClusterList;

    for (ClusterVector::const_iterator bIter = clusterVector.begin(), bIterEnd1 = clusterVector.end(); bIter != bIterEnd1; ++bIter)
    {
        if (splitClusterList.count(*bIter) > 0)
            continue;

        TwoDSlidingFitResultMap::const_iterator bFitIter = slidingFitResultMap.find(*bIter);

        if (slidingFitResultMap.end() == bFitIter)
            continue;

        const TwoDSlidingFitResult &branchSlidingFitResult(bFitIter->second);

        // Find best split position for candidate branch cluster
        CartesianVector splitPosition(0.f,0.f,0.f);
        CartesianVector splitDirection1(0.f,0.f,0.f);
        CartesianVector splitDirection2(0.f,0.f,0.f);

        if (STATUS_CODE_SUCCESS != this->FindBestSplitPosition(branchSlidingFitResult, splitPosition, splitDirection1, splitDirection2))
            continue;

        // Find candidate replacement clusters to merge into branch cluster at the split position
        TwoDSlidingFitResultMap::const_iterator bestReplacementIter1(slidingFitResultMap.end());
        TwoDSlidingFitResultMap::const_iterator bestReplacementIter2(slidingFitResultMap.end());

        float bestLengthSquared1(m_maxLongitudinalDisplacementSquared);
        float bestLengthSquared2(m_maxLongitudinalDisplacementSquared);

        for (ClusterVector::const_iterator rIter = clusterVector.begin(), rIterEnd = clusterVector.end(); rIter != rIterEnd; ++rIter)
        {
            if (splitClusterList.count(*rIter) > 0)
                continue;

            TwoDSlidingFitResultMap::const_iterator rFitIter = slidingFitResultMap.find(*rIter);

            if (slidingFitResultMap.end() == rFitIter)
                continue;

            const TwoDSlidingFitResult &replacementSlidingFitResult(rFitIter->second);

            if (branchSlidingFitResult.GetCluster() == replacementSlidingFitResult.GetCluster())
                continue;

            float lengthSquared1(std::numeric_limits<float>::max());
            float lengthSquared2(std::numeric_limits<float>::max());

            if (STATUS_CODE_SUCCESS != this->ConfirmSplitPosition(branchSlidingFitResult, replacementSlidingFitResult,
                splitPosition, splitDirection1, splitDirection2, lengthSquared1, lengthSquared2))
                continue;

            if (lengthSquared1 < bestLengthSquared1)
            {
                bestLengthSquared1 = lengthSquared1;
                bestReplacementIter1 = rFitIter;
            }

            if (lengthSquared2 < bestLengthSquared2)
            {
                bestLengthSquared2 = lengthSquared2;
                bestReplacementIter2 = rFitIter;
            }
        }

        Cluster *pReplacementCluster1 = NULL;
        Cluster *pReplacementCluster2 = NULL;

        if (slidingFitResultMap.end() != bestReplacementIter1)
            pReplacementCluster1 = const_cast<Cluster*>(bestReplacementIter1->first);

        if (slidingFitResultMap.end() != bestReplacementIter2)
            pReplacementCluster2 = const_cast<Cluster*>(bestReplacementIter2->first);

        if (NULL == pReplacementCluster1 && NULL == pReplacementCluster2)
            continue;

        Cluster *pBranchCluster = const_cast<Cluster*>(branchSlidingFitResult.GetCluster());

        // Crossed tracks
        if (pReplacementCluster1 && pReplacementCluster2)
        {
            if (!(this->IdentifyCrossedTracks(pBranchCluster, pReplacementCluster1, pReplacementCluster2, splitPosition)))
                continue;

            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->PerformDoubleSplit(pBranchCluster, pReplacementCluster1,
                pReplacementCluster2, splitPosition, splitDirection1, splitDirection2));
        }
        // Single scatter
        else if (pReplacementCluster1)
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->PerformSingleSplit(pBranchCluster, pReplacementCluster1,
                splitPosition, splitDirection1, splitDirection2));
        }
        else if (pReplacementCluster2)
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->PerformSingleSplit(pBranchCluster, pReplacementCluster2,
                splitPosition, splitDirection2, splitDirection1));
        }

        // Choose not to re-use clusters (for now)
        if (pReplacementCluster1)
            splitClusterList.insert(pReplacementCluster1);

        if (pReplacementCluster2)
            splitClusterList.insert(pReplacementCluster2);

        splitClusterList.insert(pBranchCluster);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRaySplittingAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (LArClusterHelper::GetLengthSquared(pCluster) < m_clusterMinLength * m_clusterMinLength)
            continue;

        clusterVector.push_back(pCluster);
    }

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRaySplittingAlgorithm::BuildSlidingFitResultMap(const ClusterVector &clusterVector,
    TwoDSlidingFitResultMap &slidingFitResultMap) const
{
    const float slidingFitPitch(LArGeometryHelper::GetLArTransformationPlugin(this->GetPandora())->GetWireZPitch());

    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
    {
        if (slidingFitResultMap.end() == slidingFitResultMap.find(*iter))
        {
            try
            {
                const TwoDSlidingFitResult slidingFitResult(*iter, m_halfWindowLayers, slidingFitPitch);

                if (!slidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(*iter, slidingFitResult)).second)
                    throw StatusCodeException(STATUS_CODE_FAILURE);
            }
            catch (StatusCodeException &statusCodeException)
            {
                if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                    throw statusCodeException;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRaySplittingAlgorithm::FindBestSplitPosition(const TwoDSlidingFitResult &branchSlidingFitResult, CartesianVector &splitPosition,
    CartesianVector &splitDirection1, CartesianVector &splitDirection2) const
{
    // Find position of greatest scatter for this cluster
    float splitCosTheta(m_maxCosSplittingAngle);
    bool foundSplit(false);

    const CartesianVector &minPosition(branchSlidingFitResult.GetGlobalMinLayerPosition());
    const CartesianVector &maxPosition(branchSlidingFitResult.GetGlobalMaxLayerPosition());
    const float halfWindowLength(branchSlidingFitResult.GetLayerFitHalfWindowLength());

    float minL(0.f), maxL(0.f), minT(0.f), maxT(0.f);
    branchSlidingFitResult.GetLocalPosition(minPosition, minL, minT);
    branchSlidingFitResult.GetLocalPosition(maxPosition, maxL, maxT);

    const unsigned int nSamplingPoints = static_cast<unsigned int>((maxL - minL)/ m_samplingPitch);

    for (unsigned int n = 0; n < nSamplingPoints; ++n)
    {
        const float alpha((0.5f + static_cast<float>(n)) / static_cast<float>(nSamplingPoints));
        const float rL(minL + (maxL - minL) * alpha);

        try
        {
            CartesianVector centralPosition(0.f,0.f,0.f);
            CartesianVector forwardDirection(0.f,0.f,0.f);
            CartesianVector backwardDirection(0.f,0.f,0.f);

            branchSlidingFitResult.GetGlobalFitPosition(rL, centralPosition);
            branchSlidingFitResult.GetGlobalFitDirection(rL + halfWindowLength, forwardDirection);
            branchSlidingFitResult.GetGlobalFitDirection(rL - halfWindowLength, backwardDirection);

            const float cosTheta(forwardDirection.GetDotProduct(backwardDirection));

            if (cosTheta < splitCosTheta)
            {
                CartesianVector forwardPosition(0.f,0.f,0.f);
                CartesianVector backwardPosition(0.f,0.f,0.f);

                branchSlidingFitResult.GetGlobalFitPosition(rL + halfWindowLength, forwardPosition);
                branchSlidingFitResult.GetGlobalFitPosition(rL - halfWindowLength, backwardPosition);

                CartesianVector forwardVectorOutwards(forwardPosition - centralPosition);
                CartesianVector backwardVectorOutwards(backwardPosition - centralPosition);

                splitPosition = centralPosition;
                splitDirection1 = (forwardDirection.GetDotProduct(forwardVectorOutwards) > 0.f) ? forwardDirection : forwardDirection * -1.f;
                splitDirection2 = (backwardDirection.GetDotProduct(backwardVectorOutwards) > 0.f) ? backwardDirection : backwardDirection * -1.f;
                splitCosTheta = cosTheta;
                foundSplit = true;
            }
        }
        catch (StatusCodeException &)
        {

        }
    }

    if (!foundSplit)
        return STATUS_CODE_NOT_FOUND;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRaySplittingAlgorithm::ConfirmSplitPosition(const TwoDSlidingFitResult &branchSlidingFitResult,
    const TwoDSlidingFitResult &replacementSlidingFitResult, const CartesianVector &splitPosition, const CartesianVector &splitDirection1,
    const CartesianVector &splitDirection2, float &bestLengthSquared1, float &bestLengthSquared2) const
{
    // Check if the replacement cluster points to the split position on the branch cluster
    bestLengthSquared1 = std::numeric_limits<float>::max();
    bestLengthSquared2 = std::numeric_limits<float>::max();

    bool foundMatch(false);

    const CartesianVector minPosition(replacementSlidingFitResult.GetGlobalMinLayerPosition());
    const CartesianVector maxPosition(replacementSlidingFitResult.GetGlobalMaxLayerPosition());
    const CartesianVector minDirection(replacementSlidingFitResult.GetGlobalMinLayerDirection());
    const CartesianVector maxDirection(replacementSlidingFitResult.GetGlobalMaxLayerDirection());

    const CartesianVector branchVertex1(branchSlidingFitResult.GetGlobalMinLayerPosition());
    const CartesianVector branchVertex2(branchSlidingFitResult.GetGlobalMaxLayerPosition());
    const CartesianVector branchDirection1(splitDirection1 * -1.f);
    const CartesianVector branchDirection2(splitDirection2 * -1.f);

    const float cosSplittingAngle(-splitDirection1.GetDotProduct(splitDirection2));
    const float branchLengthSquared((branchVertex2 - branchVertex1).GetMagnitudeSquared());
    const float replacementLengthSquared((maxPosition - minPosition).GetMagnitudeSquared());

    // Loop over each end of the replacement cluster
    for (unsigned int iFwd = 0; iFwd < 2; ++iFwd)
    {
        const CartesianVector pAxis((0 == iFwd) ? (maxPosition - minPosition) : (minPosition - maxPosition));
        const CartesianVector vtxPosition((0 == iFwd) ? minPosition : maxPosition);
        const CartesianVector endPosition((0 == iFwd) ? maxPosition : minPosition);
        const CartesianVector vtxDirection((0 == iFwd) ? (pAxis.GetDotProduct(minDirection) > 0 ? minDirection : minDirection * -1.f) :
            (pAxis.GetDotProduct(maxDirection) > 0 ? maxDirection : maxDirection * -1.f));

        // Choose the correct end of the replacement cluster and require proximity to the branch cluster
        const float vtxDisplacement(LArClusterHelper::GetClosestDistance(vtxPosition, branchSlidingFitResult.GetCluster()));
        const float endDisplacement(LArClusterHelper::GetClosestDistance(endPosition, branchSlidingFitResult.GetCluster()));

        const float lengthSquared((vtxPosition - splitPosition).GetMagnitudeSquared());
        const float lengthSquared1((vtxPosition - branchVertex1).GetMagnitudeSquared());
        const float lengthSquared2((vtxPosition - branchVertex2).GetMagnitudeSquared());

        if (vtxDisplacement > endDisplacement)
            continue;

        if (std::min(lengthSquared,std::min(lengthSquared1,lengthSquared2)) > std::min(branchLengthSquared, replacementLengthSquared))
            continue;

        // Require good pointing information between replacement cluster and candidate split position
        float impactL(0.f), impactT(0.f);
        LArPointingClusterHelper::GetImpactParameters(vtxPosition, vtxDirection, splitPosition, impactL, impactT);

        if (impactT > m_maxTransverseDisplacement || impactL > m_maxLongitudinalDisplacement ||
            impactL < -1.f || impactT > std::max(1.5f, 0.577f * (1.f + impactL)))
            continue;

        // Check the segment of the branch cluster above the split position
        if (vtxDirection.GetDotProduct(branchDirection1) > std::max(m_minCosMergingAngle,cosSplittingAngle))
        {
            if (lengthSquared < bestLengthSquared1)
            {
                bestLengthSquared1 = lengthSquared;
                foundMatch = true;
            }
        }

        // Check the segment of the branch cluster below the split position
        if (vtxDirection.GetDotProduct(branchDirection2) > std::max(m_minCosMergingAngle,cosSplittingAngle))
        {
            if (lengthSquared < bestLengthSquared2)
            {
                bestLengthSquared2 = lengthSquared;
                foundMatch = true;
            }
        }
    }

    if (!foundMatch)
        return STATUS_CODE_NOT_FOUND;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRaySplittingAlgorithm::PerformSingleSplit(Cluster *pBranchCluster, Cluster *pReplacementCluster,
    const CartesianVector &splitPosition, const CartesianVector &forwardDirection, const CartesianVector &backwardDirection) const
{
    CaloHitList caloHitListToMove;
    this->GetCaloHitListToMove(pBranchCluster, pReplacementCluster, splitPosition, forwardDirection, backwardDirection, caloHitListToMove);

    CaloHitList caloHitListToKeep;
    this->GetCaloHitListToKeep(pBranchCluster, caloHitListToMove, caloHitListToKeep);

    if (caloHitListToKeep.empty())
        return PandoraContentApi::MergeAndDeleteClusters(*this, pReplacementCluster, pBranchCluster);

    return this->SplitCluster(pBranchCluster, pReplacementCluster, caloHitListToMove);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRaySplittingAlgorithm::PerformDoubleSplit(Cluster *pBranchCluster, Cluster *pReplacementCluster1, Cluster *pReplacementCluster2,
    const CartesianVector &splitPosition, const CartesianVector &splitDirection1, const CartesianVector &splitDirection2) const
{
    CaloHitList caloHitListToMove1, caloHitListToMove2;
    this->GetCaloHitListsToMove(pBranchCluster, splitPosition, splitDirection1, splitDirection2, caloHitListToMove1, caloHitListToMove2);

    CaloHitList caloHitListToKeep1;
    this->GetCaloHitListToKeep(pBranchCluster, caloHitListToMove1, caloHitListToKeep1);

    if (caloHitListToKeep1.empty())
        return STATUS_CODE_FAILURE;

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->SplitCluster(pBranchCluster, pReplacementCluster1, caloHitListToMove1));

    CaloHitList caloHitListToKeep2;
    this->GetCaloHitListToKeep(pBranchCluster, caloHitListToMove2, caloHitListToKeep2);

    if (caloHitListToKeep2.empty())
        return PandoraContentApi::MergeAndDeleteClusters(*this, pReplacementCluster2, pBranchCluster);

    return this->SplitCluster(pBranchCluster, pReplacementCluster2, caloHitListToMove2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRaySplittingAlgorithm::GetCaloHitListToMove(const Cluster *const pBranchCluster, const Cluster *const pReplacementCluster,
    const CartesianVector &splitPosition, const CartesianVector &forwardDirection, const CartesianVector &backwardDirection,
    CaloHitList &caloHitListToMove) const
{
    const CartesianVector forwardPosition(LArClusterHelper::GetClosestPosition(splitPosition,pBranchCluster));
    const CartesianVector vtxPosition(LArClusterHelper::GetClosestPosition(splitPosition,pReplacementCluster));
    const CartesianVector vtxDirection((forwardPosition - vtxPosition).GetUnitVector());

    CaloHitList caloHitListToSort, caloHitListToCheck;
    pBranchCluster->GetOrderedCaloHitList().GetCaloHitList(caloHitListToSort);

    for (CaloHitList::const_iterator iter = caloHitListToSort.begin(), iterEnd = caloHitListToSort.end(); iter != iterEnd; ++iter)
    {
        CaloHit *pCaloHit = *iter;

        if (forwardDirection.GetDotProduct(pCaloHit->GetPositionVector() - forwardPosition) > 0.f)
        {
            caloHitListToMove.insert(pCaloHit);
        }
        else if(forwardDirection.GetDotProduct(pCaloHit->GetPositionVector() - vtxPosition) > -1.25f)
        {
            caloHitListToCheck.insert(pCaloHit);
        }
    }

    float closestLengthSquared(std::numeric_limits<float>::max());

    for (CaloHitList::const_iterator iter = caloHitListToCheck.begin(), iterEnd = caloHitListToCheck.end(); iter != iterEnd; ++iter)
    {
        CaloHit *pCaloHit = *iter;

        if (vtxDirection.GetCrossProduct(pCaloHit->GetPositionVector() - forwardPosition).GetMagnitude() >
            backwardDirection.GetCrossProduct(pCaloHit->GetPositionVector() - forwardPosition).GetMagnitude())
            continue;

        const float lengthSquared((pCaloHit->GetPositionVector() - vtxPosition).GetMagnitudeSquared());

        if (lengthSquared < closestLengthSquared)
            closestLengthSquared = lengthSquared;
    }

    for (CaloHitList::const_iterator iter = caloHitListToCheck.begin(), iterEnd = caloHitListToCheck.end(); iter != iterEnd; ++iter)
    {
        CaloHit *pCaloHit = *iter;

        if ((pCaloHit->GetPositionVector() - vtxPosition).GetMagnitudeSquared() >= closestLengthSquared)
            caloHitListToMove.insert(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRaySplittingAlgorithm::GetCaloHitListsToMove(const Cluster *const pBranchCluster, const CartesianVector &splitPosition,
    const CartesianVector &splitDirection1, const CartesianVector &splitDirection2, CaloHitList &caloHitListToMove1, CaloHitList &caloHitListToMove2) const
{
    CaloHitList caloHitListToSort;
    pBranchCluster->GetOrderedCaloHitList().GetCaloHitList(caloHitListToSort);

    for (CaloHitList::const_iterator iter = caloHitListToSort.begin(), iterEnd = caloHitListToSort.end(); iter != iterEnd; ++iter)
    {
        CaloHit *pCaloHit = *iter;

        const float displacement1(splitDirection1.GetDotProduct(pCaloHit->GetPositionVector() - splitPosition));
        const float displacement2(splitDirection2.GetDotProduct(pCaloHit->GetPositionVector() - splitPosition));

        if (displacement1 > displacement2)
        {
            caloHitListToMove1.insert(pCaloHit);
        }
        else
        {
            caloHitListToMove2.insert(pCaloHit);
        }
    }
}
//------------------------------------------------------------------------------------------------------------------------------------------

bool CosmicRaySplittingAlgorithm::IdentifyCrossedTracks(const Cluster *const pBranchCluster, const Cluster *const pReplacementCluster1,
    const Cluster *const pReplacementCluster2, const pandora::CartesianVector &splitPosition) const
{
    CartesianVector branchVertex1(0.f,0.f,0.f), branchVertex2(0.f,0.f,0.f);
    LArClusterHelper::GetExtremalCoordinates(pBranchCluster,branchVertex1,branchVertex2);

    const CartesianVector replacementVertex1(LArClusterHelper::GetClosestPosition(splitPosition,pReplacementCluster1));
    const CartesianVector replacementVertex2(LArClusterHelper::GetClosestPosition(splitPosition,pReplacementCluster2));

    if ((replacementVertex2 - replacementVertex1).GetMagnitudeSquared() > (branchVertex2 - branchVertex1).GetMagnitudeSquared())
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRaySplittingAlgorithm::GetCaloHitListToKeep(const Cluster *const pBranchCluster, const CaloHitList &caloHitListToMove,
    CaloHitList &caloHitListToKeep) const
{
    if (caloHitListToMove.empty())
        return STATUS_CODE_FAILURE;

    CaloHitList caloHitList;
    pBranchCluster->GetOrderedCaloHitList().GetCaloHitList(caloHitList);

    for (CaloHitList::const_iterator iter = caloHitList.begin(), iterEnd = caloHitList.end(); iter != iterEnd; ++iter)
    {
        CaloHit *pCaloHit = *iter;

        if (!caloHitListToMove.count(pCaloHit))
            caloHitListToKeep.insert(pCaloHit);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRaySplittingAlgorithm::SplitCluster(Cluster *pBranchCluster, Cluster *pReplacementCluster, const CaloHitList &caloHitListToMove) const
{
    if (caloHitListToMove.empty())
        return STATUS_CODE_FAILURE;

    for (CaloHitList::const_iterator iter = caloHitListToMove.begin(), iterEnd = caloHitListToMove.end(); iter != iterEnd; ++iter)
    {
        CaloHit *pCaloHit = *iter;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromCluster(*this, pBranchCluster, pCaloHit));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pReplacementCluster, pCaloHit));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRaySplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterMinLength", m_clusterMinLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitHalfWindow", m_halfWindowLayers));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SamplingPitch", m_samplingPitch));

    if (m_samplingPitch < std::numeric_limits<float>::epsilon())
        return STATUS_CODE_INVALID_PARAMETER;

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxCosSplittingAngle", m_maxCosSplittingAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCosMergingAngle", m_minCosMergingAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTransverseDisplacement", m_maxTransverseDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxLongitudinalDisplacement", m_maxLongitudinalDisplacement));
    m_maxLongitudinalDisplacementSquared = m_maxLongitudinalDisplacement * m_maxLongitudinalDisplacement;

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
