/**
 *  @file   larpandoracontent/LArThreeDReco/LArTrackFragments/ThreeViewTrackFragmentsAlgorithm.cc
 *
 *  @brief  Implementation of the three view fragments algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArThreeDReco/LArTrackFragments/ThreeViewTrackFragmentsAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ThreeViewTrackFragmentsAlgorithm::ThreeViewTrackFragmentsAlgorithm() :
    m_nMaxTensorToolRepeats(1000),
    m_minXOverlap(3.f),
    m_minXOverlapFraction(0.8f),
    m_maxPointDisplacementSquared(1.5f * 1.5f),
    m_minMatchedSamplingPointFraction(0.5f),
    m_minMatchedHits(5)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewTrackFragmentsAlgorithm::UpdateForNewCluster(const Cluster *const pNewCluster)
{
    try
    {
        this->AddToSlidingFitCache(pNewCluster);
    }
    catch (StatusCodeException &statusCodeException)
    {
        if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
            throw statusCodeException;

        return;
    }

    const HitType hitType(LArClusterHelper::GetClusterHitType(pNewCluster));

    if (!((TPC_VIEW_U == hitType) || (TPC_VIEW_V == hitType) || (TPC_VIEW_W == hitType)))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    ClusterList &clusterList((TPC_VIEW_U == hitType) ? m_clusterListU : (TPC_VIEW_V == hitType) ? m_clusterListV : m_clusterListW);

    if (clusterList.end() != std::find(clusterList.begin(), clusterList.end(), pNewCluster))
        throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);

    clusterList.push_back(pNewCluster);

    const ClusterList &clusterList1((TPC_VIEW_U == hitType) ? m_clusterListV : m_clusterListU);
    const ClusterList &clusterList2((TPC_VIEW_W == hitType) ? m_clusterListV : m_clusterListW);

    ClusterVector clusterVector1(clusterList1.begin(), clusterList1.end());
    ClusterVector clusterVector2(clusterList2.begin(), clusterList2.end());
    std::sort(clusterVector1.begin(), clusterVector1.end(), LArClusterHelper::SortByNHits);
    std::sort(clusterVector2.begin(), clusterVector2.end(), LArClusterHelper::SortByNHits);

    for (const Cluster *const pCluster1 : clusterVector1)
    {
        if (TPC_VIEW_U == hitType)
        {
            this->CalculateOverlapResult(pNewCluster, pCluster1, NULL);
        }
        else if (TPC_VIEW_V == hitType)
        {
            this->CalculateOverlapResult(pCluster1, pNewCluster, NULL);
        }
        else
        {
            this->CalculateOverlapResult(pCluster1, NULL, pNewCluster);
        }
    }

    for (const Cluster *const pCluster2 : clusterVector2)
    {
        if (TPC_VIEW_U == hitType)
        {
            this->CalculateOverlapResult(pNewCluster, NULL, pCluster2);
        }
        else if (TPC_VIEW_V == hitType)
        {
            this->CalculateOverlapResult(NULL, pNewCluster, pCluster2);
        }
        else
        {
            this->CalculateOverlapResult(NULL, pCluster2, pNewCluster);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewTrackFragmentsAlgorithm::RebuildClusters(const ClusterList &rebuildList, ClusterList &newClusters) const
{
    const ClusterList *pNewClusterList = NULL;
    std::string oldClusterListName, newClusterListName;

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeReclustering(*this, TrackList(), rebuildList,
        oldClusterListName));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunClusteringAlgorithm(*this, m_reclusteringAlgorithmName,
        pNewClusterList, newClusterListName));

    newClusters.insert(newClusters.end(), pNewClusterList->begin(), pNewClusterList->end());
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndReclustering(*this, newClusterListName));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewTrackFragmentsAlgorithm::PerformMainLoop()
{
    ClusterVector clusterVectorU(m_clusterListU.begin(), m_clusterListU.end());
    ClusterVector clusterVectorV(m_clusterListV.begin(), m_clusterListV.end());
    ClusterVector clusterVectorW(m_clusterListW.begin(), m_clusterListW.end());
    std::sort(clusterVectorU.begin(), clusterVectorU.end(), LArClusterHelper::SortByNHits);
    std::sort(clusterVectorV.begin(), clusterVectorV.end(), LArClusterHelper::SortByNHits);
    std::sort(clusterVectorW.begin(), clusterVectorW.end(), LArClusterHelper::SortByNHits);

    for (const Cluster *const pClusterU : clusterVectorU)
    {
        for (const Cluster *const pClusterV : clusterVectorV)
            this->CalculateOverlapResult(pClusterU, pClusterV, NULL);
    }

    for (const Cluster *const pClusterU : clusterVectorU)
    {
        for (const Cluster *const pClusterW : clusterVectorW)
            this->CalculateOverlapResult(pClusterU, NULL, pClusterW);
    }

    for (const Cluster *const pClusterV : clusterVectorV)
    {
        for (const Cluster *const pClusterW : clusterVectorW)
            this->CalculateOverlapResult(NULL, pClusterV, pClusterW);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewTrackFragmentsAlgorithm::CalculateOverlapResult(const Cluster *const pClusterU, const Cluster *const pClusterV, const Cluster *const pClusterW)
{
    const HitType missingHitType(
        ((NULL != pClusterU) && (NULL != pClusterV) && (NULL == pClusterW)) ? TPC_VIEW_W :
        ((NULL != pClusterU) && (NULL == pClusterV) && (NULL != pClusterW)) ? TPC_VIEW_V :
        ((NULL == pClusterU) && (NULL != pClusterV) && (NULL != pClusterW)) ? TPC_VIEW_U : HIT_CUSTOM);

    if (HIT_CUSTOM == missingHitType)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    // Calculate new overlap result and replace old overlap result where necessary
    FragmentOverlapResult oldOverlapResult, newOverlapResult;
    const Cluster *pMatchedClusterU(NULL), *pMatchedClusterV(NULL), *pMatchedClusterW(NULL);

    const TwoDSlidingFitResult &fitResult1((TPC_VIEW_U == missingHitType) ? this->GetCachedSlidingFitResult(pClusterV) :
        this->GetCachedSlidingFitResult(pClusterU));

    const TwoDSlidingFitResult &fitResult2((TPC_VIEW_U == missingHitType) ? this->GetCachedSlidingFitResult(pClusterW) :
        (TPC_VIEW_V == missingHitType) ? this->GetCachedSlidingFitResult(pClusterW) : this->GetCachedSlidingFitResult(pClusterV));

    const ClusterList &inputClusterList(this->GetInputClusterList(missingHitType));

    const Cluster *pBestMatchedCluster(NULL);
    const StatusCode statusCode(this->CalculateOverlapResult(fitResult1, fitResult2, inputClusterList, pBestMatchedCluster, newOverlapResult));

    if ((STATUS_CODE_SUCCESS != statusCode) && (STATUS_CODE_NOT_FOUND != statusCode))
        throw StatusCodeException(statusCode);

    if (!newOverlapResult.IsInitialized())
        return;

    if (STATUS_CODE_SUCCESS == statusCode)
    {
        pMatchedClusterU = ((NULL != pClusterU) ? pClusterU : pBestMatchedCluster);
        pMatchedClusterV = ((NULL != pClusterV) ? pClusterV : pBestMatchedCluster);
        pMatchedClusterW = ((NULL != pClusterW) ? pClusterW : pBestMatchedCluster);

        if (NULL == pMatchedClusterU || NULL == pMatchedClusterV || NULL == pMatchedClusterW)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        try
        {
            oldOverlapResult = m_overlapTensor.GetOverlapResult(pMatchedClusterU, pMatchedClusterV, pMatchedClusterW);
        }
        catch (StatusCodeException &)
        {
        }
    }

    if (!oldOverlapResult.IsInitialized())
    {
        m_overlapTensor.SetOverlapResult(pMatchedClusterU, pMatchedClusterV, pMatchedClusterW, newOverlapResult);
    }
    else if (newOverlapResult.GetFragmentCaloHitList().size() > oldOverlapResult.GetFragmentCaloHitList().size())
    {
        m_overlapTensor.ReplaceOverlapResult(pMatchedClusterU, pMatchedClusterV, pMatchedClusterW, newOverlapResult);
    }
    else if (newOverlapResult.GetFragmentCaloHitList().size() == oldOverlapResult.GetFragmentCaloHitList().size())
    {
        float newEnergySum(0.f), oldEnergySum(0.f);
        for (const CaloHit *const pCaloHit : newOverlapResult.GetFragmentCaloHitList()) newEnergySum += pCaloHit->GetHadronicEnergy();
        for (const CaloHit *const pCaloHit : oldOverlapResult.GetFragmentCaloHitList()) oldEnergySum += pCaloHit->GetHadronicEnergy();

        if (newEnergySum > oldEnergySum)
            m_overlapTensor.ReplaceOverlapResult(pMatchedClusterU, pMatchedClusterV, pMatchedClusterW, newOverlapResult);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeViewTrackFragmentsAlgorithm::CalculateOverlapResult(const TwoDSlidingFitResult &fitResult1, const TwoDSlidingFitResult &fitResult2,
    const ClusterList &inputClusterList, const Cluster *&pBestMatchedCluster, FragmentOverlapResult &fragmentOverlapResult) const
{
    const Cluster *const pCluster1(fitResult1.GetCluster());
    const Cluster *const pCluster2(fitResult2.GetCluster());

    float xMin1(0.f), xMax1(0.f), xMin2(0.f), xMax2(0.f);
    LArClusterHelper::GetClusterSpanX(pCluster1, xMin1, xMax1);
    LArClusterHelper::GetClusterSpanX(pCluster2, xMin2, xMax2);

    const float xOverlap(std::min(xMax1, xMax2) - std::max(xMin1, xMin2));

    if (xOverlap < std::numeric_limits<float>::epsilon())
        return STATUS_CODE_NOT_FOUND;

    CartesianPointVector projectedPositions;
    const StatusCode statusCode1(this->GetProjectedPositions(fitResult1, fitResult2, projectedPositions));

    if (STATUS_CODE_SUCCESS != statusCode1)
        return statusCode1;

    CaloHitList matchedHits;
    ClusterList matchedClusters;
    HitToClusterMap hitToClusterMap;
    const StatusCode statusCode2(this->GetMatchedHits(inputClusterList, projectedPositions, hitToClusterMap, matchedHits));

    if (STATUS_CODE_SUCCESS != statusCode2)
        return statusCode2;

    const StatusCode statusCode3(this->GetMatchedClusters(matchedHits, hitToClusterMap, matchedClusters, pBestMatchedCluster));

    if (STATUS_CODE_SUCCESS != statusCode3)
        return statusCode3;

    if (!this->CheckMatchedClusters(projectedPositions, matchedClusters))
        return STATUS_CODE_NOT_FOUND;

    FragmentOverlapResult overlapResult;
    this->GetFragmentOverlapResult(projectedPositions, matchedHits, matchedClusters, overlapResult);

    if (!this->CheckOverlapResult(overlapResult))
        return STATUS_CODE_NOT_FOUND;

    fragmentOverlapResult = overlapResult;
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeViewTrackFragmentsAlgorithm::GetProjectedPositions(const TwoDSlidingFitResult &fitResult1, const TwoDSlidingFitResult &fitResult2,
    CartesianPointVector &projectedPositions) const
{
    const Cluster *const pCluster1(fitResult1.GetCluster());
    const Cluster *const pCluster2(fitResult2.GetCluster());

    // Check hit types
    const HitType hitType1(LArClusterHelper::GetClusterHitType(pCluster1));
    const HitType hitType2(LArClusterHelper::GetClusterHitType(pCluster2));
    const HitType hitType3((TPC_VIEW_U != hitType1 && TPC_VIEW_U != hitType2) ? TPC_VIEW_U :
                           (TPC_VIEW_V != hitType1 && TPC_VIEW_V != hitType2) ? TPC_VIEW_V :
                           (TPC_VIEW_W != hitType1 && TPC_VIEW_W != hitType2) ? TPC_VIEW_W : HIT_CUSTOM);

    if (HIT_CUSTOM == hitType3)
        return STATUS_CODE_INVALID_PARAMETER;

    // Check absolute and fractional overlap in x coordinate
    float xMin1(0.f), xMax1(0.f), xMin2(0.f), xMax2(0.f);
    LArClusterHelper::GetClusterSpanX(pCluster1, xMin1, xMax1);
    LArClusterHelper::GetClusterSpanX(pCluster2, xMin2, xMax2);

    const float xOverlap(std::min(xMax1, xMax2) - std::max(xMin1, xMin2));
    const float xSpan(std::max(xMax1, xMax2) - std::min(xMin1, xMin2));

    if ((xOverlap < m_minXOverlap) || (xSpan < std::numeric_limits<float>::epsilon()) || ((xOverlap / xSpan) < m_minXOverlapFraction))
        return STATUS_CODE_NOT_FOUND;

    // Identify vertex and end positions (2D)
    const CartesianVector minPosition1(fitResult1.GetGlobalMinLayerPosition());
    const CartesianVector maxPosition1(fitResult1.GetGlobalMaxLayerPosition());
    const CartesianVector minPosition2(fitResult2.GetGlobalMinLayerPosition());
    const CartesianVector maxPosition2(fitResult2.GetGlobalMaxLayerPosition());

    const float dx_A(std::fabs(minPosition2.GetX() - minPosition1.GetX()));
    const float dx_B(std::fabs(maxPosition2.GetX() - maxPosition1.GetX()));
    const float dx_C(std::fabs(maxPosition2.GetX() - minPosition1.GetX()));
    const float dx_D(std::fabs(minPosition2.GetX() - maxPosition1.GetX()));

    if (std::min(dx_C,dx_D) > std::max(dx_A,dx_B) && std::min(dx_A,dx_B) > std::max(dx_C,dx_D))
       return STATUS_CODE_NOT_FOUND;

    const CartesianVector &vtxPosition1(minPosition1);
    const CartesianVector &endPosition1(maxPosition1);
    const CartesianVector &vtxPosition2((dx_A < dx_C) ? minPosition2 : maxPosition2);
    const CartesianVector &endPosition2((dx_A < dx_C) ? maxPosition2 : minPosition2);

    // Calculate vertex and end positions (3D)
    float vtxChi2(0.f);
    CartesianVector vtxPosition3D(0.f, 0.f, 0.f);
    LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), hitType1, hitType2, vtxPosition1, vtxPosition2, vtxPosition3D, vtxChi2);

    float endChi2(0.f);
    CartesianVector endPosition3D(0.f, 0.f, 0.f);
    LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), hitType1, hitType2, endPosition1, endPosition2, endPosition3D, endChi2);

    const CartesianVector vtxProjection3(LArGeometryHelper::ProjectPosition(this->GetPandora(), vtxPosition3D, hitType3));
    const CartesianVector endProjection3(LArGeometryHelper::ProjectPosition(this->GetPandora(), endPosition3D, hitType3));

    const float samplingPitch(0.5f * LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const float nSamplingPoints((endProjection3 - vtxProjection3).GetMagnitude() / samplingPitch);

    if (nSamplingPoints < 1.f)
        return STATUS_CODE_NOT_FOUND;

    // Loop over trajectory points
    bool foundLastPosition(false);
    CartesianVector lastPosition(0.f,0.f,0.f);

    for (float iSample = 0.5f; iSample < nSamplingPoints; iSample += 1.f)
    {
        const CartesianVector linearPosition3D(vtxPosition3D + (endPosition3D - vtxPosition3D) * (iSample / nSamplingPoints));
        const CartesianVector linearPosition1(LArGeometryHelper::ProjectPosition(this->GetPandora(), linearPosition3D, hitType1));
        const CartesianVector linearPosition2(LArGeometryHelper::ProjectPosition(this->GetPandora(), linearPosition3D, hitType2));

        float chi2(0.f);
        CartesianVector fitPosition1(0.f, 0.f, 0.f), fitPosition2(0.f, 0.f, 0.f);

        if ((STATUS_CODE_SUCCESS != fitResult1.GetGlobalFitProjection(linearPosition1, fitPosition1)) ||
            (STATUS_CODE_SUCCESS != fitResult2.GetGlobalFitProjection(linearPosition2, fitPosition2)))
        {
            continue;
        }

        float rL1(0.f), rL2(0.f), rT1(0.f), rT2(0.f);
        fitResult1.GetLocalPosition(fitPosition1, rL1, rT1);
        fitResult2.GetLocalPosition(fitPosition2, rL2, rT2);

        const float x(0.5 * (fitPosition1.GetX() + fitPosition2.GetX()));
        CartesianVector position1(0.f, 0.f, 0.f), position2(0.f, 0.f, 0.f), position3(0.f, 0.f, 0.f);

        try
        {
            const FitSegment &fitSegment1 = fitResult1.GetFitSegment(rL1);
            const FitSegment &fitSegment2 = fitResult2.GetFitSegment(rL2);

            if ((STATUS_CODE_SUCCESS != fitResult1.GetTransverseProjection(x, fitSegment1, position1)) ||
                (STATUS_CODE_SUCCESS != fitResult2.GetTransverseProjection(x, fitSegment2, position2)))
            {
                continue;
            }
        }
        catch (StatusCodeException &statusCodeException)
        {
            if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                throw statusCodeException;

            continue;
        }

        LArGeometryHelper::MergeTwoPositions(this->GetPandora(), hitType1, hitType2, position1, position2, position3, chi2);

        // TODO For highly multi-valued x, projected positions can be unreliable. Need to make interpolation more robust for these cases.
        if (foundLastPosition)
        {
            const float thisDisplacement((lastPosition - position3).GetMagnitude());
            if (thisDisplacement > 2.f * samplingPitch)
            {
                const float nExtraPoints(thisDisplacement / samplingPitch);
                for (float iExtra = 0.5f; iExtra < nExtraPoints; iExtra += 1.f)
                {
                    const CartesianVector extraPosition(position3 + (lastPosition - position3) * (iExtra / nExtraPoints));
                    projectedPositions.push_back(extraPosition);
                }
            }
        }

        projectedPositions.push_back(position3);

        lastPosition = position3;
        foundLastPosition = true;
    }

    // Bail out if list of projected positions is empty
    if (projectedPositions.empty())
        return STATUS_CODE_NOT_FOUND;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeViewTrackFragmentsAlgorithm::GetMatchedHits(const ClusterList &inputClusterList, const CartesianPointVector &projectedPositions,
    HitToClusterMap &hitToClusterMap, CaloHitList &matchedHits) const
{
    CaloHitVector availableCaloHits;

    for (ClusterList::const_iterator iter = inputClusterList.begin(), iterEnd = inputClusterList.end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pCluster = *iter;

        if (!pCluster->IsAvailable())
            continue;

        CaloHitList caloHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);
        availableCaloHits.insert(availableCaloHits.end(), caloHitList.begin(), caloHitList.end());

        for (CaloHitList::const_iterator hIter = caloHitList.begin(), hIterEnd = caloHitList.end(); hIter != hIterEnd; ++hIter)
            hitToClusterMap.insert(HitToClusterMap::value_type(*hIter, pCluster));
    }

    std::sort(availableCaloHits.begin(), availableCaloHits.end(), LArClusterHelper::SortHitsByPosition);

    for (const CartesianVector &projectedPosition : projectedPositions)
    {
        const CaloHit *pClosestCaloHit(NULL);
        float closestDistanceSquared(std::numeric_limits<float>::max()), tieBreakerBestEnergy(0.f);

        for (const CaloHit *const pCaloHit : availableCaloHits)
        {
            const float distanceSquared((pCaloHit->GetPositionVector() - projectedPosition).GetMagnitudeSquared());

            if ((distanceSquared < closestDistanceSquared) || ((std::fabs(distanceSquared - closestDistanceSquared) < std::numeric_limits<float>::epsilon()) && (pCaloHit->GetHadronicEnergy() > tieBreakerBestEnergy)))
            {
                pClosestCaloHit = pCaloHit;
                closestDistanceSquared = distanceSquared;
                tieBreakerBestEnergy = pCaloHit->GetHadronicEnergy();
            }
        }

        if ((closestDistanceSquared < m_maxPointDisplacementSquared) && (NULL != pClosestCaloHit) && (matchedHits.end() == std::find(matchedHits.begin(), matchedHits.end(), pClosestCaloHit)))
            matchedHits.push_back(pClosestCaloHit);
    }

    if (matchedHits.empty())
        return STATUS_CODE_NOT_FOUND;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeViewTrackFragmentsAlgorithm::GetMatchedClusters(const CaloHitList &matchedHits, const HitToClusterMap &hitToClusterMap,
    ClusterList &matchedClusters, const Cluster *&pBestMatchedCluster) const
{
    ClusterToMatchedHitsMap clusterToMatchedHitsMap;

    for (CaloHitList::const_iterator iter = matchedHits.begin(), iterEnd = matchedHits.end(); iter != iterEnd; ++iter)
    {
        const CaloHit *const pCaloHit = *iter;
        HitToClusterMap::const_iterator cIter = hitToClusterMap.find(pCaloHit);

        if (hitToClusterMap.end() == cIter)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        ++clusterToMatchedHitsMap[cIter->second];

        if (matchedClusters.end() == std::find(matchedClusters.begin(), matchedClusters.end(), cIter->second))
            matchedClusters.push_back(cIter->second);
    }

    if (matchedClusters.empty())
        return STATUS_CODE_NOT_FOUND;

    pBestMatchedCluster = NULL;
    unsigned int bestClusterMatchedHits(0);
    float tieBreakerBestEnergy(0.f);

    ClusterVector sortedClusters;
    for (const auto &mapEntry : clusterToMatchedHitsMap) sortedClusters.push_back(mapEntry.first);
    std::sort(sortedClusters.begin(), sortedClusters.end(), LArClusterHelper::SortByNHits);

    for (const Cluster *const pCluster : sortedClusters)
    {
        const unsigned int nMatchedHits(clusterToMatchedHitsMap.at(pCluster));

        if ((nMatchedHits > bestClusterMatchedHits) || ((nMatchedHits == bestClusterMatchedHits) && (pCluster->GetHadronicEnergy() > tieBreakerBestEnergy)))
        {
            pBestMatchedCluster = pCluster;
            bestClusterMatchedHits = nMatchedHits;
            tieBreakerBestEnergy = pCluster->GetHadronicEnergy();
        }
    }

    if (NULL == pBestMatchedCluster)
        return STATUS_CODE_NOT_FOUND;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewTrackFragmentsAlgorithm::GetFragmentOverlapResult(const CartesianPointVector &projectedPositions, const CaloHitList &matchedHits,
    const ClusterList &matchedClusters, FragmentOverlapResult &fragmentOverlapResult) const
{
    float chi2Sum(0.f);
    unsigned int nMatchedSamplingPoints(0);

    CaloHitVector sortedMatchedHits(matchedHits.begin(), matchedHits.end());
    std::sort(sortedMatchedHits.begin(), sortedMatchedHits.end(), LArClusterHelper::SortHitsByPosition);

    for (const CartesianVector &projectedPosition : projectedPositions)
    {
        float closestDistanceSquared(std::numeric_limits<float>::max());

        for (const CaloHit *const pCaloHit : matchedHits)
        {
            const float distanceSquared((pCaloHit->GetPositionVector() - projectedPosition).GetMagnitudeSquared());

            if (distanceSquared < closestDistanceSquared)
                closestDistanceSquared = distanceSquared;
        }

        if (closestDistanceSquared < m_maxPointDisplacementSquared)
        {
            ++nMatchedSamplingPoints;
            chi2Sum += closestDistanceSquared;
        }
    }

    fragmentOverlapResult = FragmentOverlapResult(nMatchedSamplingPoints, projectedPositions.size(), chi2Sum, matchedHits, matchedClusters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ThreeViewTrackFragmentsAlgorithm::CheckMatchedClusters(const CartesianPointVector &projectedPositions, const ClusterList &matchedClusters) const
{
    if (projectedPositions.empty() || matchedClusters.empty())
        return false;

    // Calculate X and Z span of projected positions
    float minXproj(+std::numeric_limits<float>::max());
    float maxXproj(-std::numeric_limits<float>::max());
    float minZproj(+std::numeric_limits<float>::max());
    float maxZproj(-std::numeric_limits<float>::max());

    for (const CartesianVector &projectedPosition : projectedPositions)
    {
        minXproj = std::min(minXproj, projectedPosition.GetX());
        maxXproj = std::max(maxXproj, projectedPosition.GetX());
        minZproj = std::min(minZproj, projectedPosition.GetZ());
        maxZproj = std::max(maxZproj, projectedPosition.GetZ());
    }

    const float dXproj(maxXproj - minXproj);
    const float dZproj(maxZproj - minZproj);
    const float projectedLengthSquared(dXproj * dXproj + dZproj * dZproj);

    // Calculate X and Z span of matched clusters
    float minXcluster(+std::numeric_limits<float>::max());
    float maxXcluster(-std::numeric_limits<float>::max());
    float minZcluster(+std::numeric_limits<float>::max());
    float maxZcluster(-std::numeric_limits<float>::max());

    for (const Cluster *const pCluster : matchedClusters)
    {
        CartesianVector minPosition(0.f,0.f,0.f);
        CartesianVector maxPosition(0.f,0.f,0.f);

        LArClusterHelper::GetClusterBoundingBox(pCluster, minPosition, maxPosition);

        minXcluster = std::min(minXcluster, minPosition.GetX());
        maxXcluster = std::max(maxXcluster, maxPosition.GetX());
        minZcluster = std::min(minZcluster, minPosition.GetZ());
        maxZcluster = std::max(maxZcluster, maxPosition.GetZ());
    }

    const float dXcluster(maxXcluster - minXcluster);
    const float dZcluster(maxZcluster - minZcluster);
    const float clusterLengthSquared(dXcluster * dXcluster + dZcluster * dZcluster);

    // Require that the span of the matched clusters is no larger than twice the span of the projected positions
    if (clusterLengthSquared > 4.f * projectedLengthSquared)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ThreeViewTrackFragmentsAlgorithm::CheckOverlapResult(const FragmentOverlapResult &overlapResult) const
{
    // ATTN This method is currently mirrored in ClearTrackFragments tool
    if (overlapResult.GetMatchedFraction() < m_minMatchedSamplingPointFraction)
        return false;

    if (overlapResult.GetFragmentCaloHitList().size() < m_minMatchedHits)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewTrackFragmentsAlgorithm::ExamineOverlapContainer()
{
    unsigned int repeatCounter(0);

    for (TensorToolVector::const_iterator iter = m_algorithmToolVector.begin(), iterEnd = m_algorithmToolVector.end(); iter != iterEnd; )
    {
        if ((*iter)->Run(this, m_overlapTensor))
        {
            iter = m_algorithmToolVector.begin();

            if (++repeatCounter > m_nMaxTensorToolRepeats)
                break;
        }
        else
        {
            ++iter;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeViewTrackFragmentsAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, xmlHandle,
        "ClusterRebuilding", m_reclusteringAlgorithmName));

    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle,
        "TrackTools", algorithmToolVector));

    for (AlgorithmToolVector::const_iterator iter = algorithmToolVector.begin(), iterEnd = algorithmToolVector.end(); iter != iterEnd; ++iter)
    {
        FragmentTensorTool *const pFragmentTensorTool(dynamic_cast<FragmentTensorTool*>(*iter));

        if (!pFragmentTensorTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolVector.push_back(pFragmentTensorTool);
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NMaxTensorToolRepeats", m_nMaxTensorToolRepeats));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinXOverlap", m_minXOverlap));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinXOverlapFraction", m_minXOverlapFraction));

    float maxPointDisplacement = std::sqrt(m_maxPointDisplacementSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxPointDisplacement", maxPointDisplacement));
    m_maxPointDisplacementSquared = maxPointDisplacement * maxPointDisplacement;

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedSamplingPointFraction", m_minMatchedSamplingPointFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedHits", m_minMatchedHits));

    return ThreeViewTrackMatchingAlgorithm<FragmentOverlapResult>::ReadSettings(xmlHandle);
}

} // namespace lar_content
