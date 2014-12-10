/**
 *  @file   LArContent/src/LArThreeDReco/LArTrackFragments/ThreeDTrackFragmentsAlgorithm.cc
 *
 *  @brief  Implementation of the three dimensional fragments algorithm base class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"

#include "LArPlugins/LArTransformationPlugin.h"

#include "LArThreeDReco/LArTrackFragments/ThreeDTrackFragmentsAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ThreeDTrackFragmentsAlgorithm::ThreeDTrackFragmentsAlgorithm() :
    m_nMaxTensorToolRepeats(5000),
    m_minXOverlap(3.f),
    m_minXOverlapFraction(0.8f),
    m_maxPointDisplacementSquared(1.5f * 1.5f),
    m_minMatchedSamplingPointFraction(0.5f),
    m_minMatchedHits(5)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDTrackFragmentsAlgorithm::UpdateForNewCluster(Cluster *const pNewCluster)
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

    if (!clusterList.insert(pNewCluster).second)
        throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);

    const ClusterList &clusterList1((TPC_VIEW_U == hitType) ? m_clusterListV : m_clusterListU);
    const ClusterList &clusterList2((TPC_VIEW_W == hitType) ? m_clusterListV : m_clusterListW);

    for (ClusterList::const_iterator iter1 = clusterList1.begin(), iter1End = clusterList1.end(); iter1 != iter1End; ++iter1)
    {
        if (TPC_VIEW_U == hitType)
        {
            this->CalculateOverlapResult(pNewCluster, *iter1, NULL);
        }
        else if (TPC_VIEW_V == hitType)
        {
            this->CalculateOverlapResult(*iter1, pNewCluster, NULL);
        }
        else
        {
            this->CalculateOverlapResult(*iter1, NULL, pNewCluster);
        }
    }

    for (ClusterList::const_iterator iter2 = clusterList2.begin(), iter2End = clusterList2.end(); iter2 != iter2End; ++iter2)
    {
        if (TPC_VIEW_U == hitType)
        {
            this->CalculateOverlapResult(pNewCluster, NULL, *iter2);
        }
        else if (TPC_VIEW_V == hitType)
        {
            this->CalculateOverlapResult(NULL, pNewCluster, *iter2);
        }
        else
        {
            this->CalculateOverlapResult(NULL, *iter2, pNewCluster);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDTrackFragmentsAlgorithm::RebuildClusters(const ClusterList &rebuildList, ClusterList &newClusters) const
{
    const ClusterList *pNewClusterList = NULL;
    std::string oldClusterListName, newClusterListName;

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeReclustering(*this, TrackList(), rebuildList,
        oldClusterListName));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunClusteringAlgorithm(*this, m_reclusteringAlgorithmName,
        pNewClusterList, newClusterListName));

    newClusters.insert(pNewClusterList->begin(), pNewClusterList->end());
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndReclustering(*this, newClusterListName));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDTrackFragmentsAlgorithm::PerformMainLoop()
{
    for (ClusterList::const_iterator iterU = m_clusterListU.begin(), iterUEnd = m_clusterListU.end(); iterU != iterUEnd; ++iterU)
    {
        for (ClusterList::const_iterator iterV = m_clusterListV.begin(), iterVEnd = m_clusterListV.end(); iterV != iterVEnd; ++iterV)
            this->CalculateOverlapResult(*iterU, *iterV, NULL);
    }

    for (ClusterList::const_iterator iterU = m_clusterListU.begin(), iterUEnd = m_clusterListU.end(); iterU != iterUEnd; ++iterU)
    {
        for (ClusterList::const_iterator iterW = m_clusterListW.begin(), iterWEnd = m_clusterListW.end(); iterW != iterWEnd; ++iterW)
            this->CalculateOverlapResult(*iterU, NULL, *iterW);
    }

    for (ClusterList::const_iterator iterV = m_clusterListV.begin(), iterVEnd = m_clusterListV.end(); iterV != iterVEnd; ++iterV)
    {
        for (ClusterList::const_iterator iterW = m_clusterListW.begin(), iterWEnd = m_clusterListW.end(); iterW != iterWEnd; ++iterW)
            this->CalculateOverlapResult(NULL, *iterV, *iterW);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDTrackFragmentsAlgorithm::CalculateOverlapResult(Cluster *pClusterU, Cluster *pClusterV, Cluster *pClusterW)
{
    const HitType missingHitType(
        ((NULL != pClusterU) && (NULL != pClusterV) && (NULL == pClusterW)) ? TPC_VIEW_W :
        ((NULL != pClusterU) && (NULL == pClusterV) && (NULL != pClusterW)) ? TPC_VIEW_V :
        ((NULL == pClusterU) && (NULL != pClusterV) && (NULL != pClusterW)) ? TPC_VIEW_U : HIT_CUSTOM);

    if (HIT_CUSTOM == missingHitType)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    // Calculate new overlap result and replace old overlap result where necessary
    FragmentOverlapResult oldOverlapResult, newOverlapResult;
    Cluster *pMatchedClusterU(NULL), *pMatchedClusterV(NULL), *pMatchedClusterW(NULL);

    try
    {
        const TwoDSlidingFitResult &fitResult1((TPC_VIEW_U == missingHitType) ? this->GetCachedSlidingFitResult(pClusterV) :
            (TPC_VIEW_V == missingHitType) ? this->GetCachedSlidingFitResult(pClusterU) : this->GetCachedSlidingFitResult(pClusterU));

        const TwoDSlidingFitResult &fitResult2((TPC_VIEW_U == missingHitType) ? this->GetCachedSlidingFitResult(pClusterW) :
            (TPC_VIEW_V == missingHitType) ? this->GetCachedSlidingFitResult(pClusterW) : this->GetCachedSlidingFitResult(pClusterV));

        const ClusterList &inputClusterList((TPC_VIEW_U == missingHitType) ? this->GetInputClusterListU() :
            (TPC_VIEW_V == missingHitType) ? this->GetInputClusterListV() : this->GetInputClusterListW());

        Cluster *pBestMatchedCluster(NULL);
        this->CalculateOverlapResult(fitResult1, fitResult2, inputClusterList, pBestMatchedCluster, newOverlapResult);

        pMatchedClusterU = ((NULL != pClusterU) ? pClusterU : pBestMatchedCluster);
        pMatchedClusterV = ((NULL != pClusterV) ? pClusterV : pBestMatchedCluster);
        pMatchedClusterW = ((NULL != pClusterW) ? pClusterW : pBestMatchedCluster);

        if (NULL == pMatchedClusterU || NULL == pMatchedClusterV || NULL == pMatchedClusterW)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        oldOverlapResult = m_overlapTensor.GetOverlapResult(pMatchedClusterU, pMatchedClusterV, pMatchedClusterW);
    }
    catch (StatusCodeException &statusCodeException)
    {
        if (!(STATUS_CODE_NOT_FOUND == statusCodeException.GetStatusCode() || STATUS_CODE_NOT_INITIALIZED == statusCodeException.GetStatusCode()))
            throw statusCodeException;
    }

    if (!newOverlapResult.IsInitialized())
        return;

    if (!oldOverlapResult.IsInitialized())
    {
        m_overlapTensor.SetOverlapResult(pMatchedClusterU, pMatchedClusterV, pMatchedClusterW, newOverlapResult);
    }
    else if(newOverlapResult.GetFragmentCaloHitList().size() > oldOverlapResult.GetFragmentCaloHitList().size())
    {
        m_overlapTensor.ReplaceOverlapResult(pMatchedClusterU, pMatchedClusterV, pMatchedClusterW, newOverlapResult);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDTrackFragmentsAlgorithm::CalculateOverlapResult(const TwoDSlidingFitResult &fitResult1, const TwoDSlidingFitResult &fitResult2,
    const ClusterList &inputClusterList, Cluster *&pBestMatchedCluster, FragmentOverlapResult &fragmentOverlapResult) const
{
    const Cluster *pCluster1(fitResult1.GetCluster());
    const Cluster *pCluster2(fitResult2.GetCluster());

    float xMin1(0.f), xMax1(0.f), xMin2(0.f), xMax2(0.f);
    LArClusterHelper::GetClusterSpanX(pCluster1, xMin1, xMax1);
    LArClusterHelper::GetClusterSpanX(pCluster2, xMin2, xMax2);

    const float xOverlap(std::min(xMax1, xMax2) - std::max(xMin1, xMin2));

    if (xOverlap < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    CartesianPointList projectedPositions;
    this->GetProjectedPositions(fitResult1, fitResult2, projectedPositions);

    CaloHitList matchedHits;
    ClusterList matchedClusters;
    HitToClusterMap hitToClusterMap;
    this->GetMatchedHits(inputClusterList, projectedPositions, hitToClusterMap, matchedHits);
    this->GetMatchedClusters(matchedHits, hitToClusterMap, matchedClusters, pBestMatchedCluster);

    if (!this->CheckMatchedClusters(projectedPositions, matchedClusters))
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    FragmentOverlapResult overlapResult;
    this->GetFragmentOverlapResult(projectedPositions, matchedHits, matchedClusters, overlapResult);

    if (!this->CheckOverlapResult(overlapResult))
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    fragmentOverlapResult = overlapResult;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDTrackFragmentsAlgorithm::GetProjectedPositions(const TwoDSlidingFitResult &fitResult1, const TwoDSlidingFitResult &fitResult2,
    CartesianPointList &projectedPositions) const
{
    const Cluster *pCluster1(fitResult1.GetCluster());
    const Cluster *pCluster2(fitResult2.GetCluster());

    // Check hit types
    const HitType hitType1(LArClusterHelper::GetClusterHitType(pCluster1));
    const HitType hitType2(LArClusterHelper::GetClusterHitType(pCluster2));
    const HitType hitType3((TPC_VIEW_U != hitType1 && TPC_VIEW_U != hitType2) ? TPC_VIEW_U :
                           (TPC_VIEW_V != hitType1 && TPC_VIEW_V != hitType2) ? TPC_VIEW_V :
                           (TPC_VIEW_W != hitType1 && TPC_VIEW_W != hitType2) ? TPC_VIEW_W : HIT_CUSTOM);

    if (HIT_CUSTOM == hitType3)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    // Check absolute and fractional overlap in x coordinate
    float xMin1(0.f), xMax1(0.f), xMin2(0.f), xMax2(0.f);
    LArClusterHelper::GetClusterSpanX(pCluster1, xMin1, xMax1);
    LArClusterHelper::GetClusterSpanX(pCluster2, xMin2, xMax2);

    const float xOverlap(std::min(xMax1, xMax2) - std::max(xMin1, xMin2));
    const float xSpan(std::max(xMax1, xMax2) - std::min(xMin1, xMin2));

    if ((xOverlap < m_minXOverlap) || (xSpan < std::numeric_limits<float>::epsilon()) || ((xOverlap / xSpan) < m_minXOverlapFraction))
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

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
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

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

    const float samplingPitch(0.5f * LArGeometryHelper::GetLArTransformationPlugin(this->GetPandora())->GetWireZPitch());
    const float nSamplingPoints((endProjection3 - vtxProjection3).GetMagnitude() / samplingPitch);

    if (nSamplingPoints < 1.f)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    // Loop over trajectory points
    bool foundLastPosition(false);
    CartesianVector lastPosition(0.f,0.f,0.f);

    for (float iSample = 0.5f; iSample < nSamplingPoints; iSample += 1.f)
    {
        const CartesianVector linearPosition3D(vtxPosition3D + (endPosition3D - vtxPosition3D) * (iSample / nSamplingPoints));
        const CartesianVector linearPosition1(LArGeometryHelper::ProjectPosition(this->GetPandora(), linearPosition3D, hitType1));
        const CartesianVector linearPosition2(LArGeometryHelper::ProjectPosition(this->GetPandora(), linearPosition3D, hitType2));

        try
        {
            float chi2(0.f);
            CartesianVector fitPosition1(0.f, 0.f, 0.f), fitPosition2(0.f, 0.f, 0.f);
            fitResult1.GetGlobalFitProjection(linearPosition1, fitPosition1);
            fitResult2.GetGlobalFitProjection(linearPosition2, fitPosition2);

            float rL1(0.f), rL2(0.f), rT1(0.f), rT2(0.f);
            fitResult1.GetLocalPosition(fitPosition1, rL1, rT1);
            fitResult2.GetLocalPosition(fitPosition2, rL2, rT2);

            const FitSegment &fitSegment1 = fitResult1.GetFitSegment(rL1);
            const FitSegment &fitSegment2 = fitResult2.GetFitSegment(rL2);

            const float x(0.5 * (fitPosition1.GetX() + fitPosition2.GetX()));
            CartesianVector position1(0.f, 0.f, 0.f), position2(0.f, 0.f, 0.f), position3(0.f, 0.f, 0.f);
            fitResult1.GetTransverseProjection(x, fitSegment1, position1);
            fitResult2.GetTransverseProjection(x, fitSegment2, position2);
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
        catch (StatusCodeException &)
        {
        }
    }

    // Bail out if list of projected positions is empty
    if (projectedPositions.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDTrackFragmentsAlgorithm::GetMatchedHits(const ClusterList &inputClusterList, const CartesianPointList &projectedPositions,
    HitToClusterMap &hitToClusterMap, CaloHitList &matchedHits) const
{
    CaloHitList availableCaloHits;

    for(ClusterList::const_iterator iter = inputClusterList.begin(), iterEnd = inputClusterList.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (!pCluster->IsAvailable())
            continue;

        CaloHitList caloHitList;
        pCluster->GetOrderedCaloHitList().GetCaloHitList(caloHitList);
        availableCaloHits.insert(caloHitList.begin(), caloHitList.end());

        for (CaloHitList::const_iterator hIter = caloHitList.begin(), hIterEnd = caloHitList.end(); hIter != hIterEnd; ++hIter)
            hitToClusterMap.insert(HitToClusterMap::value_type(*hIter, pCluster));
    }

    for (CartesianPointList::const_iterator iter = projectedPositions.begin(), iterEnd = projectedPositions.end(); iter != iterEnd; ++iter)
    {
        const CartesianVector &projectedPosition = *iter;
        float closestDistanceSquared(std::numeric_limits<float>::max());
        CaloHit *pClosestCaloHit(NULL);

        for (CaloHitList::const_iterator hIter = availableCaloHits.begin(), hIterEnd = availableCaloHits.end(); hIter != hIterEnd; ++hIter)
        {
            CaloHit *pCaloHit = *hIter;
            const float distanceSquared((pCaloHit->GetPositionVector() - projectedPosition).GetMagnitudeSquared());

            if (distanceSquared < closestDistanceSquared)
            {
                closestDistanceSquared = distanceSquared;
                pClosestCaloHit = pCaloHit;
            }
        }

        if ((closestDistanceSquared < m_maxPointDisplacementSquared) && (NULL != pClosestCaloHit))
            matchedHits.insert(pClosestCaloHit);
    }

    if (matchedHits.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDTrackFragmentsAlgorithm::GetMatchedClusters(const CaloHitList &matchedHits, const HitToClusterMap &hitToClusterMap,
    ClusterList &matchedClusters, Cluster *&pBestMatchedCluster) const
{
    ClusterToMatchedHitsMap clusterToMatchedHitsMap;

    for (CaloHitList::const_iterator iter = matchedHits.begin(), iterEnd = matchedHits.end(); iter != iterEnd; ++iter)
    {
        CaloHit *pCaloHit = *iter;
        HitToClusterMap::const_iterator cIter = hitToClusterMap.find(pCaloHit);

        if (hitToClusterMap.end() == cIter)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        matchedClusters.insert(cIter->second);
        ++clusterToMatchedHitsMap[cIter->second];
    }

    if (matchedClusters.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    pBestMatchedCluster = NULL;
    unsigned int bestClusterMatchedHits(0);

    for (ClusterToMatchedHitsMap::const_iterator iter = clusterToMatchedHitsMap.begin(), iterEnd = clusterToMatchedHitsMap.end(); iter != iterEnd; ++iter)
    {
        if (iter->second > bestClusterMatchedHits)
        {
            pBestMatchedCluster = iter->first;
            bestClusterMatchedHits = iter->second;
        }
    }

    if (NULL == pBestMatchedCluster)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDTrackFragmentsAlgorithm::GetFragmentOverlapResult(const CartesianPointList &projectedPositions, const CaloHitList &matchedHits,
    const ClusterList &matchedClusters, FragmentOverlapResult &fragmentOverlapResult) const
{
    float chi2Sum(0.f);
    unsigned int nMatchedSamplingPoints(0);

    for (CartesianPointList::const_iterator iter = projectedPositions.begin(), iterEnd = projectedPositions.end(); iter != iterEnd; ++iter)
    {
        const CartesianVector &projectedPosition = *iter;
        float closestDistanceSquared(std::numeric_limits<float>::max());

        for (CaloHitList::const_iterator hIter = matchedHits.begin(), hIterEnd = matchedHits.end(); hIter != hIterEnd; ++hIter)
        {
            CaloHit *pCaloHit = *hIter;
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

bool ThreeDTrackFragmentsAlgorithm::CheckMatchedClusters(const CartesianPointList &projectedPositions, const ClusterList &matchedClusters) const
{
    if (projectedPositions.empty() || matchedClusters.empty())
        return false;

    // Calculate X and Z span of projected positions
    float minXproj(+std::numeric_limits<float>::max());
    float maxXproj(-std::numeric_limits<float>::max());
    float minZproj(+std::numeric_limits<float>::max());
    float maxZproj(-std::numeric_limits<float>::max());

    for (CartesianPointList::const_iterator iter = projectedPositions.begin(), iterEnd = projectedPositions.end(); iter != iterEnd; ++iter)
    {
        const CartesianVector &projectedPosition = *iter;

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

    for(ClusterList::const_iterator iter = matchedClusters.begin(), iterEnd = matchedClusters.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

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

bool ThreeDTrackFragmentsAlgorithm::CheckOverlapResult(const FragmentOverlapResult &overlapResult) const
{
    // ATTN This method is currently mirrored in ClearTrackFragments tool
    if (overlapResult.GetMatchedFraction() < m_minMatchedSamplingPointFraction)
        return false;

    if (overlapResult.GetFragmentCaloHitList().size() < m_minMatchedHits)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDTrackFragmentsAlgorithm::ExamineTensor()
{
    unsigned int repeatCounter(0);

    for (TensorToolList::const_iterator iter = m_algorithmToolList.begin(), iterEnd = m_algorithmToolList.end(); iter != iterEnd; )
    {
        if ((*iter)->Run(this, m_overlapTensor))
        {
            iter = m_algorithmToolList.begin();

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

StatusCode ThreeDTrackFragmentsAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, xmlHandle,
        "ClusterRebuilding", m_reclusteringAlgorithmName));

    AlgorithmToolList algorithmToolList;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle,
        "TrackTools", algorithmToolList));

    for (AlgorithmToolList::const_iterator iter = algorithmToolList.begin(), iterEnd = algorithmToolList.end(); iter != iterEnd; ++iter)
    {
        FragmentTensorTool *pFragmentTensorTool(dynamic_cast<FragmentTensorTool*>(*iter));

        if (NULL == pFragmentTensorTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolList.push_back(pFragmentTensorTool);
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

    return ThreeDTracksBaseAlgorithm<FragmentOverlapResult>::ReadSettings(xmlHandle);
}

} // namespace lar_content
