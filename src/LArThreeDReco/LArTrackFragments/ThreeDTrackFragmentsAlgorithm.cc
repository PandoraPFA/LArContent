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
#include "LArHelpers/LArThreeDHelper.h"

#include "LArThreeDReco/LArTrackFragments/ThreeDTrackFragmentsAlgorithm.h"

using namespace pandora;

namespace lar
{

void ThreeDTrackFragmentsAlgorithm::UpdateForNewCluster(Cluster *const pNewCluster)
{
    this->AddToSlidingFitCache(pNewCluster);

    const HitType hitType(LArThreeDHelper::GetClusterHitType(pNewCluster));

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

void ThreeDTrackFragmentsAlgorithm::RebuildClusters(const CaloHitList &caloHitList, ClusterList &newClusters) const
{
    const ClusterList *pTemporaryList = NULL;
    std::string currentListName, temporaryListName;

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*this, currentListName));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent<ClusterList>(*this,
        pTemporaryList, temporaryListName));

    this->BuildNewClusters(caloHitList, newClusters);

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, temporaryListName, currentListName))
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, currentListName));
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
        ((NULL == pClusterU) && (NULL != pClusterV) && (NULL != pClusterW)) ? TPC_VIEW_U : CUSTOM);

    if (CUSTOM == missingHitType)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const TwoDSlidingFitResult &fitResult1((TPC_VIEW_U == missingHitType) ? this->GetCachedSlidingFitResult(pClusterV) :
        (TPC_VIEW_V == missingHitType) ? this->GetCachedSlidingFitResult(pClusterU) : this->GetCachedSlidingFitResult(pClusterU));

    const TwoDSlidingFitResult &fitResult2((TPC_VIEW_U == missingHitType) ? this->GetCachedSlidingFitResult(pClusterW) :
        (TPC_VIEW_V == missingHitType) ? this->GetCachedSlidingFitResult(pClusterW) : this->GetCachedSlidingFitResult(pClusterV));

    const ClusterList &inputClusterList((TPC_VIEW_U == missingHitType) ? this->GetInputClusterListU() :
        (TPC_VIEW_V == missingHitType) ? this->GetInputClusterListV() : this->GetInputClusterListW());


    // Calculate new overlap result and replace old overlap result where necessary
    FragmentOverlapResult oldOverlapResult, newOverlapResult;
    Cluster *pMatchedClusterU(NULL), *pMatchedClusterV(NULL), *pMatchedClusterW(NULL);

    try
    {
        Cluster *pBestMatchedCluster(NULL);
        this->CalculateOverlapResult(fitResult1, fitResult2, inputClusterList, pBestMatchedCluster, newOverlapResult);

        pMatchedClusterU = ((NULL != pClusterU) ? pClusterU : pBestMatchedCluster);
        pMatchedClusterV = ((NULL != pClusterV) ? pClusterV : pBestMatchedCluster);
        pMatchedClusterW = ((NULL != pClusterW) ? pClusterW : pBestMatchedCluster);

        oldOverlapResult = m_overlapTensor.GetOverlapResult(pMatchedClusterU, pMatchedClusterV, pMatchedClusterW);
    }
    catch (StatusCodeException &statusCodeException)
    {
        if (!(STATUS_CODE_NOT_FOUND == statusCodeException.GetStatusCode() || STATUS_CODE_NOT_INITIALIZED == statusCodeException.GetStatusCode()))
            throw statusCodeException;
    }

    if (!newOverlapResult.IsInitialized())
        return;

    if (NULL == pMatchedClusterU || NULL == pMatchedClusterV || NULL == pMatchedClusterW)
        throw StatusCodeException(STATUS_CODE_FAILURE);

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

    CaloHitList associatedHits;
    HitToClusterMap hitToClusterMap;
    this->GetAssociatedHits(inputClusterList, projectedPositions, hitToClusterMap, associatedHits);

    CaloHitList matchedHits;
    this->GetMatchedHits(associatedHits, matchedHits);

    ClusterList matchedClusters;
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
    const HitType hitType1(LArThreeDHelper::GetClusterHitType(pCluster1));
    const HitType hitType2(LArThreeDHelper::GetClusterHitType(pCluster2));
    const HitType hitType3((TPC_VIEW_U != hitType1 && TPC_VIEW_U != hitType2) ? TPC_VIEW_U :
                           (TPC_VIEW_V != hitType1 && TPC_VIEW_V != hitType2) ? TPC_VIEW_V :
                           (TPC_VIEW_W != hitType1 && TPC_VIEW_W != hitType2) ? TPC_VIEW_W : CUSTOM);

    if (CUSTOM == hitType3)
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
    LArGeometryHelper::MergeTwoPositions3D(hitType1, hitType2, vtxPosition1, vtxPosition2, vtxPosition3D, vtxChi2);

    float endChi2(0.f);
    CartesianVector endPosition3D(0.f, 0.f, 0.f);
    LArGeometryHelper::MergeTwoPositions3D(hitType1, hitType2, endPosition1, endPosition2, endPosition3D, endChi2);

    const CartesianVector vtxProjection1(LArGeometryHelper::ProjectPosition(vtxPosition3D, hitType1));
    const CartesianVector vtxProjection2(LArGeometryHelper::ProjectPosition(vtxPosition3D, hitType2));
    const CartesianVector vtxProjection3(LArGeometryHelper::ProjectPosition(vtxPosition3D, hitType3));

    const CartesianVector endProjection1(LArGeometryHelper::ProjectPosition(endPosition3D, hitType1));
    const CartesianVector endProjection2(LArGeometryHelper::ProjectPosition(endPosition3D, hitType2));
    const CartesianVector endProjection3(LArGeometryHelper::ProjectPosition(endPosition3D, hitType3));

    const float nSamplingPoints(3.f * (endProjection3 - vtxProjection3).GetMagnitude() / m_maxPointDisplacement);

    if (nSamplingPoints < 1.f)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    // Loop over trajectory points
    for (float iSample = 0.f; iSample < nSamplingPoints; iSample += 1.f)
    {
        const CartesianVector linearPosition3D(vtxPosition3D + (endPosition3D - vtxPosition3D) * ((0.5f + iSample) / nSamplingPoints));
        const CartesianVector linearPosition1(LArGeometryHelper::ProjectPosition(linearPosition3D, hitType1));
        const CartesianVector linearPosition2(LArGeometryHelper::ProjectPosition(linearPosition3D, hitType2));

        try
        {
            float chi2(0.f);
            CartesianVector fitPosition1(0.f, 0.f, 0.f), fitPosition2(0.f, 0.f, 0.f);
            fitResult1.GetGlobalFitProjection(linearPosition1, fitPosition1);
            fitResult2.GetGlobalFitProjection(linearPosition2, fitPosition2);

            float rL1(0.f), rL2(0.f), rT1(0.f), rT2(0.f);
            fitResult1.GetLocalPosition(fitPosition1, rL1, rT1);
            fitResult2.GetLocalPosition(fitPosition2, rL2, rT2);

            const TwoDSlidingFitResult::FitSegment& fitSegment1 = fitResult1.GetFitSegment(rL1);
            const TwoDSlidingFitResult::FitSegment& fitSegment2 = fitResult2.GetFitSegment(rL2);

            const float x(0.5 * (fitPosition1.GetX() + fitPosition2.GetX()));
            CartesianVector position1(0.f, 0.f, 0.f), position2(0.f, 0.f, 0.f), position3(0.f, 0.f, 0.f);
            fitResult1.GetTransverseProjection(x, fitSegment1, position1);
            fitResult2.GetTransverseProjection(x, fitSegment2, position2);
            LArGeometryHelper::MergeTwoPositions(hitType1, hitType2, position1, position2, position3, chi2);
            projectedPositions.push_back(position3);
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

void ThreeDTrackFragmentsAlgorithm::GetAssociatedHits(const ClusterList &inputClusterList, const CartesianPointList &projectedPositions,
    HitToClusterMap &hitToClusterMap, CaloHitList &associatedHits) const
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
            associatedHits.insert(pClosestCaloHit);
    }

    if (associatedHits.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDTrackFragmentsAlgorithm::GetMatchedHits(const CaloHitList &associatedHits, CaloHitList &matchedHits) const
{
    for (CaloHitList::const_iterator hIter1 = associatedHits.begin(), hIterEnd1 = associatedHits.end(); hIter1 != hIterEnd1; ++hIter1)
    {
        CaloHit *pCaloHit1 = *hIter1;
        bool isGoodHit(false);

        for (CaloHitList::const_iterator hIter2 = associatedHits.begin(), hIterEnd2 = associatedHits.end(); hIter2 != hIterEnd2; ++hIter2)
        {
            CaloHit *pCaloHit2 = *hIter2;

            if ((pCaloHit1 != pCaloHit2) && ((pCaloHit1->GetPositionVector() - pCaloHit2->GetPositionVector()).GetMagnitudeSquared() < m_maxHitDisplacementSquared))
            {
                isGoodHit = true;
                break;
            }
        }

        if (isGoodHit)
            matchedHits.insert(pCaloHit1);
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

// --- BEGIN EVENT DISPLAY ---
// std::cout << " MatchedSamplingPoints=" << fragmentOverlapResult.GetNMatchedSamplingPoints() << " MatchedFraction=" << fragmentOverlapResult.GetMatchedFraction() << " MatchedCaloHits=" << fragmentOverlapResult.GetFragmentCaloHitList().size() << std::endl;
// PandoraMonitoringApi::SetEveDisplayParameters(false, DETECTOR_VIEW_XZ);
// PandoraMonitoringApi::VisualizeCaloHits(&matchedHits, "MatchedHits", BLUE);
// for (unsigned int p=0; p < projectedPositions.size(); ++p)
// {
// CartesianVector projectedPosition = projectedPositions.at(p);
// PANDORA_MONITORING_API(AddMarkerToVisualization(&projectedPosition, "FitPosition", GREEN, 2.5));
// }
// PandoraMonitoringApi::ViewEvent();
// --- END EVENT DISPLAY

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

        LArClusterHelper::GetClusterSpanXZ(pCluster, minPosition, maxPosition);

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
    // ATTN: This method is currently mirrored in ClearTrackFragments tool

    if (overlapResult.GetMatchedFraction() < m_minMatchedSamplingPointFraction)
        return false;

    if (overlapResult.GetFragmentCaloHitList().size() < m_minMatchedHits)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDTrackFragmentsAlgorithm::BuildNewClusters(const CaloHitList &inputCaloHitList, ClusterList &newClusters) const
{
    // TODO: This code is cribbed from TwoDSlidingFitConsolidationAlgorithm, so create a common algorithm
    if (inputCaloHitList.empty())
        return;

    // Form simple associations between residual hits from deleted cluster
    HitToHitMap hitAssociationMap;

    for (CaloHitList::const_iterator iterI = inputCaloHitList.begin(), iterEndI = inputCaloHitList.end(); iterI != iterEndI; ++iterI)
    {
        CaloHit* pCaloHitI = *iterI;

        for (CaloHitList::const_iterator iterJ = iterI, iterEndJ = iterEndI; iterJ != iterEndJ; ++iterJ)
        {
            CaloHit* pCaloHitJ = *iterJ;

            if (pCaloHitI == pCaloHitJ)
                continue;

            if ((pCaloHitI->GetPositionVector() - pCaloHitJ->GetPositionVector()).GetMagnitudeSquared() < m_reclusteringWindow * m_reclusteringWindow)
            {
                hitAssociationMap[pCaloHitI].insert(pCaloHitJ);
                hitAssociationMap[pCaloHitJ].insert(pCaloHitI);
            }
        }
    }

    // Collect up associations and build new clusters
    CaloHitList vetoList;

    for (CaloHitList::const_iterator iterI = inputCaloHitList.begin(), iterEndI = inputCaloHitList.end(); iterI != iterEndI; ++iterI)
    {
        CaloHit* pSeedCaloHit = *iterI;

        if (vetoList.count(pSeedCaloHit))
            continue;

        CaloHitList mergeList;
        this->CollectAssociatedHits(pSeedCaloHit, pSeedCaloHit, hitAssociationMap, vetoList, mergeList);

        Cluster *pCluster = NULL;
        PandoraContentApi::Cluster::Parameters parameters;
        parameters.m_caloHitList.insert(pSeedCaloHit);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pCluster));
        newClusters.insert(pCluster);
        vetoList.insert(pSeedCaloHit);

        for (CaloHitList::const_iterator iterJ = mergeList.begin(), iterEndJ = mergeList.end(); iterJ != iterEndJ; ++iterJ)
        {
            CaloHit* pAssociatedCaloHit = *iterJ;

            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pCluster, pAssociatedCaloHit));
            vetoList.insert(pAssociatedCaloHit);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDTrackFragmentsAlgorithm::CollectAssociatedHits(CaloHit *pSeedCaloHit, CaloHit *pCurrentCaloHit,
    const HitToHitMap &hitAssociationMap, const CaloHitList &vetoList, CaloHitList &mergeList) const
{
    if (vetoList.count(pCurrentCaloHit))
        return;

    HitToHitMap::const_iterator iter1 = hitAssociationMap.find(pCurrentCaloHit);
    if (iter1 == hitAssociationMap.end())
        return;

    for (CaloHitList::const_iterator iter2 = iter1->second.begin(), iterEnd2 = iter1->second.end(); iter2 != iterEnd2; ++iter2)
    {
        CaloHit* pAssociatedCaloHit = *iter2;

        if (pAssociatedCaloHit == pSeedCaloHit)
            continue;

        if (!mergeList.insert(pAssociatedCaloHit).second)
            continue;

        this->CollectAssociatedHits(pSeedCaloHit, pAssociatedCaloHit, hitAssociationMap, vetoList, mergeList);
    }

    return;
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
    AlgorithmToolList algorithmToolList;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle,
        "TrackTools", algorithmToolList));

    for (AlgorithmToolList::const_iterator iter = algorithmToolList.begin(), iterEnd = algorithmToolList.end(); iter != iterEnd; ++iter)
    {
        FragmentTensorTool *pTensorManipulationTool(dynamic_cast<FragmentTensorTool*>(*iter));

        if (NULL == pTensorManipulationTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolList.push_back(pTensorManipulationTool);
    }

    m_nMaxTensorToolRepeats = 5000;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NMaxTensorToolRepeats", m_nMaxTensorToolRepeats));

    m_minXOverlap = 3.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinXOverlap", m_minXOverlap));

    m_minXOverlapFraction = 0.8f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinXOverlapFraction", m_minXOverlapFraction));

    m_maxPointDisplacement = 1.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxPointDisplacement", m_maxPointDisplacement));
    m_maxPointDisplacementSquared = m_maxPointDisplacement * m_maxPointDisplacement;

    float maxHitDisplacement = 5.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxHitDisplacement", maxHitDisplacement));
    m_maxHitDisplacementSquared = maxHitDisplacement * maxHitDisplacement;

    m_minMatchedSamplingPointFraction = 0.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedSamplingPointFraction", m_minMatchedSamplingPointFraction));

    m_minMatchedHits = 5;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedHits", m_minMatchedHits));

    m_reclusteringWindow = 2.5f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ReclusteringWindow", m_reclusteringWindow));

    return ThreeDTracksBaseAlgorithm<FragmentOverlapResult>::ReadSettings(xmlHandle);
}

} // namespace lar
