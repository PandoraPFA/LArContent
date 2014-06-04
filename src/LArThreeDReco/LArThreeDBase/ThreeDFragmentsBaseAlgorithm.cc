/**
 *  @file   LArContent/src/LArThreeDReco/LArThreeDBase/ThreeDFragmentsBaseAlgorithm.cc
 *
 *  @brief  Implementation of the three dimensional fragments algorithm base class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArThreeDHelper.h"

#include "LArThreeDReco/LArThreeDBase/ThreeDFragmentsBaseAlgorithm.h"

using namespace pandora;

namespace lar
{

void ThreeDFragmentsBaseAlgorithm::UpdateForNewCluster(Cluster *const pNewCluster)
{
    TwoDSlidingFitResult slidingFitResult;
    LArClusterHelper::LArTwoDSlidingFit(pNewCluster, m_slidingFitWindow, slidingFitResult);

    if (!m_slidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pNewCluster, slidingFitResult)).second)
        throw StatusCodeException(STATUS_CODE_FAILURE);

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

void ThreeDFragmentsBaseAlgorithm::PerformMainLoop()
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

void ThreeDFragmentsBaseAlgorithm::CalculateOverlapResult(Cluster *pClusterU, Cluster *pClusterV, Cluster *pClusterW)
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

    try
    {
        Cluster *pBestMatchedCluster(NULL);
        FragmentOverlapResult fragmentOverlapResult;
        this->CalculateOverlapResult(fitResult1, fitResult2, inputClusterList, pBestMatchedCluster, fragmentOverlapResult);

        if (fragmentOverlapResult.IsInitialized())
        {
            m_overlapTensor.SetOverlapResult((NULL != pClusterU) ? pClusterU : pBestMatchedCluster,
                (NULL != pClusterV) ? pClusterV : pBestMatchedCluster,
                (NULL != pClusterW) ? pClusterW : pBestMatchedCluster, fragmentOverlapResult);
        }
    }
    catch (StatusCodeException &statusCodeException)
    {
        if (!(STATUS_CODE_NOT_FOUND == statusCodeException.GetStatusCode() || STATUS_CODE_NOT_INITIALIZED == statusCodeException.GetStatusCode()))
            throw statusCodeException;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDFragmentsBaseAlgorithm::CalculateOverlapResult(const TwoDSlidingFitResult &fitResult1, const TwoDSlidingFitResult &fitResult2,
    const ClusterList &inputClusterList, Cluster *&pBestMatchedCluster, FragmentOverlapResult &fragmentOverlapResult) const
{
    const Cluster *pCluster1(fitResult1.GetCluster());
    const Cluster *pCluster2(fitResult2.GetCluster());

    float xMin1(0.f), xMax1(0.f), xMin2(0.f), xMax2(0.f);
    LArClusterHelper::GetClusterSpanX(pCluster1, xMin1, xMax1);
    LArClusterHelper::GetClusterSpanX(pCluster2, xMin2, xMax2);

    const float xOverlap(std::min(xMax1, xMax2) - std::max(xMin1, xMin2));
    const float xSpan(std::max(xMax1, xMax2) - std::min(xMin1, xMin2));

    if ((xOverlap < m_minXOverlap) || (xSpan < std::numeric_limits<float>::epsilon()) || ((xOverlap / xSpan) < m_minXOverlapFraction))
        return;

    CartesianPointList projectedPositions;
    this->GetProjectedPositions(fitResult1, fitResult2, projectedPositions);

    CaloHitList associatedHits;
    HitToClusterMap hitToClusterMap;
    this->GetAssociatedHits(inputClusterList, projectedPositions, hitToClusterMap, associatedHits);

    CaloHitList matchedHits;
    this->GetMatchedHits(associatedHits, matchedHits);

    ClusterList matchedClusters;
    this->GetMatchedClusters(matchedHits, hitToClusterMap, matchedClusters, pBestMatchedCluster);

    this->GetFragmentOverlapResult(projectedPositions, matchedHits, matchedClusters, fragmentOverlapResult);

    if (!this->PassesChecks(fitResult1, fitResult2, fragmentOverlapResult))
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDFragmentsBaseAlgorithm::GetAssociatedHits(const ClusterList &inputClusterList, const CartesianPointList &projectedPositions,
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
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDFragmentsBaseAlgorithm::GetMatchedHits(const CaloHitList &associatedHits, CaloHitList &matchedHits) const
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

void ThreeDFragmentsBaseAlgorithm::GetMatchedClusters(const CaloHitList &matchedHits, const HitToClusterMap &hitToClusterMap,
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

void ThreeDFragmentsBaseAlgorithm::GetFragmentOverlapResult(const CartesianPointList &projectedPositions, const CaloHitList &matchedHits,
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

    fragmentOverlapResult = FragmentOverlapResult(nMatchedSamplingPoints, m_nSamplingPoints, chi2Sum, matchedHits, matchedClusters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ThreeDFragmentsBaseAlgorithm::PassesChecks(const TwoDSlidingFitResult &fitResult1, const TwoDSlidingFitResult &fitResult2,
    FragmentOverlapResult &fragmentOverlapResult) const
{
    if (fragmentOverlapResult.GetNMatchedSamplingPoints() < m_minMatchedSamplingPoints)
        return false;

    if (fragmentOverlapResult.GetMatchedFraction() < m_minMatchedSamplingPointFraction)
        return false;

    if (fragmentOverlapResult.GetFragmentCaloHitList().size() < m_minMatchedHits)
        return false;

    const ClusterList &matchedClusters(fragmentOverlapResult.GetFragmentClusterList());

    if (matchedClusters.empty())
        return false;

    const Cluster *pCluster1(fitResult1.GetCluster());
    const Cluster *pCluster2(fitResult2.GetCluster());
    float xMin1(0.f), xMax1(0.f), xMin2(0.f), xMax2(0.f);
    LArClusterHelper::GetClusterSpanX(pCluster1, xMin1, xMax1);
    LArClusterHelper::GetClusterSpanX(pCluster2, xMin2, xMax2);

    for (ClusterList::const_iterator cIter = matchedClusters.begin(), cIterEnd = matchedClusters.end(); cIter != cIterEnd; ++cIter)
    {
        Cluster *pCluster = *cIter;

        float xMin(0.f), xMax(0.f);
        LArClusterHelper::GetClusterSpanX(pCluster, xMin, xMax);

        if ((xMax - xMin) > std::min((xMax1 - xMin1), (xMax2 - xMin2)))
            return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDFragmentsBaseAlgorithm::ExamineTensor()
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

StatusCode ThreeDFragmentsBaseAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
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

    m_nSamplingPoints = 100;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NSamplingPoints", m_nSamplingPoints));

    float maxPointDisplacement = 1.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxPointDisplacement", maxPointDisplacement));
    m_maxPointDisplacementSquared = maxPointDisplacement * maxPointDisplacement;

    float maxHitDisplacement = 5.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxHitDisplacement", maxHitDisplacement));
    m_maxHitDisplacementSquared = maxHitDisplacement * maxHitDisplacement;

    m_minMatchedSamplingPoints = 10;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedSamplingPoints", m_minMatchedSamplingPoints));

    m_minMatchedSamplingPointFraction = 0.6f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedSamplingPointFraction", m_minMatchedSamplingPointFraction));

    m_minMatchedHits = 5;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedHits", m_minMatchedHits));

    return ThreeDTracksBaseAlgorithm<FragmentOverlapResult>::ReadSettings(xmlHandle);
}

} // namespace lar
