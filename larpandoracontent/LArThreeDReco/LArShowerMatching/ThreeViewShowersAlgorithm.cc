/**
 *  @file   larpandoracontent/LArThreeDReco/LArShowerMatching/ThreeViewShowersAlgorithm.cc
 *
 *  @brief  Implementation of the three view showers algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArThreeDReco/LArShowerMatching/ThreeViewShowersAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ThreeViewShowersAlgorithm::ThreeViewShowersAlgorithm() :
    m_nMaxTensorToolRepeats(1000),
    m_slidingFitWindow(20),
    m_ignoreUnavailableClusters(true),
    m_minClusterCaloHits(5),
    m_minClusterLengthSquared(3.f * 3.f),
    m_minShowerMatchedFraction(0.2f),
    m_minShowerMatchedPoints(20),
    m_visualize(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

const TwoDSlidingShowerFitResult &ThreeViewShowersAlgorithm::GetCachedSlidingFitResult(const Cluster *const pCluster) const
{
    TwoDSlidingShowerFitResultMap::const_iterator iter = m_slidingFitResultMap.find(pCluster);

    if (m_slidingFitResultMap.end() == iter)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewShowersAlgorithm::UpdateForNewCluster(const Cluster *const pNewCluster)
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

    BaseAlgorithm::UpdateForNewCluster(pNewCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewShowersAlgorithm::UpdateUponDeletion(const Cluster *const pDeletedCluster)
{
    this->RemoveFromSlidingFitCache(pDeletedCluster);
    BaseAlgorithm::UpdateUponDeletion(pDeletedCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewShowersAlgorithm::SelectInputClusters(const ClusterList *const pInputClusterList, ClusterList &selectedClusterList) const
{
    for (ClusterList::const_iterator iter = pInputClusterList->begin(), iterEnd = pInputClusterList->end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pCluster = *iter;

        if (m_ignoreUnavailableClusters && !pCluster->IsAvailable())
            continue;

        if (pCluster->GetNCaloHits() < m_minClusterCaloHits)
            continue;

        if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLengthSquared)
            continue;

        selectedClusterList.push_back(pCluster);
    }
    if (m_visualize)
    {
        PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &selectedClusterList, "Selected Clusters", RED));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewShowersAlgorithm::PrepareInputClusters(ClusterList &preparedClusterList)
{
    for (ClusterList::iterator iter = preparedClusterList.begin(), iterEnd = preparedClusterList.end(); iter != iterEnd;)
    {
        const Cluster *const pCluster(*iter);

        try
        {
            this->AddToSlidingFitCache(pCluster);
            ++iter;
        }
        catch (StatusCodeException &statusCodeException)
        {
            preparedClusterList.erase(iter++);

            if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                throw statusCodeException;
        }
    }
    if (m_visualize)
    {
        PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &preparedClusterList, "Prepared Clusters", GREEN));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewShowersAlgorithm::TidyUp()
{
    m_slidingFitResultMap.clear();
    return BaseAlgorithm::TidyUp();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewShowersAlgorithm::AddToSlidingFitCache(const Cluster *const pCluster)
{
    const float slidingFitPitch(LArGeometryHelper::GetWirePitch(this->GetPandora(), LArClusterHelper::GetClusterHitType(pCluster)));
    const TwoDSlidingShowerFitResult slidingShowerFitResult(pCluster, m_slidingFitWindow, slidingFitPitch);

    if (!m_slidingFitResultMap.insert(TwoDSlidingShowerFitResultMap::value_type(pCluster, slidingShowerFitResult)).second)
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewShowersAlgorithm::RemoveFromSlidingFitCache(const Cluster *const pCluster)
{
    TwoDSlidingShowerFitResultMap::iterator iter = m_slidingFitResultMap.find(pCluster);

    if (m_slidingFitResultMap.end() != iter)
        m_slidingFitResultMap.erase(iter);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewShowersAlgorithm::CalculateOverlapResult(const Cluster *const pClusterU, const Cluster *const pClusterV, const Cluster *const pClusterW)
{
    ShowerOverlapResult overlapResult;
    PANDORA_THROW_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, this->CalculateOverlapResult(pClusterU, pClusterV, pClusterW, overlapResult));

    if (overlapResult.IsInitialized())
        this->GetMatchingControl().GetOverlapTensor().SetOverlapResult(pClusterU, pClusterV, pClusterW, overlapResult);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeViewShowersAlgorithm::CalculateOverlapResult(
    const Cluster *const pClusterU, const Cluster *const pClusterV, const Cluster *const pClusterW, ShowerOverlapResult &overlapResult)
{
    const TwoDSlidingShowerFitResult &fitResultU(this->GetCachedSlidingFitResult(pClusterU));
    const TwoDSlidingShowerFitResult &fitResultV(this->GetCachedSlidingFitResult(pClusterV));
    const TwoDSlidingShowerFitResult &fitResultW(this->GetCachedSlidingFitResult(pClusterW));

    const XSampling xSampling(fitResultU.GetShowerFitResult(), fitResultV.GetShowerFitResult(), fitResultW.GetShowerFitResult());

    if (xSampling.m_xOverlapSpan < std::numeric_limits<float>::epsilon())
        return STATUS_CODE_NOT_FOUND;

    ShowerPositionMapPair positionMapsU, positionMapsV, positionMapsW;
    this->GetShowerPositionMaps(fitResultU, fitResultV, fitResultW, xSampling, positionMapsU, positionMapsV, positionMapsW);

    unsigned int nSampledHitsU(0), nMatchedHitsU(0);
    this->GetBestHitOverlapFraction(pClusterU, xSampling, positionMapsU, nSampledHitsU, nMatchedHitsU);

    unsigned int nSampledHitsV(0), nMatchedHitsV(0);
    this->GetBestHitOverlapFraction(pClusterV, xSampling, positionMapsV, nSampledHitsV, nMatchedHitsV);

    unsigned int nSampledHitsW(0), nMatchedHitsW(0);
    this->GetBestHitOverlapFraction(pClusterW, xSampling, positionMapsW, nSampledHitsW, nMatchedHitsW);

    const unsigned int nMatchedHits(nMatchedHitsU + nMatchedHitsV + nMatchedHitsW);
    const unsigned int nSampledHits(nSampledHitsU + nSampledHitsV + nSampledHitsW);

    if (0 == nSampledHits)
        return STATUS_CODE_NOT_FOUND;

    const XOverlap xOverlapObject(xSampling.m_uMinX, xSampling.m_uMaxX, xSampling.m_vMinX, xSampling.m_vMaxX, xSampling.m_wMinX,
        xSampling.m_wMaxX, xSampling.m_xOverlapSpan);
    const ShowerOverlapResult showerOverlapResult(nMatchedHits, nSampledHits, xOverlapObject);

    if ((showerOverlapResult.GetMatchedFraction() < m_minShowerMatchedFraction) || (showerOverlapResult.GetNMatchedSamplingPoints() < m_minShowerMatchedPoints))
        return STATUS_CODE_NOT_FOUND;

    if (m_visualize)
    {
        ClusterList clusterList{pClusterU, pClusterV, pClusterW};
        PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clusterList, "Overlapping Clusters", BLUE));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    overlapResult = showerOverlapResult;
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewShowersAlgorithm::GetShowerPositionMaps(const TwoDSlidingShowerFitResult &fitResultU,
    const TwoDSlidingShowerFitResult &fitResultV, const TwoDSlidingShowerFitResult &fitResultW, const XSampling &xSampling,
    ShowerPositionMapPair &positionMapsU, ShowerPositionMapPair &positionMapsV, ShowerPositionMapPair &positionMapsW) const
{
    const unsigned int nPoints(static_cast<unsigned int>(xSampling.m_nPoints));

    for (unsigned n = 0; n <= nPoints; ++n)
    {
        const float x(xSampling.m_minX + (xSampling.m_maxX - xSampling.m_minX) * static_cast<float>(n) / static_cast<float>(nPoints));

        int xBin(-1);
        if (STATUS_CODE_SUCCESS != xSampling.GetBin(x, xBin))
            continue;

        FloatVector uValues, vValues, wValues;
        fitResultU.GetShowerEdges(x, true, uValues);
        fitResultV.GetShowerEdges(x, true, vValues);
        fitResultW.GetShowerEdges(x, true, wValues);

        std::sort(uValues.begin(), uValues.end());
        std::sort(vValues.begin(), vValues.end());
        std::sort(wValues.begin(), wValues.end());

        if ((uValues.size() > 1) && (vValues.size() > 1))
        {
            const float uMin(uValues.front()), uMax(uValues.back());
            const float vMin(vValues.front()), vMax(vValues.back());
            const float uv2wMinMin(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, uMin, vMin));
            const float uv2wMaxMax(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, uMax, vMax));
            const float uv2wMinMax(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, uMin, vMax));
            const float uv2wMaxMin(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, uMax, vMin));
            positionMapsW.first.insert(ShowerPositionMap::value_type(xBin, ShowerExtent(x, uv2wMinMin, uv2wMaxMax)));
            positionMapsW.second.insert(ShowerPositionMap::value_type(xBin, ShowerExtent(x, uv2wMinMax, uv2wMaxMin)));
        }

        if ((uValues.size() > 1) && (wValues.size() > 1))
        {
            const float uMin(uValues.front()), uMax(uValues.back());
            const float wMin(wValues.front()), wMax(wValues.back());
            const float uw2vMinMin(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_W, uMin, wMin));
            const float uw2vMaxMax(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_W, uMax, wMax));
            const float uw2vMinMax(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_W, uMin, wMax));
            const float uw2vMaxMin(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_W, uMax, wMin));
            positionMapsV.first.insert(ShowerPositionMap::value_type(xBin, ShowerExtent(x, uw2vMinMin, uw2vMaxMax)));
            positionMapsV.second.insert(ShowerPositionMap::value_type(xBin, ShowerExtent(x, uw2vMinMax, uw2vMaxMin)));
        }

        if ((vValues.size() > 1) && (wValues.size() > 1))
        {
            const float vMin(vValues.front()), vMax(vValues.back());
            const float wMin(wValues.front()), wMax(wValues.back());
            const float vw2uMinMin(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_V, TPC_VIEW_W, vMin, wMin));
            const float vw2uMaxMax(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_V, TPC_VIEW_W, vMax, wMax));
            const float vw2uMinMax(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_V, TPC_VIEW_W, vMin, wMax));
            const float vw2uMaxMin(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_V, TPC_VIEW_W, vMax, wMin));
            positionMapsU.first.insert(ShowerPositionMap::value_type(xBin, ShowerExtent(x, vw2uMinMin, vw2uMaxMax)));
            positionMapsU.second.insert(ShowerPositionMap::value_type(xBin, ShowerExtent(x, vw2uMinMax, vw2uMaxMin)));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewShowersAlgorithm::GetBestHitOverlapFraction(const Cluster *const pCluster, const XSampling &xSampling,
    const ShowerPositionMapPair &positionMaps, unsigned int &nSampledHits, unsigned int &nMatchedHits) const
{
    if ((xSampling.m_maxX - xSampling.m_minX) < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    nSampledHits = 0;
    nMatchedHits = 0;
    unsigned int nMatchedHits1(0), nMatchedHits2(0);
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hIter = iter->second->begin(), hIterEnd = iter->second->end(); hIter != hIterEnd; ++hIter)
        {
            const CaloHit *const pCaloHit = *hIter;
            const float x(pCaloHit->GetPositionVector().GetX());
            const float z(pCaloHit->GetPositionVector().GetZ());

            int xBin(-1);
            if (STATUS_CODE_SUCCESS != xSampling.GetBin(x, xBin))
                continue;

            ++nSampledHits;

            ShowerPositionMap::const_iterator positionIter1 = positionMaps.first.find(xBin);
            ShowerPositionMap::const_iterator positionIter2 = positionMaps.second.find(xBin);

            if ((positionMaps.first.end() != positionIter1) && (z > positionIter1->second.GetLowEdgeZ()) && (z < positionIter1->second.GetHighEdgeZ()))
                ++nMatchedHits1;

            if ((positionMaps.second.end() != positionIter2) && (z > positionIter2->second.GetLowEdgeZ()) && (z < positionIter2->second.GetHighEdgeZ()))
                ++nMatchedHits2;
        }
    }

    nMatchedHits = std::max(nMatchedHits1, nMatchedHits2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewShowersAlgorithm::ExamineOverlapContainer()
{
    unsigned int repeatCounter(0);

    for (TensorToolVector::const_iterator iter = m_algorithmToolVector.begin(), iterEnd = m_algorithmToolVector.end(); iter != iterEnd;)
    {
        if ((*iter)->Run(this, this->GetMatchingControl().GetOverlapTensor()))
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
//------------------------------------------------------------------------------------------------------------------------------------------

ThreeViewShowersAlgorithm::XSampling::XSampling(
    const TwoDSlidingFitResult &fitResultU, const TwoDSlidingFitResult &fitResultV, const TwoDSlidingFitResult &fitResultW)
{
    fitResultU.GetMinAndMaxX(m_uMinX, m_uMaxX);
    fitResultV.GetMinAndMaxX(m_vMinX, m_vMaxX);
    fitResultW.GetMinAndMaxX(m_wMinX, m_wMaxX);
    m_minX = std::max(m_uMinX, std::max(m_vMinX, m_wMinX));
    m_maxX = std::min(m_uMaxX, std::min(m_vMaxX, m_wMaxX));
    m_xOverlapSpan = (m_maxX - m_minX);
    m_nPoints = 1.f;

    if (m_xOverlapSpan > std::numeric_limits<float>::epsilon())
    {
        const float nPointsU(
            std::fabs((m_xOverlapSpan / (m_uMaxX - m_uMinX)) * static_cast<float>(fitResultU.GetMaxLayer() - fitResultU.GetMinLayer())));
        const float nPointsV(
            std::fabs((m_xOverlapSpan / (m_vMaxX - m_vMinX)) * static_cast<float>(fitResultV.GetMaxLayer() - fitResultV.GetMinLayer())));
        const float nPointsW(
            std::fabs((m_xOverlapSpan / (m_wMaxX - m_wMinX)) * static_cast<float>(fitResultW.GetMaxLayer() - fitResultW.GetMinLayer())));
        m_nPoints = 1.f + ((nPointsU + nPointsV + nPointsW) / 3.f);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeViewShowersAlgorithm::XSampling::GetBin(const float x, int &xBin) const
{
    if (((x - m_minX) < -std::numeric_limits<float>::epsilon()) || ((x - m_maxX) > +std::numeric_limits<float>::epsilon()))
        return STATUS_CODE_NOT_FOUND;

    xBin = static_cast<int>(0.5f + m_nPoints * (x - m_minX) / (m_maxX - m_minX));
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeViewShowersAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "ShowerTools", algorithmToolVector));

    for (AlgorithmToolVector::const_iterator iter = algorithmToolVector.begin(), iterEnd = algorithmToolVector.end(); iter != iterEnd; ++iter)
    {
        ShowerTensorTool *const pShowerTensorTool(dynamic_cast<ShowerTensorTool *>(*iter));

        if (!pShowerTensorTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolVector.push_back(pShowerTensorTool);
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NMaxTensorToolRepeats", m_nMaxTensorToolRepeats));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "IgnoreUnavailableClusters", m_ignoreUnavailableClusters));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterCaloHits", m_minClusterCaloHits));

    float minClusterLength = std::sqrt(m_minClusterLengthSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterLength", minClusterLength));
    m_minClusterLengthSquared = minClusterLength * minClusterLength;

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinShowerMatchedFraction", m_minShowerMatchedFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinShowerMatchedPoints", m_minShowerMatchedPoints));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualize", m_visualize));

    return BaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
