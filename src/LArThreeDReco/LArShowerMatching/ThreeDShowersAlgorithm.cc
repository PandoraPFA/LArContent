/**
 *  @file   LArContent/src/LArThreeDReco/LArShowerMatching/ThreeDShowersAlgorithm.cc
 * 
 *  @brief  Implementation of the three dimensional showers algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"

#include "LArThreeDReco/LArShowerMatching/ThreeDShowersAlgorithm.h"

using namespace pandora;

namespace lar
{

const ThreeDShowersAlgorithm::SlidingShowerFitResult &ThreeDShowersAlgorithm::GetCachedSlidingFitResult(Cluster *const pCluster) const
{
    SlidingShowerFitResultMap::const_iterator iter = m_slidingFitResultMap.find(pCluster);

    if (m_slidingFitResultMap.end() == iter)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDShowersAlgorithm::UpdateForNewCluster(Cluster *const pNewCluster)
{
    this->AddToSlidingFitCache(pNewCluster);
    ThreeDBaseAlgorithm<ShowerOverlapResult>::UpdateForNewCluster(pNewCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDShowersAlgorithm::UpdateUponDeletion(Cluster *const pDeletedCluster)
{
    this->RemoveFromSlidingFitCache(pDeletedCluster);
    ThreeDBaseAlgorithm<ShowerOverlapResult>::UpdateUponDeletion(pDeletedCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDShowersAlgorithm::SelectInputClusters(const ClusterList *const pInputClusterList, ClusterList &selectedClusterList) const
{
    for (ClusterList::const_iterator iter = pInputClusterList->begin(), iterEnd = pInputClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (m_ignoreUnavailableClusters && !pCluster->IsAvailable())
            continue;

        if (pCluster->GetNCaloHits() < m_minClusterCaloHits)
            continue;

        if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLengthSquared)
            continue;

        selectedClusterList.insert(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDShowersAlgorithm::PreparationStep()
{
    ClusterList allClustersList;
    allClustersList.insert(this->m_clusterListU.begin(), this->m_clusterListU.end());
    allClustersList.insert(this->m_clusterListV.begin(), this->m_clusterListV.end());
    allClustersList.insert(this->m_clusterListW.begin(), this->m_clusterListW.end());

    for (ClusterList::const_iterator iter = allClustersList.begin(), iterEnd = allClustersList.end(); iter != iterEnd; ++iter)
    {
        this->AddToSlidingFitCache(*iter);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDShowersAlgorithm::TidyUp()
{
    m_slidingFitResultMap.clear();
    return ThreeDBaseAlgorithm<ShowerOverlapResult>::TidyUp();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDShowersAlgorithm::AddToSlidingFitCache(Cluster *const pCluster)
{
    SlidingShowerFitResult slidingShowerFitResult;
    LArClusterHelper::LArTwoDSlidingFit(pCluster, m_slidingFitWindow, slidingShowerFitResult.m_showerFitResult);
    LArClusterHelper::LArTwoDShowerEdgeFit(slidingShowerFitResult.m_showerFitResult, NEGATIVE_SHOWER_EDGE, slidingShowerFitResult.m_negativeEdgeFitResult);
    LArClusterHelper::LArTwoDShowerEdgeFit(slidingShowerFitResult.m_showerFitResult, POSITIVE_SHOWER_EDGE, slidingShowerFitResult.m_positiveEdgeFitResult);

    if (!m_slidingFitResultMap.insert(SlidingShowerFitResultMap::value_type(pCluster, slidingShowerFitResult)).second)
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDShowersAlgorithm::RemoveFromSlidingFitCache(Cluster *const pCluster)
{
    SlidingShowerFitResultMap::iterator iter = m_slidingFitResultMap.find(pCluster);

    if (m_slidingFitResultMap.end() != iter)
        m_slidingFitResultMap.erase(iter);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDShowersAlgorithm::CalculateOverlapResult(Cluster *pClusterU, Cluster *pClusterV, Cluster *pClusterW)
{
    try
    {
        ShowerOverlapResult overlapResult;
        this->CalculateOverlapResult(pClusterU, pClusterV, pClusterW, overlapResult);

        if (overlapResult.IsInitialized())
            m_overlapTensor.SetOverlapResult(pClusterU, pClusterV, pClusterW, overlapResult);
    }
    catch (StatusCodeException &statusCodeException)
    {
        if (!(STATUS_CODE_NOT_FOUND == statusCodeException.GetStatusCode() || STATUS_CODE_NOT_INITIALIZED == statusCodeException.GetStatusCode()))
            throw statusCodeException;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDShowersAlgorithm::CalculateOverlapResult(Cluster *pClusterU, Cluster *pClusterV, Cluster *pClusterW, ShowerOverlapResult &overlapResult)
{
    const SlidingShowerFitResult &fitResultU(this->GetCachedSlidingFitResult(pClusterU));
    const SlidingShowerFitResult &fitResultV(this->GetCachedSlidingFitResult(pClusterV));
    const SlidingShowerFitResult &fitResultW(this->GetCachedSlidingFitResult(pClusterW));

    const XSampling xSampling(fitResultU.m_showerFitResult, fitResultV.m_showerFitResult, fitResultW.m_showerFitResult);

    ShowerPositionMap uMap, vMap, wMap;
    this->GetShowerPositionMaps(fitResultU.m_showerFitResult, fitResultV.m_showerFitResult, fitResultW.m_showerFitResult, xSampling, uMap, vMap, wMap);

    ShowerPositionMap uPosMap, vPosMap, wPosMap;
    this->GetShowerPositionMaps(fitResultU.m_positiveEdgeFitResult, fitResultV.m_positiveEdgeFitResult, fitResultW.m_positiveEdgeFitResult, xSampling, uPosMap, vPosMap, wPosMap);

    ShowerPositionMap uNegMap, vNegMap, wNegMap;
    this->GetShowerPositionMaps(fitResultU.m_negativeEdgeFitResult, fitResultV.m_negativeEdgeFitResult, fitResultW.m_negativeEdgeFitResult, xSampling, uNegMap, vNegMap, wNegMap);

    unsigned int nSampledHitsU(0), nMatchedHitsU(0);
    this->GetHitOverlapFraction(pClusterU, xSampling, uPosMap, uNegMap, nSampledHitsU, nMatchedHitsU);

    unsigned int nSampledHitsV(0), nMatchedHitsV(0);
    this->GetHitOverlapFraction(pClusterV, xSampling, vPosMap, vNegMap, nSampledHitsV, nMatchedHitsV);

    unsigned int nSampledHitsW(0), nMatchedHitsW(0);
    this->GetHitOverlapFraction(pClusterW, xSampling, wPosMap, wNegMap, nSampledHitsW, nMatchedHitsW);

    const unsigned int nMatchedHits(nMatchedHitsU + nMatchedHitsV + nMatchedHitsW);
    const unsigned int nSampledHits(nSampledHitsU + nSampledHitsV + nSampledHitsW);

    if (0 == nSampledHits)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    const XOverlap xOverlapObject(xSampling.m_uMinX, xSampling.m_uMaxX, xSampling.m_vMinX, xSampling.m_vMaxX, xSampling.m_wMinX, xSampling.m_wMaxX, xSampling.m_xOverlapSpan);
    const ShowerOverlapResult showerOverlapResult(nMatchedHits, nSampledHits, xOverlapObject);

    if ((showerOverlapResult.GetMatchedFraction() < m_minShowerMatchedFraction) || (showerOverlapResult.GetNMatchedSamplingPoints() < m_minShowerMatchedPoints))
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    overlapResult = showerOverlapResult;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDShowersAlgorithm::GetShowerPositionMaps(const TwoDSlidingFitResult &fitResultU, const TwoDSlidingFitResult &fitResultV,
    const TwoDSlidingFitResult &fitResultW, const XSampling &xSampling, ShowerPositionMap &positionMapU, ShowerPositionMap &positionMapV,
    ShowerPositionMap &positionMapW) const
{
    for (float x = xSampling.m_minX; x < xSampling.m_maxX; x += xSampling.m_xPitch)
    {
        try
        {
            CartesianVector fitUVector(0.f, 0.f, 0.f), fitVVector(0.f, 0.f, 0.f), fitWVector(0.f, 0.f, 0.f);
            fitResultU.GetGlobalFitPositionAtX(x, fitUVector);
            fitResultV.GetGlobalFitPositionAtX(x, fitVVector);
            fitResultW.GetGlobalFitPositionAtX(x, fitWVector);

            const float u(fitUVector.GetZ()), v(fitVVector.GetZ()), w(fitWVector.GetZ());
            const float uv2w(LArGeometryHelper::MergeTwoPositions(TPC_VIEW_U, TPC_VIEW_V, u, v));
            const float uw2v(LArGeometryHelper::MergeTwoPositions(TPC_VIEW_U, TPC_VIEW_W, u, w));
            const float vw2u(LArGeometryHelper::MergeTwoPositions(TPC_VIEW_V, TPC_VIEW_W, v, w));

            const unsigned int xBin((x - xSampling.m_minX) / xSampling.m_xPitch);
            positionMapU.insert(ShowerPositionMap::value_type(xBin, CartesianVector(x, 0., vw2u)));
            positionMapV.insert(ShowerPositionMap::value_type(xBin, CartesianVector(x, 0., uw2v)));
            positionMapW.insert(ShowerPositionMap::value_type(xBin, CartesianVector(x, 0., uv2w)));
        }
        catch (StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDShowersAlgorithm::GetHitOverlapFraction(const Cluster *const pCluster, const XSampling &xSampling, const ShowerPositionMap &edgeMap1,
    const ShowerPositionMap &edgeMap2, unsigned int &nSampledHits, unsigned int &nMatchedHits) const
{
    if (((xSampling.m_maxX - xSampling.m_minX) < std::numeric_limits<float>::epsilon()) || (xSampling.m_xPitch < std::numeric_limits<float>::epsilon()))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    nSampledHits = 0; nMatchedHits = 0;
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hIter = iter->second->begin(), hIterEnd = iter->second->end(); hIter != hIterEnd; ++hIter)
        {
            CaloHit *pCaloHit = *hIter;
            const float x(pCaloHit->GetPositionVector().GetX());
            const float z(pCaloHit->GetPositionVector().GetZ());

            if ((x < xSampling.m_minX) || (x > xSampling.m_maxX))
                continue;

            ++nSampledHits;
            const unsigned int xBin((x - xSampling.m_minX) / xSampling.m_xPitch);

            ShowerPositionMap::const_iterator edgeIter1 = edgeMap1.find(xBin);
            ShowerPositionMap::const_iterator edgeIter2 = edgeMap2.find(xBin);

            if ((edgeMap1.end() == edgeIter1) || (edgeMap2.end() == edgeIter2))
                continue;

            const float maxZ(std::max(edgeIter1->second.GetZ(), edgeIter2->second.GetZ()));
            const float minZ(std::min(edgeIter1->second.GetZ(), edgeIter2->second.GetZ()));

            if ((z > minZ) && (z < maxZ))
            {
                ++nMatchedHits;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDShowersAlgorithm::ExamineTensor()
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
//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDShowersAlgorithm::XSampling::XSampling(const TwoDSlidingFitResult &fitResultU, const TwoDSlidingFitResult &fitResultV,
        const TwoDSlidingFitResult &fitResultW) :
    m_uMinX(std::min(fitResultU.GetGlobalMinLayerPosition().GetX(), fitResultU.GetGlobalMaxLayerPosition().GetX())),
    m_uMaxX(std::max(fitResultU.GetGlobalMinLayerPosition().GetX(), fitResultU.GetGlobalMaxLayerPosition().GetX())),
    m_vMinX(std::min(fitResultV.GetGlobalMinLayerPosition().GetX(), fitResultV.GetGlobalMaxLayerPosition().GetX())),
    m_vMaxX(std::max(fitResultV.GetGlobalMinLayerPosition().GetX(), fitResultV.GetGlobalMaxLayerPosition().GetX())),
    m_wMinX(std::min(fitResultW.GetGlobalMinLayerPosition().GetX(), fitResultW.GetGlobalMaxLayerPosition().GetX())),
    m_wMaxX(std::max(fitResultW.GetGlobalMinLayerPosition().GetX(), fitResultW.GetGlobalMaxLayerPosition().GetX())),
    m_minX(std::max(m_uMinX, std::max(m_vMinX, m_wMinX))),
    m_maxX(std::min(m_uMaxX, std::min(m_vMaxX, m_wMaxX))),
    m_xOverlapSpan(m_maxX - m_minX)
{
    const float nPointsU(std::fabs((m_xOverlapSpan / (m_uMaxX - m_uMinX)) * static_cast<float>(fitResultU.GetMaxLayer() - fitResultU.GetMinLayer())));
    const float nPointsV(std::fabs((m_xOverlapSpan / (m_vMaxX - m_vMinX)) * static_cast<float>(fitResultV.GetMaxLayer() - fitResultV.GetMinLayer())));
    const float nPointsW(std::fabs((m_xOverlapSpan / (m_wMaxX - m_wMinX)) * static_cast<float>(fitResultW.GetMaxLayer() - fitResultW.GetMinLayer())));
    const float nPoints(nPointsU + nPointsV + nPointsW);

    if ((m_xOverlapSpan < std::numeric_limits<float>::epsilon()) || (nPoints < std::numeric_limits<float>::epsilon()))
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    m_xPitch = 3.f * m_xOverlapSpan / nPoints;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDShowersAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    AlgorithmToolList algorithmToolList;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle,
        "ShowerTools", algorithmToolList));

    for (AlgorithmToolList::const_iterator iter = algorithmToolList.begin(), iterEnd = algorithmToolList.end(); iter != iterEnd; ++iter)
    {
        ShowerTensorTool *pShowerTensorTool(dynamic_cast<ShowerTensorTool*>(*iter));

        if (NULL == pShowerTensorTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolList.push_back(pShowerTensorTool);
    }

    m_nMaxTensorToolRepeats = 5000;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NMaxTensorToolRepeats", m_nMaxTensorToolRepeats));

    m_slidingFitWindow = 20;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));

    m_ignoreUnavailableClusters = false;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "IgnoreUnavailableClusters", m_ignoreUnavailableClusters));

    m_minClusterCaloHits = 5;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterCaloHits", m_minClusterCaloHits));

    float minClusterLength = 3.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLength", minClusterLength));
    m_minClusterLengthSquared = minClusterLength * minClusterLength;

    m_minShowerMatchedFraction = 0.2f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinShowerMatchedFraction", m_minShowerMatchedFraction));

    m_minShowerMatchedPoints = 20;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinShowerMatchedPoints", m_minShowerMatchedPoints));

    return ThreeDBaseAlgorithm<ShowerOverlapResult>::ReadSettings(xmlHandle);
}

} // namespace lar
