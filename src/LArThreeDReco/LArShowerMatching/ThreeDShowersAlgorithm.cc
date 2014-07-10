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

void ThreeDShowersAlgorithm::SelectInputClusters(const ClusterList *const pInputClusterList, ClusterList &selectedClusterList) const
{
    for (ClusterList::const_iterator iter = pInputClusterList->begin(), iterEnd = pInputClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (!pCluster->IsAvailable())
            continue;

        if (pCluster->GetNCaloHits() < m_minClusterCaloHits)
            continue;

        if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLengthSquared)
            continue;

        selectedClusterList.insert(pCluster);
    }
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
    TwoDSlidingFitResult slidingFitResultU, slidingFitResultV, slidingFitResultW;
    LArClusterHelper::LArTwoDSlidingFit(pClusterU, m_slidingFitWindow, slidingFitResultU);
    LArClusterHelper::LArTwoDSlidingFit(pClusterV, m_slidingFitWindow, slidingFitResultV);
    LArClusterHelper::LArTwoDSlidingFit(pClusterW, m_slidingFitWindow, slidingFitResultW);

    // Assess x-overlap
    const float minXU(std::min(slidingFitResultU.GetGlobalMinLayerPosition().GetX(), slidingFitResultU.GetGlobalMaxLayerPosition().GetX()));
    const float maxXU(std::max(slidingFitResultU.GetGlobalMinLayerPosition().GetX(), slidingFitResultU.GetGlobalMaxLayerPosition().GetX()));
    const float xSpanU(maxXU - minXU);

    const float minXV(std::min(slidingFitResultV.GetGlobalMinLayerPosition().GetX(), slidingFitResultV.GetGlobalMaxLayerPosition().GetX()));
    const float maxXV(std::max(slidingFitResultV.GetGlobalMinLayerPosition().GetX(), slidingFitResultV.GetGlobalMaxLayerPosition().GetX()));
    const float xSpanV(maxXV - minXV);

    const float minXW(std::min(slidingFitResultW.GetGlobalMinLayerPosition().GetX(), slidingFitResultW.GetGlobalMaxLayerPosition().GetX()));
    const float maxXW(std::max(slidingFitResultW.GetGlobalMinLayerPosition().GetX(), slidingFitResultW.GetGlobalMaxLayerPosition().GetX()));
    const float xSpanW(maxXW - minXW);

    const float minX(std::max(minXU, std::max(minXV, minXW)));
    const float maxX(std::min(maxXU, std::min(maxXV, maxXW)));
    const float xOverlap(maxX - minX);

    if (xOverlap < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    // Sampling in x
    const float nPointsU(std::fabs((xOverlap / xSpanU) * static_cast<float>(slidingFitResultU.GetMaxLayer() - slidingFitResultU.GetMinLayer())));
    const float nPointsV(std::fabs((xOverlap / xSpanV) * static_cast<float>(slidingFitResultV.GetMaxLayer() - slidingFitResultV.GetMinLayer())));
    const float nPointsW(std::fabs((xOverlap / xSpanW) * static_cast<float>(slidingFitResultW.GetMaxLayer() - slidingFitResultW.GetMinLayer())));
    const float xPitch(3.f * xOverlap / (nPointsU + nPointsV + nPointsW));

    ShowerPositionMap uMap, vMap, wMap;
    this->GetShowerPositionMaps(slidingFitResultU, slidingFitResultV, slidingFitResultW, minX, maxX, xPitch, uMap, vMap, wMap);

    ShowerPositionMap uPosMap, vPosMap, wPosMap;
    TwoDSlidingFitResult posShowerEdgeFitU, posShowerEdgeFitV, posShowerEdgeFitW;
    LArClusterHelper::LArTwoDShowerEdgeFit(slidingFitResultU, POSITIVE_SHOWER_EDGE, posShowerEdgeFitU);
    LArClusterHelper::LArTwoDShowerEdgeFit(slidingFitResultV, POSITIVE_SHOWER_EDGE, posShowerEdgeFitV);
    LArClusterHelper::LArTwoDShowerEdgeFit(slidingFitResultW, POSITIVE_SHOWER_EDGE, posShowerEdgeFitW);
    this->GetShowerPositionMaps(posShowerEdgeFitU, posShowerEdgeFitV, posShowerEdgeFitW, minX, maxX, xPitch, uPosMap, vPosMap, wPosMap);

    ShowerPositionMap uNegMap, vNegMap, wNegMap;
    TwoDSlidingFitResult negShowerEdgeFitU, negShowerEdgeFitV, negShowerEdgeFitW;
    LArClusterHelper::LArTwoDShowerEdgeFit(slidingFitResultU, NEGATIVE_SHOWER_EDGE, negShowerEdgeFitU);
    LArClusterHelper::LArTwoDShowerEdgeFit(slidingFitResultV, NEGATIVE_SHOWER_EDGE, negShowerEdgeFitV);
    LArClusterHelper::LArTwoDShowerEdgeFit(slidingFitResultW, NEGATIVE_SHOWER_EDGE, negShowerEdgeFitW);
    this->GetShowerPositionMaps(negShowerEdgeFitU, negShowerEdgeFitV, negShowerEdgeFitW, minX, maxX, xPitch, uNegMap, vNegMap, wNegMap);

    // Assess hit overlap fractions
    unsigned int nSampledHitsU(0), nMatchedHitsU(0);
    this->GetHitOverlapFraction(pClusterU, minX, maxX, xPitch, uPosMap, uNegMap, nSampledHitsU, nMatchedHitsU);

    unsigned int nSampledHitsV(0), nMatchedHitsV(0);
    this->GetHitOverlapFraction(pClusterV, minX, maxX, xPitch, vPosMap, vNegMap, nSampledHitsV, nMatchedHitsV);

    unsigned int nSampledHitsW(0), nMatchedHitsW(0);
    this->GetHitOverlapFraction(pClusterW, minX, maxX, xPitch, wPosMap, wNegMap, nSampledHitsW, nMatchedHitsW);

    const unsigned int nMatchedHits(nMatchedHitsU + nMatchedHitsV + nMatchedHitsW);
    const unsigned int nSampledHits(nSampledHitsU + nSampledHitsV + nSampledHitsW);

    if (0 == nSampledHits)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    const ShowerOverlapResult showerOverlapResult(nMatchedHits, nSampledHits);

    if ((showerOverlapResult.GetMatchedFraction() < m_minShowerMatchedFraction) || (showerOverlapResult.GetNMatchedSamplingPoints() < m_minShowerMatchedPoints))
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    overlapResult = showerOverlapResult;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDShowersAlgorithm::GetShowerPositionMaps(const TwoDSlidingFitResult &slidingFitResultU, const TwoDSlidingFitResult &slidingFitResultV,
    const TwoDSlidingFitResult &slidingFitResultW, const float minX, const float maxX, const float xPitch, ShowerPositionMap &showerPositionMapU,
    ShowerPositionMap &showerPositionMapV, ShowerPositionMap &showerPositionMapW) const
{
    for (float x = minX; x < maxX; x += xPitch)
    {
        try
        {
            CartesianVector fitUVector(0.f, 0.f, 0.f), fitVVector(0.f, 0.f, 0.f), fitWVector(0.f, 0.f, 0.f);
            slidingFitResultU.GetGlobalFitPositionAtX(x, fitUVector);
            slidingFitResultV.GetGlobalFitPositionAtX(x, fitVVector);
            slidingFitResultW.GetGlobalFitPositionAtX(x, fitWVector);

            const float u(fitUVector.GetZ()), v(fitVVector.GetZ()), w(fitWVector.GetZ());
            const float uv2w(LArGeometryHelper::MergeTwoPositions(TPC_VIEW_U, TPC_VIEW_V, u, v));
            const float uw2v(LArGeometryHelper::MergeTwoPositions(TPC_VIEW_U, TPC_VIEW_W, u, w));
            const float vw2u(LArGeometryHelper::MergeTwoPositions(TPC_VIEW_V, TPC_VIEW_W, v, w));

            const unsigned int xBin((x - minX) / xPitch);
            showerPositionMapU.insert(ShowerPositionMap::value_type(xBin, CartesianVector(x, 0., vw2u)));
            showerPositionMapV.insert(ShowerPositionMap::value_type(xBin, CartesianVector(x, 0., uw2v)));
            showerPositionMapW.insert(ShowerPositionMap::value_type(xBin, CartesianVector(x, 0., uv2w)));
        }
        catch (StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDShowersAlgorithm::GetHitOverlapFraction(const Cluster *const pCluster, const float minX, const float maxX, const float xPitch,
    const ShowerPositionMap &edgeMap1, const ShowerPositionMap &edgeMap2, unsigned int &nSampledHits, unsigned int &nMatchedHits) const
{
    if (((maxX - minX) < std::numeric_limits<float>::epsilon()) || (xPitch < std::numeric_limits<float>::epsilon()))
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

            if ((x < minX) || (x > maxX))
                continue;

            ++nSampledHits;
            const unsigned int xBin((x - minX) / xPitch);

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
