/**
 *  @file   LArContent/src/LArThreeDSeed/ThreeDTracksAlgorithm.cc
 * 
 *  @brief  Implementation of the three dimensional tracks algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"

#include "LArThreeDSeed/ThreeDTracksAlgorithm.h"

using namespace pandora;

namespace lar
{

void ThreeDTracksAlgorithm::CalculateOverlapResult(Cluster *pClusterU, Cluster *pClusterV, Cluster *pClusterW)
{
    LArClusterHelper::TwoDSlidingFitResult slidingFitResultU, slidingFitResultV, slidingFitResultW;
    LArClusterHelper::LArTwoDSlidingFit(pClusterU, 20, slidingFitResultU);
    LArClusterHelper::LArTwoDSlidingFit(pClusterV, 20, slidingFitResultV);
    LArClusterHelper::LArTwoDSlidingFit(pClusterW, 20, slidingFitResultW);

    FitSegmentList fitSegmentListU, fitSegmentListV, fitSegmentListW;
    this->GetFitSegmentList(slidingFitResultU, fitSegmentListU);
    this->GetFitSegmentList(slidingFitResultV, fitSegmentListV);
    this->GetFitSegmentList(slidingFitResultW, fitSegmentListW);

    SegmentComparisonList segmentComparisonList;
    this->GetSegmentComparisonList(fitSegmentListU, fitSegmentListV, fitSegmentListW, segmentComparisonList);

    unsigned int nSamplingPoints(0), nMatchedSamplingPoints(0);

    for (SegmentComparisonList::const_iterator iter = segmentComparisonList.begin(), iterEnd = segmentComparisonList.end(); iter != iterEnd; ++iter)
    {
        try
        {
            const TrackOverlapResult trackOverlapResult(this->CalculateOverlapResult(*iter, slidingFitResultU, slidingFitResultV, slidingFitResultW));
            nSamplingPoints += trackOverlapResult.GetNSamplingPoints();
            nMatchedSamplingPoints += trackOverlapResult.GetNMatchedSamplingPoints();
        }
        catch (StatusCodeException &)
        {
            std::cout << "Failed to calculate track overlap result for a list of segments " << std::endl;
        }
    }

    if (0 == nSamplingPoints)
    {
        std::cout << "ThreeDTracksAlgorithm: Cannot calculate overlap result, nSamplingPoints " << nSamplingPoints << std::endl;
        return;
    }

//const float matchedSamplingFraction(static_cast<float>(nMatchedSamplingPoints) / static_cast<float>(nSamplingPoints));
//std::cout << " POPULATE TENSOR: xOverlap " << xOverlap << ", xOverlapU " << (xOverlap / xSpanU) << ", xOverlapV " << (xOverlap / xSpanV) << ", xOverlapW " << (xOverlap / xSpanW) << ", nMatchedSamplingPoints " << nMatchedSamplingPoints << ", nSamplingPoints " << nSamplingPoints << ", matchedSamplingFraction " << matchedSamplingFraction << std::endl;
//PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
//ClusterList clusterListU; clusterListU.insert(pClusterU);
//PandoraMonitoringApi::VisualizeClusters(&clusterListU, "ClusterListU", RED);
//ClusterList clusterListV; clusterListV.insert(pClusterV);
//PandoraMonitoringApi::VisualizeClusters(&clusterListV, "ClusterListV", GREEN);
//ClusterList clusterListW; clusterListW.insert(pClusterW);
//PandoraMonitoringApi::VisualizeClusters(&clusterListW, "ClusterListW", BLUE);
//PandoraMonitoringApi::ViewEvent();
    m_overlapTensor.SetOverlapResult(pClusterU, pClusterV, pClusterW, TrackOverlapResult(nMatchedSamplingPoints, nSamplingPoints));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDTracksAlgorithm::GetFitSegmentList(const LArClusterHelper::TwoDSlidingFitResult &slidingFitResult, FitSegmentList &fitSegmentList) const
{
    // TODO - proper logic needed here
    const LayerFitResultMap &layerFitResultMap(slidingFitResult.GetLayerFitResultMap());

    if (layerFitResultMap.size() < 2)
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    fitSegmentList.push_back(FitSegment(slidingFitResult, layerFitResultMap.begin(), --(layerFitResultMap.end())));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDTracksAlgorithm::GetSegmentComparisonList(const FitSegmentList &fitSegmentListU, const FitSegmentList &fitSegmentListV,
    const FitSegmentList &fitSegmentListW, SegmentComparisonList &segmentComparisonList) const
{
    // TODO - proper logic needed here
    if ((1 != fitSegmentListU.size()) || (1 != fitSegmentListV.size()) || (1 != fitSegmentListW.size()))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    segmentComparisonList.push_back(SegmentComparison(*(fitSegmentListU.begin()), *(fitSegmentListV.begin()), *(fitSegmentListW.begin())));
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackOverlapResult ThreeDTracksAlgorithm::CalculateOverlapResult(const SegmentComparison &segmentComparison,
    const LArClusterHelper::TwoDSlidingFitResult &slidingFitResultU, const LArClusterHelper::TwoDSlidingFitResult &slidingFitResultV,
    const LArClusterHelper::TwoDSlidingFitResult &slidingFitResultW) const
{
    const float minXU(segmentComparison.GetFitSegmentU().GetMinX());
    const float maxXU(segmentComparison.GetFitSegmentU().GetMaxX());

    const float minXV(segmentComparison.GetFitSegmentV().GetMinX());
    const float maxXV(segmentComparison.GetFitSegmentV().GetMaxX());

    const float minXW(segmentComparison.GetFitSegmentW().GetMinX());
    const float maxXW(segmentComparison.GetFitSegmentW().GetMaxX());

    const float xSpanU(maxXU - minXU);
    const float xSpanV(maxXV - minXV);
    const float xSpanW(maxXW - minXW);

    // Assess x-overlap
    const float minX(std::max(minXU, std::max(minXV, minXW)));
    const float maxX(std::min(maxXU, std::min(maxXV, maxXW)));
    const float xOverlap(maxX - minX);

    if ((xOverlap < 0.f) || ((xOverlap / xSpanU) < 0.3f) || ((xOverlap / xSpanV) < 0.3f) || ((xOverlap / xSpanW) < 0.3f))
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    // Sampling in x, TODO, work out x pitch for segment comparison
    const float nPointsU((xOverlap / xSpanU) * slidingFitResultU.GetCluster()->GetNCaloHits());
    const float nPointsV((xOverlap / xSpanV) * slidingFitResultV.GetCluster()->GetNCaloHits());
    const float nPointsW((xOverlap / xSpanW) * slidingFitResultW.GetCluster()->GetNCaloHits());

    const float xPitch(3.f * xOverlap / (nPointsU + nPointsV + nPointsW));

    // Chi2 calculations
    unsigned int nSamplingPoints(0), nMatchedSamplingPoints(0);

    for (float x = minX; x < maxX; x += xPitch)
    {
        try
        {
            CartesianVector fitUVector(0.f, 0.f, 0.f), fitVVector(0.f, 0.f, 0.f), fitWVector(0.f, 0.f, 0.f);
            slidingFitResultU.GetGlobalFitPosition(x, true, fitUVector);
            slidingFitResultV.GetGlobalFitPosition(x, true, fitVVector);
            slidingFitResultW.GetGlobalFitPosition(x, true, fitWVector);

            const float u(fitUVector.GetZ()), v(fitVVector.GetZ()), w(fitWVector.GetZ());
            const float uv2w(LArGeometryHelper::MergeTwoPositions(VIEW_U, VIEW_V, u, v));
            const float uw2v(LArGeometryHelper::MergeTwoPositions(VIEW_U, VIEW_W, u, w));
            const float vw2u(LArGeometryHelper::MergeTwoPositions(VIEW_V, VIEW_W, v, w));

            ++nSamplingPoints;
            const float deltaW(uv2w - w), deltaV(uw2v - v), deltaU(vw2u - u);
            const float pseudoChi2(deltaW * deltaW + deltaV * deltaV + deltaU * deltaU);

            if (pseudoChi2 < m_pseudoChi2Cut)
                ++nMatchedSamplingPoints;
//const CartesianVector expU(x, 0., vw2u); PANDORA_MONITORING_API(AddMarkerToVisualization(&expU, "expU", RED, 1.));
//const CartesianVector expV(x, 0., uw2v); PANDORA_MONITORING_API(AddMarkerToVisualization(&expV, "expV", GREEN, 1.));
//const CartesianVector expW(x, 0., uv2w); PANDORA_MONITORING_API(AddMarkerToVisualization(&expW, "expW", BLUE, 1.));
//std::cout << " TRK pseudoChi2 " << pseudoChi2 << " deltaW " << deltaW << " deltaV " << deltaV << " deltaU " << deltaU << std::endl;
        }
        catch (StatusCodeException &)
        {
        }
    }

    if (0 == nSamplingPoints)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return TrackOverlapResult(nMatchedSamplingPoints, nSamplingPoints);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDTracksAlgorithm::CalculateConstantXOverlapResult(const LArClusterHelper::TwoDSlidingFitResult &slidingFitResultU,
    const LArClusterHelper::TwoDSlidingFitResult &slidingFitResultV, const LArClusterHelper::TwoDSlidingFitResult &slidingFitResultW)
{
    // TODO - Blake code to go here... fasten your seatbelt!!!
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ThreeDTracksAlgorithm::ExamineTensor()
{
    float bestOverlapResult(m_minMatchedFraction);
    Cluster *pBestClusterU(NULL), *pBestClusterV(NULL), *pBestClusterW(NULL);

    const ClusterList &clusterListU(m_overlapTensor.GetClusterListU());
    const ClusterList &clusterListV(m_overlapTensor.GetClusterListV());
    const ClusterList &clusterListW(m_overlapTensor.GetClusterListW());

    for (ClusterList::const_iterator iterU = clusterListU.begin(), iterUEnd = clusterListU.end(); iterU != iterUEnd; ++iterU)
    {
        for (ClusterList::const_iterator iterV = clusterListV.begin(), iterVEnd = clusterListV.end(); iterV != iterVEnd; ++iterV)
        {
            for (ClusterList::const_iterator iterW = clusterListW.begin(), iterWEnd = clusterListW.end(); iterW != iterWEnd; ++iterW)
            {
                try
                {
                    const TrackOverlapResult &overlapResult(m_overlapTensor.GetOverlapResult(*iterU, *iterV, *iterW));

                    if (overlapResult.GetNMatchedSamplingPoints() < m_minMatchedPoints)
                        continue;

                    if (overlapResult.GetMatchedFraction() > bestOverlapResult)
                    {
                        bestOverlapResult = overlapResult.GetMatchedFraction();
                        pBestClusterU = *iterU;
                        pBestClusterV = *iterV;
                        pBestClusterW = *iterW;
                    }
                }
                catch (StatusCodeException &)
                {
                }
            }
        }
    }

    if (!pBestClusterU || !pBestClusterV || !pBestClusterW)
        return false;

    ProtoParticle protoParticle;
    protoParticle.m_clusterVectorU.push_back(pBestClusterU);
    protoParticle.m_clusterVectorV.push_back(pBestClusterV);
    protoParticle.m_clusterVectorW.push_back(pBestClusterW);
    m_protoParticleVector.push_back(protoParticle);

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDTracksAlgorithm::FitSegment::FitSegment(const LArClusterHelper::TwoDSlidingFitResult &twoDSlidingFitResult,
        LayerFitResultMap::const_iterator startLayerIter, LayerFitResultMap::const_iterator endLayerIter) :
    m_startLayer(startLayerIter->first),
    m_endLayer(endLayerIter->first)
{
    CartesianVector startLayerPosition(0.f, 0.f, 0.f);
    twoDSlidingFitResult.GetGlobalPosition(startLayerIter->second.GetL(), startLayerIter->second.GetFitT(), startLayerPosition);

    CartesianVector endLayerPosition(0.f, 0.f, 0.f);
    twoDSlidingFitResult.GetGlobalPosition(endLayerIter->second.GetL(), endLayerIter->second.GetFitT(), endLayerPosition);

    m_minX = std::min(startLayerPosition.GetX(), endLayerPosition.GetX());
    m_maxX = std::max(startLayerPosition.GetX(), endLayerPosition.GetX());
    m_startValue = startLayerPosition.GetZ();
    m_endValue = endLayerPosition.GetZ();
    m_isIncreasingX = (endLayerPosition.GetX() > startLayerPosition.GetX());
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDTracksAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_pseudoChi2Cut = 3.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PseudoChi2Cut", m_pseudoChi2Cut));

    m_minMatchedFraction = 0.8f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedFraction", m_minMatchedFraction));

    m_minMatchedPoints = 0;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedPoints", m_minMatchedPoints));

    return ThreeDBaseAlgorithm<TrackOverlapResult>::ReadSettings(xmlHandle);
}

} // namespace lar
