/**
 *  @file   LArContent/src/LArThreeDSeed/ThreeDShowersAlgorithm.cc
 * 
 *  @brief  Implementation of the three dimensional showers algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"

#include "LArThreeDSeed/ThreeDShowersAlgorithm.h"

using namespace pandora;

namespace lar
{

void ThreeDShowersAlgorithm::CalculateOverlapResult(Cluster *pClusterU, Cluster *pClusterV, Cluster *pClusterW)
{
    static const unsigned int m_layerFitHalfWindow = 100;

    // U
    LArClusterHelper::TwoDSlidingFitResult slidingFitResultU;
    LArClusterHelper::LArTwoDSlidingFit(pClusterU, m_layerFitHalfWindow, slidingFitResultU);
    LArClusterHelper::TwoDSlidingFitResult positiveShowerEdgeFitU;
    LArClusterHelper::LArTwoDShowerEdgeFit(pClusterU, m_layerFitHalfWindow, slidingFitResultU.GetAxisIntercept(), slidingFitResultU.GetAxisDirection(), POSITIVE_SHOWER_EDGE, positiveShowerEdgeFitU);
    LArClusterHelper::TwoDSlidingFitResult negativeShowerEdgeFitU;
    LArClusterHelper::LArTwoDShowerEdgeFit(pClusterU, m_layerFitHalfWindow, slidingFitResultU.GetAxisIntercept(), slidingFitResultU.GetAxisDirection(), NEGATIVE_SHOWER_EDGE, negativeShowerEdgeFitU);

    const float innerXU(slidingFitResultU.GetGlobalMinLayerPosition().GetX());
    const float outerXU(slidingFitResultU.GetGlobalMaxLayerPosition().GetX());
    const float minXU(std::min(innerXU, outerXU));
    const float maxXU(std::max(innerXU, outerXU));
    const float xSpanU(maxXU - minXU);

    // V
    LArClusterHelper::TwoDSlidingFitResult slidingFitResultV;
    LArClusterHelper::LArTwoDSlidingFit(pClusterV, m_layerFitHalfWindow, slidingFitResultV);
    LArClusterHelper::TwoDSlidingFitResult positiveShowerEdgeFitV;
    LArClusterHelper::LArTwoDShowerEdgeFit(pClusterV, m_layerFitHalfWindow, slidingFitResultV.GetAxisIntercept(), slidingFitResultV.GetAxisDirection(), POSITIVE_SHOWER_EDGE, positiveShowerEdgeFitV);
    LArClusterHelper::TwoDSlidingFitResult negativeShowerEdgeFitV;
    LArClusterHelper::LArTwoDShowerEdgeFit(pClusterV, m_layerFitHalfWindow, slidingFitResultV.GetAxisIntercept(), slidingFitResultV.GetAxisDirection(), NEGATIVE_SHOWER_EDGE, negativeShowerEdgeFitV);

    const float innerXV(slidingFitResultV.GetGlobalMinLayerPosition().GetX());
    const float outerXV(slidingFitResultV.GetGlobalMaxLayerPosition().GetX());
    const float minXV(std::min(innerXV, outerXV));
    const float maxXV(std::max(innerXV, outerXV));
    const float xSpanV(maxXV - minXV);

    // W
    LArClusterHelper::TwoDSlidingFitResult slidingFitResultW;
    LArClusterHelper::LArTwoDSlidingFit(pClusterW, m_layerFitHalfWindow, slidingFitResultW);
    LArClusterHelper::TwoDSlidingFitResult positiveShowerEdgeFitW;
    LArClusterHelper::LArTwoDShowerEdgeFit(pClusterW, m_layerFitHalfWindow, slidingFitResultW.GetAxisIntercept(), slidingFitResultW.GetAxisDirection(), POSITIVE_SHOWER_EDGE, positiveShowerEdgeFitW);
    LArClusterHelper::TwoDSlidingFitResult negativeShowerEdgeFitW;
    LArClusterHelper::LArTwoDShowerEdgeFit(pClusterW, m_layerFitHalfWindow, slidingFitResultW.GetAxisIntercept(), slidingFitResultW.GetAxisDirection(), NEGATIVE_SHOWER_EDGE, negativeShowerEdgeFitW);

    const float innerXW(slidingFitResultW.GetGlobalMinLayerPosition().GetX());
    const float outerXW(slidingFitResultW.GetGlobalMaxLayerPosition().GetX());
    const float minXW(std::min(innerXW, outerXW));
    const float maxXW(std::max(innerXW, outerXW));
    const float xSpanW(maxXW - minXW);

    // Assess x-overlap
    const float minX(std::max(minXU, std::max(minXV, minXW)));
    const float maxX(std::min(maxXU, std::min(maxXV, maxXW)));
    const float xOverlap(maxX - minX);

    if ((xOverlap < 0.f) || ((xOverlap / xSpanU) < 0.3f) || ((xOverlap / xSpanV) < 0.3f) || ((xOverlap / xSpanW) < 0.3f))
    {
std::cout << "Fail overlap criteria, xOverlap " << xOverlap << " xSpanU " << (xOverlap / xSpanU) << " xSpanV " << (xOverlap / xSpanV) << " xSpanW " << (xOverlap / xSpanW) << std::endl;
//PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
//ClusterList clusterListU; clusterListU.insert(pClusterU);
//PandoraMonitoringApi::VisualizeClusters(&clusterListU, "pClusterU", RED);
//ClusterList clusterListV; clusterListV.insert(pClusterV);
//PandoraMonitoringApi::VisualizeClusters(&clusterListV, "pClusterV", GREEN);
//ClusterList clusterListW; clusterListW.insert(pClusterW);
//PandoraMonitoringApi::VisualizeClusters(&clusterListW, "pClusterW", BLUE);
//PandoraMonitoringApi::ViewEvent();
        return;
    }
/*
// DEBUG
//------------------------------------------------------------------------------------------------------------------------------------------
std::cout << " ClusterU and sliding shower fits " << std::endl;
const LArClusterHelper::TwoDSlidingFitResult::LayerFitResultMap layerFitResultMapU1(slidingFitResultU.GetLayerFitResultMap());
for (LArClusterHelper::TwoDSlidingFitResult::LayerFitResultMap::const_iterator iter = layerFitResultMapU1.begin(); iter != layerFitResultMapU1.end(); ++iter)
{
const int layer(iter->first);const float rL(iter->second.GetL());const float rT(iter->second.GetFitT());
CartesianVector position(0.f, 0.f, 0.f);slidingFitResultU.GetGlobalCoordinates(rL, rT, position);
PandoraMonitoringApi::AddMarkerToVisualization(&position, "posU", CYAN, 1.);
}
const LArClusterHelper::TwoDSlidingFitResult::LayerFitResultMap layerFitResultMapU2(positiveShowerEdgeFitU.GetLayerFitResultMap());
for (LArClusterHelper::TwoDSlidingFitResult::LayerFitResultMap::const_iterator iter = layerFitResultMapU2.begin(); iter != layerFitResultMapU2.end(); ++iter)
{
const int layer(iter->first);const float rL(iter->second.GetL());const float rT(iter->second.GetFitT());
CartesianVector position(0.f, 0.f, 0.f);positiveShowerEdgeFitU.GetGlobalCoordinates(rL, rT, position);
PandoraMonitoringApi::AddMarkerToVisualization(&position, "posU", GRAY, 1.);
}
const LArClusterHelper::TwoDSlidingFitResult::LayerFitResultMap layerFitResultMapU3(negativeShowerEdgeFitU.GetLayerFitResultMap());
for (LArClusterHelper::TwoDSlidingFitResult::LayerFitResultMap::const_iterator iter = layerFitResultMapU3.begin(); iter != layerFitResultMapU3.end(); ++iter)
{
const int layer(iter->first);const float rL(iter->second.GetL());const float rT(iter->second.GetFitT());
CartesianVector position(0.f, 0.f, 0.f);negativeShowerEdgeFitU.GetGlobalCoordinates(rL, rT, position);
PandoraMonitoringApi::AddMarkerToVisualization(&position, "posU", GRAY, 1.);
}
PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
ClusterList clusterListU1; clusterListU1.insert(pClusterU);
PandoraMonitoringApi::VisualizeClusters(&clusterListU1, "ClusterListU", RED);
PandoraMonitoringApi::ViewEvent();
//------------------------------------------------------------------------------------------------------------------------------------------
std::cout << " ClusterV and sliding shower fits " << std::endl;
const LArClusterHelper::TwoDSlidingFitResult::LayerFitResultMap layerFitResultMapV1(slidingFitResultV.GetLayerFitResultMap());
for (LArClusterHelper::TwoDSlidingFitResult::LayerFitResultMap::const_iterator iter = layerFitResultMapV1.begin(); iter != layerFitResultMapV1.end(); ++iter)
{
const int layer(iter->first);const float rL(iter->second.GetL());const float rT(iter->second.GetFitT());
CartesianVector position(0.f, 0.f, 0.f);slidingFitResultV.GetGlobalCoordinates(rL, rT, position);
PandoraMonitoringApi::AddMarkerToVisualization(&position, "posV", CYAN, 1.);
}
const LArClusterHelper::TwoDSlidingFitResult::LayerFitResultMap layerFitResultMapV2(positiveShowerEdgeFitV.GetLayerFitResultMap());
for (LArClusterHelper::TwoDSlidingFitResult::LayerFitResultMap::const_iterator iter = layerFitResultMapV2.begin(); iter != layerFitResultMapV2.end(); ++iter)
{
const int layer(iter->first);const float rL(iter->second.GetL());const float rT(iter->second.GetFitT());
CartesianVector position(0.f, 0.f, 0.f);positiveShowerEdgeFitV.GetGlobalCoordinates(rL, rT, position);
PandoraMonitoringApi::AddMarkerToVisualization(&position, "posV", GRAY, 1.);
}
const LArClusterHelper::TwoDSlidingFitResult::LayerFitResultMap layerFitResultMapV3(negativeShowerEdgeFitV.GetLayerFitResultMap());
for (LArClusterHelper::TwoDSlidingFitResult::LayerFitResultMap::const_iterator iter = layerFitResultMapV3.begin(); iter != layerFitResultMapV3.end(); ++iter)
{
const int layer(iter->first);const float rL(iter->second.GetL());const float rT(iter->second.GetFitT());
CartesianVector position(0.f, 0.f, 0.f);negativeShowerEdgeFitV.GetGlobalCoordinates(rL, rT, position);
PandoraMonitoringApi::AddMarkerToVisualization(&position, "posV", GRAY, 1.);
}
PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
ClusterList clusterListV1; clusterListV1.insert(pClusterV);
PandoraMonitoringApi::VisualizeClusters(&clusterListV1, "ClusterListV", RED);
PandoraMonitoringApi::ViewEvent();
//------------------------------------------------------------------------------------------------------------------------------------------
std::cout << " ClusterW and sliding shower fits " << std::endl;
const LArClusterHelper::TwoDSlidingFitResult::LayerFitResultMap layerFitResultMapW1(slidingFitResultW.GetLayerFitResultMap());
for (LArClusterHelper::TwoDSlidingFitResult::LayerFitResultMap::const_iterator iter = layerFitResultMapW1.begin(); iter != layerFitResultMapW1.end(); ++iter)
{
const int layer(iter->first);const float rL(iter->second.GetL());const float rT(iter->second.GetFitT());
CartesianVector position(0.f, 0.f, 0.f);slidingFitResultW.GetGlobalCoordinates(rL, rT, position);
PandoraMonitoringApi::AddMarkerToVisualization(&position, "posW", CYAN, 1.);
}
const LArClusterHelper::TwoDSlidingFitResult::LayerFitResultMap layerFitResultMapW2(positiveShowerEdgeFitW.GetLayerFitResultMap());
for (LArClusterHelper::TwoDSlidingFitResult::LayerFitResultMap::const_iterator iter = layerFitResultMapW2.begin(); iter != layerFitResultMapW2.end(); ++iter)
{
const int layer(iter->first);const float rL(iter->second.GetL());const float rT(iter->second.GetFitT());
CartesianVector position(0.f, 0.f, 0.f);positiveShowerEdgeFitW.GetGlobalCoordinates(rL, rT, position);
PandoraMonitoringApi::AddMarkerToVisualization(&position, "posW", GRAY, 1.);
}
const LArClusterHelper::TwoDSlidingFitResult::LayerFitResultMap layerFitResultMapW3(negativeShowerEdgeFitW.GetLayerFitResultMap());
for (LArClusterHelper::TwoDSlidingFitResult::LayerFitResultMap::const_iterator iter = layerFitResultMapW3.begin(); iter != layerFitResultMapW3.end(); ++iter)
{
const int layer(iter->first);const float rL(iter->second.GetL());const float rT(iter->second.GetFitT());
CartesianVector position(0.f, 0.f, 0.f);negativeShowerEdgeFitW.GetGlobalCoordinates(rL, rT, position);
PandoraMonitoringApi::AddMarkerToVisualization(&position, "posW", GRAY, 1.);
}
PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
ClusterList clusterListW1; clusterListW1.insert(pClusterW);
PandoraMonitoringApi::VisualizeClusters(&clusterListW1, "ClusterListW", RED);
PandoraMonitoringApi::ViewEvent();
//------------------------------------------------------------------------------------------------------------------------------------------
*/
    ShowerEdgeMap uMap, uPosMap, uNegMap;
    ShowerEdgeMap vMap, vPosMap, vNegMap;
    ShowerEdgeMap wMap, wPosMap, wNegMap;

    // Sampling in x
    const float nPointsU((xOverlap / xSpanU) * pClusterU->GetNCaloHits());
    const float nPointsV((xOverlap / xSpanV) * pClusterV->GetNCaloHits());
    const float nPointsW((xOverlap / xSpanW) * pClusterW->GetNCaloHits());
    const float xPitch(3.f * xOverlap / (nPointsU + nPointsV + nPointsW));

    for (float x = minX; x < maxX; x += xPitch)
    {
        // Fit to all hits in shower
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

            const unsigned int xBin((x - minX) / xPitch);
            uMap.insert(ShowerEdgeMap::value_type(xBin, CartesianVector(x, 0., vw2u)));
            vMap.insert(ShowerEdgeMap::value_type(xBin, CartesianVector(x, 0., uw2v)));
            wMap.insert(ShowerEdgeMap::value_type(xBin, CartesianVector(x, 0., uv2w)));
        }
        catch (StatusCodeException &)
        {
        }

        // Positive shower edge
        try
        {
            CartesianVector fitUVector(0.f, 0.f, 0.f), fitVVector(0.f, 0.f, 0.f), fitWVector(0.f, 0.f, 0.f);
            positiveShowerEdgeFitU.GetGlobalFitPosition(x, true, fitUVector);
            positiveShowerEdgeFitV.GetGlobalFitPosition(x, true, fitVVector);
            positiveShowerEdgeFitW.GetGlobalFitPosition(x, true, fitWVector);

            const float u(fitUVector.GetZ()), v(fitVVector.GetZ()), w(fitWVector.GetZ());
            const float uv2w(LArGeometryHelper::MergeTwoPositions(VIEW_U, VIEW_V, u, v));
            const float uw2v(LArGeometryHelper::MergeTwoPositions(VIEW_U, VIEW_W, u, w));
            const float vw2u(LArGeometryHelper::MergeTwoPositions(VIEW_V, VIEW_W, v, w));

            const unsigned int xBin((x - minX) / xPitch);
            uPosMap.insert(ShowerEdgeMap::value_type(xBin, CartesianVector(x, 0., vw2u)));
            vPosMap.insert(ShowerEdgeMap::value_type(xBin, CartesianVector(x, 0., uw2v)));
            wPosMap.insert(ShowerEdgeMap::value_type(xBin, CartesianVector(x, 0., uv2w)));
        }
        catch (StatusCodeException &)
        {
        }

        // Negative shower edge
        try
        {
            CartesianVector fitUVector(0.f, 0.f, 0.f), fitVVector(0.f, 0.f, 0.f), fitWVector(0.f, 0.f, 0.f);
            negativeShowerEdgeFitU.GetGlobalFitPosition(x, true, fitUVector);
            negativeShowerEdgeFitV.GetGlobalFitPosition(x, true, fitVVector);
            negativeShowerEdgeFitW.GetGlobalFitPosition(x, true, fitWVector);

            const float u(fitUVector.GetZ()), v(fitVVector.GetZ()), w(fitWVector.GetZ());
            const float uv2w(LArGeometryHelper::MergeTwoPositions(VIEW_U, VIEW_V, u, v));
            const float uw2v(LArGeometryHelper::MergeTwoPositions(VIEW_U, VIEW_W, u, w));
            const float vw2u(LArGeometryHelper::MergeTwoPositions(VIEW_V, VIEW_W, v, w));

            const unsigned int xBin((x - minX) / xPitch);
            uNegMap.insert(ShowerEdgeMap::value_type(xBin, CartesianVector(x, 0., vw2u)));
            vNegMap.insert(ShowerEdgeMap::value_type(xBin, CartesianVector(x, 0., uw2v)));
            wNegMap.insert(ShowerEdgeMap::value_type(xBin, CartesianVector(x, 0., uv2w)));
        }
        catch (StatusCodeException &)
        {
        }
    }

//------------------------------------------------------------------------------------------------------------------------------------------
    const float includedFractionU(this->GetIncludedHitFraction(pClusterU, minX, maxX, xPitch, uPosMap, uNegMap));
//    std::cout << "ClusterU: predicted shower shape, includedFractionU " << includedFractionU << std::endl;
//for (ShowerEdgeMap::const_iterator iter = uMap.begin(); iter != uMap.end(); ++iter) PandoraMonitoringApi::AddMarkerToVisualization(&(iter->second), "ctrU", CYAN, 1.);
//for (ShowerEdgeMap::const_iterator iter = uPosMap.begin(); iter != uPosMap.end(); ++iter) PandoraMonitoringApi::AddMarkerToVisualization(&(iter->second), "posU", GRAY, 1.);
//for (ShowerEdgeMap::const_iterator iter = uNegMap.begin(); iter != uNegMap.end(); ++iter) PandoraMonitoringApi::AddMarkerToVisualization(&(iter->second), "negU", GRAY, 1.);
//PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
//ClusterList clusterListU2; clusterListU2.insert(pClusterU);
//PandoraMonitoringApi::VisualizeClusters(&clusterListU2, "ClusterListU", RED);
//PandoraMonitoringApi::ViewEvent();

    const float includedFractionV(this->GetIncludedHitFraction(pClusterV, minX, maxX, xPitch, vPosMap, vNegMap));
//    std::cout << "ClusterV: predicted shower shape, includedFractionV " << includedFractionV << std::endl;
//for (ShowerEdgeMap::const_iterator iter = vMap.begin(); iter != vMap.end(); ++iter) PandoraMonitoringApi::AddMarkerToVisualization(&(iter->second), "ctrV", CYAN, 1.);
//for (ShowerEdgeMap::const_iterator iter = vPosMap.begin(); iter != vPosMap.end(); ++iter) PandoraMonitoringApi::AddMarkerToVisualization(&(iter->second), "posV", GRAY, 1.);
//for (ShowerEdgeMap::const_iterator iter = vNegMap.begin(); iter != vNegMap.end(); ++iter) PandoraMonitoringApi::AddMarkerToVisualization(&(iter->second), "negV", GRAY, 1.);
//PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
//ClusterList clusterListV2; clusterListV2.insert(pClusterV);
//PandoraMonitoringApi::VisualizeClusters(&clusterListV2, "ClusterListV", RED);
//PandoraMonitoringApi::ViewEvent();

    const float includedFractionW(this->GetIncludedHitFraction(pClusterW, minX, maxX, xPitch, wPosMap, wNegMap));
//    std::cout << "ClusterW: predicted shower shape, includedFractionW " << includedFractionW << std::endl;
//for (ShowerEdgeMap::const_iterator iter = wMap.begin(); iter != wMap.end(); ++iter) PandoraMonitoringApi::AddMarkerToVisualization(&(iter->second), "ctrW", CYAN, 1.);
//for (ShowerEdgeMap::const_iterator iter = wPosMap.begin(); iter != wPosMap.end(); ++iter) PandoraMonitoringApi::AddMarkerToVisualization(&(iter->second), "posW", GRAY, 1.);
//for (ShowerEdgeMap::const_iterator iter = wNegMap.begin(); iter != wNegMap.end(); ++iter) PandoraMonitoringApi::AddMarkerToVisualization(&(iter->second), "negW", GRAY, 1.);
//PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
//ClusterList clusterListW2; clusterListW2.insert(pClusterW);
//PandoraMonitoringApi::VisualizeClusters(&clusterListW2, "ClusterListW", RED);
//PandoraMonitoringApi::ViewEvent();
//------------------------------------------------------------------------------------------------------------------------------------------

    if ((includedFractionU < 0.5f) || (includedFractionV < 0.5f) || (includedFractionW < 0.5f))
        return;

    m_overlapTensor.SetOverlapResult(pClusterU, pClusterV, pClusterW, (includedFractionU + includedFractionV + includedFractionW) / 3.f);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ThreeDShowersAlgorithm::ExamineTensor()
{
    float bestOverlapResult(0.5f); // TODO Min overlap result for PFO creation
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
                    const float overlapResult(m_overlapTensor.GetOverlapResult(*iterU, *iterV, *iterW));

                    if (overlapResult > bestOverlapResult)
                    {
                        bestOverlapResult = overlapResult;
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
    protoParticle.m_clusterListU.insert(pBestClusterU);
    protoParticle.m_clusterListV.insert(pBestClusterV);
    protoParticle.m_clusterListW.insert(pBestClusterW);
    m_protoParticleVector.push_back(protoParticle);

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ThreeDShowersAlgorithm::GetIncludedHitFraction(const Cluster *const pCluster, const float minX, const float maxX, const float xPitch,
    const ShowerEdgeMap &edgeMap1, const ShowerEdgeMap &edgeMap2) const
{
    if (((maxX - minX) < std::numeric_limits<float>::epsilon()) || (xPitch < std::numeric_limits<float>::epsilon()))
        {std::cout << "maxX " << maxX << " minX " << minX << " xPitch " << xPitch << std::endl; throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);}

    unsigned int nHitsInRange(0), nIncludedHits(0);
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
//CaloHitList included;
    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hIter = iter->second->begin(), hIterEnd = iter->second->end(); hIter != hIterEnd; ++hIter)
        {
            CaloHit *pCaloHit = *hIter;
            const float x(pCaloHit->GetPositionVector().GetX());
            const float z(pCaloHit->GetPositionVector().GetZ());

            if ((x < minX) || (x > maxX))
                continue;

            ++nHitsInRange;
            const unsigned int xBin((x - minX) / xPitch);

            ShowerEdgeMap::const_iterator edgeIter1 = edgeMap1.find(xBin);
            ShowerEdgeMap::const_iterator edgeIter2 = edgeMap2.find(xBin);

            if ((edgeMap1.end() == edgeIter1) || (edgeMap2.end() == edgeIter2))
                continue;

            const float maxZ(std::max(edgeIter1->second.GetZ(), edgeIter2->second.GetZ()));
            const float minZ(std::min(edgeIter1->second.GetZ(), edgeIter2->second.GetZ()));

            if ((z > minZ) && (z < maxZ))
                {/*included.insert(pCaloHit);*/ ++nIncludedHits;}
        }
    }

//PandoraMonitoringApi::VisualizeCaloHits(&included, "included", BLUE);

    if (0 == nHitsInRange)
        return 0.f;

    return (static_cast<float>(nIncludedHits) / static_cast<float>(nHitsInRange));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDShowersAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    return ThreeDBaseAlgorithm<float>::ReadSettings(xmlHandle);
}

} // namespace lar
