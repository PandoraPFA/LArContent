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

void ThreeDShowersAlgorithm::SelectInputClusters()
{
    this->SelectInputClusters(m_pInputClusterListU, m_clusterVectorU);
    this->SelectInputClusters(m_pInputClusterListV, m_clusterVectorV);
    this->SelectInputClusters(m_pInputClusterListW, m_clusterVectorW);

    std::sort(m_clusterVectorU.begin(), m_clusterVectorU.end(), LArClusterHelper::SortByNOccupiedLayers);
    std::sort(m_clusterVectorV.begin(), m_clusterVectorV.end(), LArClusterHelper::SortByNOccupiedLayers);
    std::sort(m_clusterVectorW.begin(), m_clusterVectorW.end(), LArClusterHelper::SortByNOccupiedLayers);

std::cout << "Clusters for 2D->3D matching " << std::endl;
PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
ClusterList clusterListU; clusterListU.insert(m_clusterVectorU.begin(), m_clusterVectorU.end());
PandoraMonitoringApi::VisualizeClusters(&clusterListU, "ClusterListU", RED);
ClusterList clusterListV; clusterListV.insert(m_clusterVectorV.begin(), m_clusterVectorV.end());
PandoraMonitoringApi::VisualizeClusters(&clusterListV, "ClusterListV", GREEN);
ClusterList clusterListW; clusterListW.insert(m_clusterVectorW.begin(), m_clusterVectorW.end());
PandoraMonitoringApi::VisualizeClusters(&clusterListW, "ClusterListW", BLUE);
PandoraMonitoringApi::ViewEvent();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDShowersAlgorithm::SelectInputClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (!pCluster->IsAvailable())
            continue;

        clusterVector.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDShowersAlgorithm::ModifyInputClusters()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDShowersAlgorithm::InitializeTensor()
{
    m_overlapTensor.Clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

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
        return;

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

typedef std::vector<CartesianVector> CartesianVectorList;
CartesianVectorList uList, uPosList, uNegList;
CartesianVectorList vList, vPosList, vNegList;
CartesianVectorList wList, wPosList, wNegList;
//------------------------------------------------------------------------------------------------------------------------------------------

    // TODO
    // if (slidingFitResultU.IsMultivaluedInX() || slidingFitResultV.IsMultivaluedInX() || slidingFitResultW.IsMultivaluedInX())
    //     return this->CalculateConstantXOverlapResult(slidingFitResultU, slidingFitResultV, slidingFitResultW);

    // Sampling in x
    const float nPointsU((xOverlap / xSpanU) * pClusterU->GetNCaloHits());
    const float nPointsV((xOverlap / xSpanV) * pClusterV->GetNCaloHits());
    const float nPointsW((xOverlap / xSpanW) * pClusterW->GetNCaloHits());
    const float xPitch(3.f * xOverlap / (nPointsU + nPointsV + nPointsW));

    // Chi2 calculations
    unsigned int nSamplingPoints(0), nMatchedSamplingPoints(0);
    unsigned int nPosSamplingPoints(0), nPosMatchedSamplingPoints(0);
    unsigned int nNegSamplingPoints(0), nNegMatchedSamplingPoints(0);

    for (float x = minX; x < maxX; x += xPitch)
    {
        // Fit to all hits in shower
        try
        {
            CartesianVector fitUVector(0.f, 0.f, 0.f), fitVVector(0.f, 0.f, 0.f), fitWVector(0.f, 0.f, 0.f);
            slidingFitResultU.GetGlobalFitCoordinates(x, fitUVector);
            slidingFitResultV.GetGlobalFitCoordinates(x, fitVVector);
            slidingFitResultW.GetGlobalFitCoordinates(x, fitWVector);

            const float u(fitUVector.GetZ()), v(fitVVector.GetZ()), w(fitWVector.GetZ());
            const float uv2w(LArGeometryHelper::MergeTwoPositions(VIEW_U, VIEW_V, u, v));
            const float uw2v(LArGeometryHelper::MergeTwoPositions(VIEW_U, VIEW_W, u, w));
            const float vw2u(LArGeometryHelper::MergeTwoPositions(VIEW_V, VIEW_W, v, w));

            ++nSamplingPoints;
            const float deltaW(uv2w - w), deltaV(uw2v - v), deltaU(vw2u - u);
            const float pseudoChi2(deltaW * deltaW + deltaV * deltaV + deltaU * deltaU);

            if (pseudoChi2 < 6.f)
                ++nMatchedSamplingPoints;
uList.push_back(CartesianVector(x, 0., vw2u));
vList.push_back(CartesianVector(x, 0., uw2v));
wList.push_back(CartesianVector(x, 0., uv2w));
        }
        catch (StatusCodeException &)
        {
        }

        // Positive shower edge
        try
        {
            CartesianVector fitUVector(0.f, 0.f, 0.f), fitVVector(0.f, 0.f, 0.f), fitWVector(0.f, 0.f, 0.f);
            positiveShowerEdgeFitU.GetGlobalFitCoordinates(x, fitUVector);
            positiveShowerEdgeFitV.GetGlobalFitCoordinates(x, fitVVector);
            positiveShowerEdgeFitW.GetGlobalFitCoordinates(x, fitWVector);

            const float u(fitUVector.GetZ()), v(fitVVector.GetZ()), w(fitWVector.GetZ());
            const float uv2w(LArGeometryHelper::MergeTwoPositions(VIEW_U, VIEW_V, u, v));
            const float uw2v(LArGeometryHelper::MergeTwoPositions(VIEW_U, VIEW_W, u, w));
            const float vw2u(LArGeometryHelper::MergeTwoPositions(VIEW_V, VIEW_W, v, w));

            ++nPosSamplingPoints;
            const float deltaW(uv2w - w), deltaV(uw2v - v), deltaU(vw2u - u);
            const float pseudoChi2(deltaW * deltaW + deltaV * deltaV + deltaU * deltaU);

            if (pseudoChi2 < 6.f)
                ++nPosMatchedSamplingPoints;
uPosList.push_back(CartesianVector(x, 0., vw2u));
vPosList.push_back(CartesianVector(x, 0., uw2v));
wPosList.push_back(CartesianVector(x, 0., uv2w));
        }
        catch (StatusCodeException &)
        {
        }

        // Negative shower edge
        try
        {
            CartesianVector fitUVector(0.f, 0.f, 0.f), fitVVector(0.f, 0.f, 0.f), fitWVector(0.f, 0.f, 0.f);
            negativeShowerEdgeFitU.GetGlobalFitCoordinates(x, fitUVector);
            negativeShowerEdgeFitV.GetGlobalFitCoordinates(x, fitVVector);
            negativeShowerEdgeFitW.GetGlobalFitCoordinates(x, fitWVector);

            const float u(fitUVector.GetZ()), v(fitVVector.GetZ()), w(fitWVector.GetZ());
            const float uv2w(LArGeometryHelper::MergeTwoPositions(VIEW_U, VIEW_V, u, v));
            const float uw2v(LArGeometryHelper::MergeTwoPositions(VIEW_U, VIEW_W, u, w));
            const float vw2u(LArGeometryHelper::MergeTwoPositions(VIEW_V, VIEW_W, v, w));

            ++nNegSamplingPoints;
            const float deltaW(uv2w - w), deltaV(uw2v - v), deltaU(vw2u - u);
            const float pseudoChi2(deltaW * deltaW + deltaV * deltaV + deltaU * deltaU);

            if (pseudoChi2 < 6.f)
                ++nNegMatchedSamplingPoints;
uNegList.push_back(CartesianVector(x, 0., vw2u));
vNegList.push_back(CartesianVector(x, 0., uw2v));
wNegList.push_back(CartesianVector(x, 0., uv2w));
        }
        catch (StatusCodeException &)
        {
        }
    }

//------------------------------------------------------------------------------------------------------------------------------------------
std::cout << " ClusterU: predicted shower shape"  << std::endl;
for (CartesianVectorList::const_iterator iter = uList.begin(); iter != uList.end(); ++iter) PandoraMonitoringApi::AddMarkerToVisualization(&(*iter), "ctrU", CYAN, 1.);
for (CartesianVectorList::const_iterator iter = uPosList.begin(); iter != uPosList.end(); ++iter) PandoraMonitoringApi::AddMarkerToVisualization(&(*iter), "posU", GRAY, 1.);
for (CartesianVectorList::const_iterator iter = uNegList.begin(); iter != uNegList.end(); ++iter) PandoraMonitoringApi::AddMarkerToVisualization(&(*iter), "negU", GRAY, 1.);
PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
ClusterList clusterListU2; clusterListU2.insert(pClusterU);
PandoraMonitoringApi::VisualizeClusters(&clusterListU2, "ClusterListU", RED);
PandoraMonitoringApi::ViewEvent();
std::cout << " ClusterV: predicted shower shape"  << std::endl;
for (CartesianVectorList::const_iterator iter = vList.begin(); iter != vList.end(); ++iter) PandoraMonitoringApi::AddMarkerToVisualization(&(*iter), "ctrV", CYAN, 1.);
for (CartesianVectorList::const_iterator iter = vPosList.begin(); iter != vPosList.end(); ++iter) PandoraMonitoringApi::AddMarkerToVisualization(&(*iter), "posV", GRAY, 1.);
for (CartesianVectorList::const_iterator iter = vNegList.begin(); iter != vNegList.end(); ++iter) PandoraMonitoringApi::AddMarkerToVisualization(&(*iter), "negV", GRAY, 1.);
PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
ClusterList clusterListV2; clusterListV2.insert(pClusterV);
PandoraMonitoringApi::VisualizeClusters(&clusterListV2, "ClusterListV", RED);
PandoraMonitoringApi::ViewEvent();
std::cout << " ClusterW: predicted shower shape"  << std::endl;
for (CartesianVectorList::const_iterator iter = wList.begin(); iter != wList.end(); ++iter) PandoraMonitoringApi::AddMarkerToVisualization(&(*iter), "ctrW", CYAN, 1.);
for (CartesianVectorList::const_iterator iter = wPosList.begin(); iter != wPosList.end(); ++iter) PandoraMonitoringApi::AddMarkerToVisualization(&(*iter), "posW", GRAY, 1.);
for (CartesianVectorList::const_iterator iter = wNegList.begin(); iter != wNegList.end(); ++iter) PandoraMonitoringApi::AddMarkerToVisualization(&(*iter), "negW", GRAY, 1.);
PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
ClusterList clusterListW2; clusterListW2.insert(pClusterW);
PandoraMonitoringApi::VisualizeClusters(&clusterListW2, "ClusterListW", RED);
PandoraMonitoringApi::ViewEvent();
//------------------------------------------------------------------------------------------------------------------------------------------

    const float matchedFraction((nSamplingPoints > 0) ? static_cast<float>(nMatchedSamplingPoints) / static_cast<float>(nSamplingPoints) : 0.f);
    const float posMatchedFraction((nPosSamplingPoints > 0) ? static_cast<float>(nPosMatchedSamplingPoints) / static_cast<float>(nPosSamplingPoints) : 0.f);
    const float negMatchedFraction((nNegSamplingPoints > 0) ? static_cast<float>(nNegMatchedSamplingPoints) / static_cast<float>(nNegSamplingPoints) : 0.f);

std::cout << " All clusters: matchedFraction " << matchedFraction << " posMatchedFraction " << posMatchedFraction << " negMatchedFraction " << negMatchedFraction << std::endl;
PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
ClusterList clusterListU; clusterListU.insert(pClusterU);
PandoraMonitoringApi::VisualizeClusters(&clusterListU, "ClusterListU", RED);
ClusterList clusterListV; clusterListV.insert(pClusterV);
PandoraMonitoringApi::VisualizeClusters(&clusterListV, "ClusterListV", GREEN);
ClusterList clusterListW; clusterListW.insert(pClusterW);
PandoraMonitoringApi::VisualizeClusters(&clusterListW, "ClusterListW", BLUE);
PandoraMonitoringApi::ViewEvent();

    // m_overlapTensor.SetOverlapResult(pClusterU, pClusterV, pClusterW, matchedFraction);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ThreeDShowersAlgorithm::ExamineTensor()
{
    float bestOverlapResult(std::numeric_limits<float>::max()); // TODO Min overlap result for PFO creation
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

                    if (overlapResult < bestOverlapResult)
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
    protoParticle.m_clusterVectorU.push_back(pBestClusterU);
    protoParticle.m_clusterVectorV.push_back(pBestClusterV);
    protoParticle.m_clusterVectorW.push_back(pBestClusterW);
    m_protoParticleVector.push_back(protoParticle);

    // DEBUG
std::cout << " Best particle, overlapResult " << bestOverlapResult << std::endl;
//ClusterList tempListU; tempListU.insert(pBestClusterU);
//ClusterList tempListV; tempListV.insert(pBestClusterV);
//ClusterList tempListW; tempListW.insert(pBestClusterW);
//PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
//PandoraMonitoringApi::VisualizeClusters(&tempListU, "ClusterListU", RED);
//PandoraMonitoringApi::VisualizeClusters(&tempListV, "ClusterListV", GREEN);
//PandoraMonitoringApi::VisualizeClusters(&tempListW, "ClusterListW", BLUE);
//PandoraMonitoringApi::ViewEvent();

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDShowersAlgorithm::CalculateConstantXOverlapResult(const LArClusterHelper::TwoDSlidingFitResult &slidingFitResultU,
    const LArClusterHelper::TwoDSlidingFitResult &slidingFitResultV, const LArClusterHelper::TwoDSlidingFitResult &slidingFitResultW)
{
    std::cout << "TODO - ThreeDShowersAlgorithm::CalculateConstantXOverlapResult " << std::endl;

PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
ClusterList clusterListU; clusterListU.insert(const_cast<Cluster*>(slidingFitResultU.GetCluster()));
PandoraMonitoringApi::VisualizeClusters(&clusterListU, "ClusterListU", RED);
ClusterList clusterListV; clusterListV.insert(const_cast<Cluster*>(slidingFitResultV.GetCluster()));
PandoraMonitoringApi::VisualizeClusters(&clusterListV, "ClusterListV", GREEN);
ClusterList clusterListW; clusterListW.insert(const_cast<Cluster*>(slidingFitResultW.GetCluster()));
PandoraMonitoringApi::VisualizeClusters(&clusterListW, "ClusterListW", BLUE);
PandoraMonitoringApi::ViewEvent();

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDShowersAlgorithm::UpdateTensor()
{
    ClusterList usedClusters;

    for (ProtoParticleVector::const_iterator iter = m_protoParticleVector.begin(), iterEnd = m_protoParticleVector.end(); iter != iterEnd; ++iter)
    {
        usedClusters.insert(iter->m_clusterVectorU.begin(), iter->m_clusterVectorU.end());
        usedClusters.insert(iter->m_clusterVectorV.begin(), iter->m_clusterVectorV.end());
        usedClusters.insert(iter->m_clusterVectorW.begin(), iter->m_clusterVectorW.end());
    }

    for (ClusterList::const_iterator iter = usedClusters.begin(), iterEnd = usedClusters.end(); iter != iterEnd; ++iter)
    {
        m_overlapTensor.RemoveCluster(*iter);
    }

    ThreeDBaseAlgorithm::UpdateTensor();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDShowersAlgorithm::TidyUp()
{
    m_overlapTensor.Clear();
    ThreeDBaseAlgorithm::TidyUp();
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDShowersAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    return ThreeDBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
