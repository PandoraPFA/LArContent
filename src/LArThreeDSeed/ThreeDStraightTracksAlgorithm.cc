/**
 *  @file   LArContent/src/LArThreeDSeed/ThreeDStraightTracksAlgorithm.cc
 * 
 *  @brief  Implementation of the three dimension straight tracks algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"

#include "LArThreeDSeed/ThreeDStraightTracksAlgorithm.h"

using namespace pandora;

namespace lar
{

void ThreeDStraightTracksAlgorithm::SelectInputClusters()
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

void ThreeDStraightTracksAlgorithm::SelectInputClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (!pCluster->IsAvailable())
            continue;

        if (LArClusterHelper::GetLayerOccupancy(pCluster) < 0.75f)
            continue;

        LArClusterHelper::TwoDSlidingFitResult twoDSlidingFitResult;
        LArClusterHelper::LArTwoDSlidingFit(pCluster, 20, twoDSlidingFitResult);

        FloatVector residuals;
        unsigned int nHitsOnTrack(0);

        const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
        const LArClusterHelper::TwoDSlidingFitResult::LayerFitResultMap &layerFitResultMap(twoDSlidingFitResult.GetLayerFitResultMap());

        for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(); iter != orderedCaloHitList.end(); ++iter)
        {
            for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
            {
                float rL(0.f), rT(0.f);
                twoDSlidingFitResult.GetLocalCoordinates((*hitIter)->GetPositionVector(), rL, rT);
                const int layer(twoDSlidingFitResult.GetLayer(rL));

                LArClusterHelper::TwoDSlidingFitResult::LayerFitResultMap::const_iterator fitResultIter = layerFitResultMap.find(layer);

                if (layerFitResultMap.end() == fitResultIter)
                    continue;

                const double fitT(fitResultIter->second.GetFitT());
                const double gradient(fitResultIter->second.GetGradient());
                const double residualSquared((fitT - rT) * (fitT - rT) / (1. + gradient * gradient)); // angular correction (note: this is cheating!)
                residuals.push_back(residualSquared);

                if (residualSquared < 1.f)
                    ++nHitsOnTrack;
            }
        }

        if (residuals.empty())
           continue;

        std::sort(residuals.begin(), residuals.end());
        static const float m_trackResidualQuantile(0.8f);
        const float theQuantile(residuals[m_trackResidualQuantile * residuals.size()]);
        const float slidingFitWidth(std::sqrt(theQuantile));
std::cout << " SELECT? slidingFitWidth " << slidingFitWidth << " hitsOnTrackFraction " << (static_cast<float>(nHitsOnTrack) / (static_cast<float>(pCluster->GetNCaloHits()))) << std::endl;

        if (slidingFitWidth > 1.5f) // TODO
            continue;

        const float hitsOnTrackFraction(static_cast<float>(nHitsOnTrack) / (static_cast<float>(pCluster->GetNCaloHits())));

        if (hitsOnTrackFraction < 0.9f) // TODO
            continue;

        clusterVector.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDStraightTracksAlgorithm::ModifyInputClusters()
{
    // TODO, see ShowerMipSeparationAlgorithm for some basic algorithm mechanics
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDStraightTracksAlgorithm::InitializeTensor()
{
    m_overlapTensor.Clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDStraightTracksAlgorithm::CalculateOverlapResult(Cluster *pClusterU, Cluster *pClusterV, Cluster *pClusterW)
{
    // U
    LArClusterHelper::TwoDSlidingFitResult slidingFitResultU;
    LArClusterHelper::LArTwoDSlidingFit(pClusterU, 20, slidingFitResultU);

    const float innerXU(slidingFitResultU.GetGlobalMinLayerPosition().GetX());
    const float outerXU(slidingFitResultU.GetGlobalMaxLayerPosition().GetX());
    const float minXU(std::min(innerXU, outerXU));
    const float maxXU(std::max(innerXU, outerXU));
    const float xSpanU(maxXU - minXU);

    // V
    LArClusterHelper::TwoDSlidingFitResult slidingFitResultV;
    LArClusterHelper::LArTwoDSlidingFit(pClusterV, 20, slidingFitResultV);

    const float innerXV(slidingFitResultV.GetGlobalMinLayerPosition().GetX());
    const float outerXV(slidingFitResultV.GetGlobalMaxLayerPosition().GetX());
    const float minXV(std::min(innerXV, outerXV));
    const float maxXV(std::max(innerXV, outerXV));
    const float xSpanV(maxXV - minXV);

    // W
    LArClusterHelper::TwoDSlidingFitResult slidingFitResultW;
    LArClusterHelper::LArTwoDSlidingFit(pClusterW, 20, slidingFitResultW);

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

    if (slidingFitResultU.IsMultivaluedInX() || slidingFitResultV.IsMultivaluedInX() || slidingFitResultW.IsMultivaluedInX())
        return this->CalculateConstantXOverlapResult(slidingFitResultU, slidingFitResultV, slidingFitResultW);

    // Sampling in x
    const float nPointsU((xOverlap / xSpanU) * pClusterU->GetNCaloHits());
    const float nPointsV((xOverlap / xSpanV) * pClusterV->GetNCaloHits());
    const float nPointsW((xOverlap / xSpanW) * pClusterW->GetNCaloHits());
    const float xPitch(3.f * xOverlap / (nPointsU + nPointsV + nPointsW));

    // Chi2 calculations
    unsigned int nSamplingPoints(0), nMatchedSamplingPoints(0);

    for (float x = minX; x < maxX; x += xPitch)
    {
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

            if (pseudoChi2 < 3.f) // TODO
                ++nMatchedSamplingPoints;

//PANDORA_MONITORING_API(AddMarkerToVisualization(&fitUVector, "fitUVector", RED, 1.));
//PANDORA_MONITORING_API(AddMarkerToVisualization(&fitVVector, "fitVVector", GREEN, 1.));
//PANDORA_MONITORING_API(AddMarkerToVisualization(&fitWVector, "fitWVector", BLUE, 1.));
const CartesianVector expU(x, 0., vw2u); PANDORA_MONITORING_API(AddMarkerToVisualization(&expU, "expU", RED, 1.));
const CartesianVector expV(x, 0., uw2v); PANDORA_MONITORING_API(AddMarkerToVisualization(&expV, "expV", GREEN, 1.));
const CartesianVector expW(x, 0., uv2w); PANDORA_MONITORING_API(AddMarkerToVisualization(&expW, "expW", BLUE, 1.));
        }
        catch (StatusCodeException &)
        {
        }
    }

    if (0 == nSamplingPoints)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const float matchedSamplingFraction(static_cast<float>(nMatchedSamplingPoints) / static_cast<float>(nSamplingPoints));

PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
ClusterList clusterListU; clusterListU.insert(pClusterU);
PandoraMonitoringApi::VisualizeClusters(&clusterListU, "ClusterListU", RED);
ClusterList clusterListV; clusterListV.insert(pClusterV);
PandoraMonitoringApi::VisualizeClusters(&clusterListV, "ClusterListV", GREEN);
ClusterList clusterListW; clusterListW.insert(pClusterW);
PandoraMonitoringApi::VisualizeClusters(&clusterListW, "ClusterListW", BLUE);

    if (matchedSamplingFraction < 0.9f) // TODO
{
std::cout << " VETO: matchedSamplingFraction " << matchedSamplingFraction << std::endl;
PandoraMonitoringApi::ViewEvent();
        return;
}

std::cout << " POPULATE TENSOR: xOverlap " << xOverlap << ", xOverlapU " << (xOverlap / xSpanU) << ", xOverlapV " << (xOverlap / xSpanV) << ", xOverlapW " << (xOverlap / xSpanW) << ", nMatchedSamplingPoints " << nMatchedSamplingPoints << ", nSamplingPoints " << nSamplingPoints << ", matchedSamplingFraction " << matchedSamplingFraction << std::endl;
PandoraMonitoringApi::ViewEvent();
    m_overlapTensor.SetOverlapResult(pClusterU, pClusterV, pClusterW, matchedSamplingFraction);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ThreeDStraightTracksAlgorithm::ExamineTensor()
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
ClusterList tempListU; tempListU.insert(pBestClusterU);
ClusterList tempListV; tempListV.insert(pBestClusterV);
ClusterList tempListW; tempListW.insert(pBestClusterW);
PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
PandoraMonitoringApi::VisualizeClusters(&tempListU, "ClusterListU", RED);
PandoraMonitoringApi::VisualizeClusters(&tempListV, "ClusterListV", GREEN);
PandoraMonitoringApi::VisualizeClusters(&tempListW, "ClusterListW", BLUE);
PandoraMonitoringApi::ViewEvent();

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDStraightTracksAlgorithm::CalculateConstantXOverlapResult(const LArClusterHelper::TwoDSlidingFitResult &slidingFitResultU,
    const LArClusterHelper::TwoDSlidingFitResult &slidingFitResultV, const LArClusterHelper::TwoDSlidingFitResult &slidingFitResultW)
{
    std::cout << "TODO - ThreeDStraightTracksAlgorithm::CalculateConstantXOverlapResult " << std::endl;

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

void ThreeDStraightTracksAlgorithm::UpdateTensor()
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

void ThreeDStraightTracksAlgorithm::TidyUp()
{
    m_overlapTensor.Clear();
    ThreeDBaseAlgorithm::TidyUp();
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDStraightTracksAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    return ThreeDBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
