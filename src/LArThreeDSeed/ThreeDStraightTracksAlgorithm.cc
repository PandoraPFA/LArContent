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
#include "LArHelpers/LArParticleIdHelper.h"

#include "LArThreeDSeed/ThreeDStraightTracksAlgorithm.h"

using namespace pandora;

namespace lar
{

void ThreeDStraightTracksAlgorithm::SelectInputClusters()
{
    for (ClusterList::const_iterator iter = m_pInputClusterListU->begin(), iterEnd = m_pInputClusterListU->end(); iter != iterEnd; ++iter)
    {
        if ((*iter)->IsAvailable() && (LArParticleIdHelper::LArTrackWidth(*iter) < 0.75f))
            m_clusterVectorU.push_back(*iter);
    }

    for (ClusterList::const_iterator iter = m_pInputClusterListV->begin(), iterEnd = m_pInputClusterListV->end(); iter != iterEnd; ++iter)
    {
        if ((*iter)->IsAvailable() && (LArParticleIdHelper::LArTrackWidth(*iter) < 0.75f))
            m_clusterVectorV.push_back(*iter);
    }

    for (ClusterList::const_iterator iter = m_pInputClusterListW->begin(), iterEnd = m_pInputClusterListW->end(); iter != iterEnd; ++iter)
    {
        if ((*iter)->IsAvailable() && (LArParticleIdHelper::LArTrackWidth(*iter) < 0.75f))
            m_clusterVectorW.push_back(*iter);
    }

    std::sort(m_clusterVectorU.begin(), m_clusterVectorU.end(), LArClusterHelper::SortByNOccupiedLayers);
    std::sort(m_clusterVectorV.begin(), m_clusterVectorV.end(), LArClusterHelper::SortByNOccupiedLayers);
    std::sort(m_clusterVectorW.begin(), m_clusterVectorW.end(), LArClusterHelper::SortByNOccupiedLayers);


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
    ClusterHelper::ClusterFitResult centroidFitResultU;
    ClusterHelper::FitLayerCentroids(pClusterU, pClusterU->GetInnerPseudoLayer(), pClusterU->GetOuterPseudoLayer(), centroidFitResultU);

    const float innerXU(pClusterU->GetCentroid(pClusterU->GetInnerPseudoLayer()).GetX());
    const float outerXU(pClusterU->GetCentroid(pClusterU->GetOuterPseudoLayer()).GetX());
    const float minXU(std::min(innerXU, outerXU));
    const float maxXU(std::max(innerXU, outerXU));

    // V
    ClusterHelper::ClusterFitResult centroidFitResultV;
    ClusterHelper::FitLayerCentroids(pClusterV, pClusterV->GetInnerPseudoLayer(), pClusterV->GetOuterPseudoLayer(), centroidFitResultV);

    const float innerXV(pClusterV->GetCentroid(pClusterV->GetInnerPseudoLayer()).GetX());
    const float outerXV(pClusterV->GetCentroid(pClusterV->GetOuterPseudoLayer()).GetX());
    const float minXV(std::min(innerXV, outerXV));
    const float maxXV(std::max(innerXV, outerXV));

    // W
    ClusterHelper::ClusterFitResult centroidFitResultW;
    ClusterHelper::FitLayerCentroids(pClusterW, pClusterW->GetInnerPseudoLayer(), pClusterW->GetOuterPseudoLayer(), centroidFitResultW);

    const float innerXW(pClusterW->GetCentroid(pClusterW->GetInnerPseudoLayer()).GetX());
    const float outerXW(pClusterW->GetCentroid(pClusterW->GetOuterPseudoLayer()).GetX());
    const float minXW(std::min(innerXW, outerXW));
    const float maxXW(std::max(innerXW, outerXW));

    // Chi2 calculations
    float chi2(std::numeric_limits<float>::max());

    const float minX(std::max(minXU, std::max(minXV, minXW)));
    const float maxX(std::min(maxXU, std::min(maxXV, maxXW)));

    if (minX < maxX)
    {
        chi2 = 0.f;
        static const float xPitch(0.35f);
        unsigned int nSamplingPoints(0);

        for (float x = minX; x < maxX; x += xPitch)
        {
            const float u((centroidFitResultU.GetIntercept() + centroidFitResultU.GetDirection() *
                (x - centroidFitResultU.GetIntercept().GetX()) * (1.f / centroidFitResultU.GetDirection().GetX())).GetZ());
            const float v((centroidFitResultV.GetIntercept() + centroidFitResultV.GetDirection() *
                (x - centroidFitResultV.GetIntercept().GetX()) * (1.f / centroidFitResultV.GetDirection().GetX())).GetZ());
            const float w((centroidFitResultW.GetIntercept() + centroidFitResultW.GetDirection() *
                (x - centroidFitResultW.GetIntercept().GetX()) * (1.f / centroidFitResultW.GetDirection().GetX())).GetZ());

            const float uv2w(LArGeometryHelper::MergeTwoPositions(VIEW_U, VIEW_V, u, v));
            const float uw2v(LArGeometryHelper::MergeTwoPositions(VIEW_U, VIEW_W, u, w));
            const float vw2u(LArGeometryHelper::MergeTwoPositions(VIEW_V, VIEW_W, v, w));

            const float deltaW(uv2w - w), deltaV(uw2v - v), deltaU(vw2u - u);
            chi2 += deltaW * deltaW;
            ++nSamplingPoints;

const CartesianVector expU(x, 0., vw2u); PANDORA_MONITORING_API(AddMarkerToVisualization(&expU, "expU", RED, 1.));
const CartesianVector expV(x, 0., uw2v); PANDORA_MONITORING_API(AddMarkerToVisualization(&expV, "expV", GREEN, 1.));
const CartesianVector expW(x, 0., uv2w); PANDORA_MONITORING_API(AddMarkerToVisualization(&expW, "expW", BLUE, 1.));
        }

        if (nSamplingPoints == 0)
            chi2 = std::numeric_limits<float>::max();

std::cout << " xOverlap " << (maxX - minX) << ", chi2 " << chi2 << std::endl;
PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
ClusterList clusterListU; clusterListU.insert(pClusterU);
PandoraMonitoringApi::VisualizeClusters(&clusterListU, "ClusterListU", RED);
ClusterList clusterListV; clusterListV.insert(pClusterV);
PandoraMonitoringApi::VisualizeClusters(&clusterListV, "ClusterListV", GREEN);
ClusterList clusterListW; clusterListW.insert(pClusterW);
PandoraMonitoringApi::VisualizeClusters(&clusterListW, "ClusterListW", BLUE);
PandoraMonitoringApi::ViewEvent();
    }

    m_overlapTensor.SetOverlapResult(pClusterU, pClusterV, pClusterW, chi2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ThreeDStraightTracksAlgorithm::ExamineTensor()
{
    float bestChi2(std::numeric_limits<float>::max());
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
                const float chi2(m_overlapTensor.GetOverlapResult(*iterU, *iterV, *iterW));

                if (chi2 < bestChi2)
                {
                    bestChi2 = chi2;
                    pBestClusterU = *iterU;
                    pBestClusterV = *iterV;
                    pBestClusterW = *iterW;
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
std::cout << " Best particle, chi2 " << bestChi2 << std::endl;
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
