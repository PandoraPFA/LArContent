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

    if (m_constantXTreatment && (slidingFitResultU.IsMultivaluedInX() || slidingFitResultV.IsMultivaluedInX() || slidingFitResultW.IsMultivaluedInX()))
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

void ThreeDTracksAlgorithm::CalculateConstantXOverlapResult(const LArClusterHelper::TwoDSlidingFitResult &slidingFitResultU,
    const LArClusterHelper::TwoDSlidingFitResult &slidingFitResultV, const LArClusterHelper::TwoDSlidingFitResult &slidingFitResultW)
{
    std::cout << "TODO - ThreeDTracksAlgorithm::CalculateConstantXOverlapResult " << std::endl;
    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDTracksAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_constantXTreatment = false;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ConstantXTreatment", m_constantXTreatment));

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
