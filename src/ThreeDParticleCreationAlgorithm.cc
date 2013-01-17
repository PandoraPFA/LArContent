/**
 *  @file   ThreeDParticleCreationAlgorithm.cc
 * 
 *  @brief  Implementation of the three dimensional particle creation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArClusterHelper.h"
#include "LArGeometryHelper.h"

#include "ThreeDParticleCreationAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode ThreeDParticleCreationAlgorithm::Run()
{
    // Obtain sorted vectors of seed clusters in U, V and W views
    const ClusterList *pClusterListU(NULL), *pClusterListV(NULL), *pClusterListW(NULL);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetClusterList(*this, m_inputClusterListNameU, pClusterListU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetClusterList(*this, m_inputClusterListNameV, pClusterListV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetClusterList(*this, m_inputClusterListNameW, pClusterListW));

    ClusterVector clusterVectorU(pClusterListU->begin(), pClusterListU->end());
    ClusterVector clusterVectorV(pClusterListV->begin(), pClusterListV->end());
    ClusterVector clusterVectorW(pClusterListW->begin(), pClusterListW->end());

    std::sort(clusterVectorU.begin(), clusterVectorU.end(), LArClusterHelper::SortByNOccupiedLayers);
    std::sort(clusterVectorV.begin(), clusterVectorV.end(), LArClusterHelper::SortByNOccupiedLayers);
    std::sort(clusterVectorW.begin(), clusterVectorW.end(), LArClusterHelper::SortByNOccupiedLayers);

    // Investigate combinations of seed clusters in U, V and W views
    ClusterOverlapTensor clusterOverlapTensor;

    for (ClusterVector::const_iterator iterU = clusterVectorU.begin(), iterUEnd = clusterVectorU.end(); iterU != iterUEnd; ++iterU)
    {
        for (ClusterVector::const_iterator iterV = clusterVectorV.begin(), iterVEnd = clusterVectorV.end(); iterV != iterVEnd; ++iterV)
        {
            for (ClusterVector::const_iterator iterW = clusterVectorW.begin(), iterWEnd = clusterVectorW.end(); iterW != iterWEnd; ++iterW)
            {
                clusterOverlapTensor.SetClusterOverlapInfo(*iterU, *iterV, *iterW);
            }
        }
    }

    // Examine output
//    for (ClusterVector::const_iterator iterU = clusterVectorU.begin(), iterUEnd = clusterVectorU.end(); iterU != iterUEnd; ++iterU)
//    {
//        for (ClusterVector::const_iterator iterV = clusterVectorV.begin(), iterVEnd = clusterVectorV.end(); iterV != iterVEnd; ++iterV)
//        {
//            for (ClusterVector::const_iterator iterW = clusterVectorW.begin(), iterWEnd = clusterVectorW.end(); iterW != iterWEnd; ++iterW)
//            {
//                const ClusterOverlapInfo &clusterOverlapInfo(clusterOverlapTensor.GetClusterOverlapInfo(*iterU, *iterV, *iterW));
//                std::cout << " Clusters " <<  *iterU << ", " << *iterV << ", " << *iterW
//                          << ", chi2UVW " << clusterOverlapInfo.GetChi2UVW()
//                          << ", chi2UWV " << clusterOverlapInfo.GetChi2UWV()
//                          << ", chi2VWU " << clusterOverlapInfo.GetChi2VWU() << std::endl;
//                ClusterList clusterListU; clusterListU.insert(*iterU);
//                ClusterList clusterListV; clusterListV.insert(*iterV);
//                ClusterList clusterListW; clusterListW.insert(*iterW);
//                PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
//                PandoraMonitoringApi::VisualizeClusters(&clusterListU, "ClusterListU", RED);
//                PandoraMonitoringApi::VisualizeClusters(&clusterListV, "ClusterListV", BLUE);
//                PandoraMonitoringApi::VisualizeClusters(&clusterListW, "ClusterListW", GREEN);
//                PandoraMonitoringApi::ViewEvent();
//            }
//        }
//    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDParticleCreationAlgorithm::ClusterOverlapInfo::ClusterOverlapInfo(Cluster *pClusterU, Cluster *pClusterV, Cluster *pClusterW) :
    m_chi2UVW(-1.f),
    m_chi2UWV(-1.f),
    m_chi2VWU(-1.f),
    m_xOverlap(0.f)
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
    const float minX(std::max(minXU, std::max(minXV, minXW)));
    const float maxX(std::min(maxXU, std::min(maxXV, maxXW)));

    if (minX < maxX)
    {
        static const float xPitch(0.35f);
        unsigned int nSamplingPoints(0);
        m_xOverlap = maxX - minX;

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
            m_chi2UVW += deltaW * deltaW;
            m_chi2UWV += deltaV * deltaV;
            m_chi2VWU += deltaU * deltaU;
            ++nSamplingPoints;

            // Debug display
            const CartesianVector fitU(x, 0., vw2u); PANDORA_MONITORING_API(AddMarkerToVisualization(&fitU, "fitU", GREEN, 1.));
            const CartesianVector fitV(x, 0., uw2v); PANDORA_MONITORING_API(AddMarkerToVisualization(&fitV, "fitV", RED, 1.));
            const CartesianVector fitW(x, 0., uv2w); PANDORA_MONITORING_API(AddMarkerToVisualization(&fitW, "fitW", BLUE, 1.));
        }

        if (nSamplingPoints > 0)
        {
            m_chi2UVW /= static_cast<float>(nSamplingPoints);
            m_chi2UWV /= static_cast<float>(nSamplingPoints);
            m_chi2VWU /= static_cast<float>(nSamplingPoints);
        }

        std::cout << " chi2UVW " << m_chi2UVW << " chi2UWV " << m_chi2UWV << " chi2VWU " << m_chi2VWU << " sum " << m_chi2UVW + m_chi2UWV + m_chi2VWU << std::endl;

        PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
        ClusterList clusterListU; clusterListU.insert(pClusterU);
        PandoraMonitoringApi::VisualizeClusters(&clusterListU, "ClusterListU", GREEN);
        ClusterList clusterListV; clusterListV.insert(pClusterV);
        PandoraMonitoringApi::VisualizeClusters(&clusterListV, "ClusterListV", RED);
        ClusterList clusterListW; clusterListW.insert(pClusterW);
        PandoraMonitoringApi::VisualizeClusters(&clusterListW, "ClusterListW", BLUE);
        PandoraMonitoringApi::ViewEvent();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

const ThreeDParticleCreationAlgorithm::ClusterOverlapMatrix &ThreeDParticleCreationAlgorithm::ClusterOverlapTensor::GetClusterOverlapMatrix(
    Cluster *pClusterU) const
{
    OverlapTensor::const_iterator iter = m_overlapTensor.find(pClusterU);

    if (m_overlapTensor.end() == iter)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const ThreeDParticleCreationAlgorithm::ClusterOverlapList &ThreeDParticleCreationAlgorithm::ClusterOverlapTensor::GetClusterOverlapList(
    Cluster *pClusterU, Cluster *pClusterV) const
{
    const ClusterOverlapMatrix &overlapMatrix(this->GetClusterOverlapMatrix(pClusterU));
    ClusterOverlapMatrix::const_iterator iter = overlapMatrix.find(pClusterV);

    if (overlapMatrix.end() == iter)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const ThreeDParticleCreationAlgorithm::ClusterOverlapInfo &ThreeDParticleCreationAlgorithm::ClusterOverlapTensor::GetClusterOverlapInfo(
    Cluster *pClusterU, Cluster *pClusterV, Cluster *pClusterW) const
{
    const ClusterOverlapList &overlapList(this->GetClusterOverlapList(pClusterU, pClusterV));
    ClusterOverlapList::const_iterator iter = overlapList.find(pClusterW);

    if (overlapList.end() == iter)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDParticleCreationAlgorithm::ClusterOverlapTensor::SetClusterOverlapInfo(Cluster *pClusterU, Cluster *pClusterV, Cluster *pClusterW)
{
    ClusterOverlapList &clusterOverlapList = m_overlapTensor[pClusterU][pClusterV];
    ClusterOverlapList::const_iterator iter = clusterOverlapList.find(pClusterW);

    if (clusterOverlapList.end() != iter)
        throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);

    if (!clusterOverlapList.insert(ClusterOverlapList::value_type(pClusterW, ClusterOverlapInfo(pClusterU, pClusterV, pClusterW))).second)
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDParticleCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameU", m_inputClusterListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameV", m_inputClusterListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameW", m_inputClusterListNameW));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputPfoListName", m_outputPfoListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
