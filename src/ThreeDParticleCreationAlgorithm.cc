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

//    // Examine output
//    for (ClusterOverlapTensor::const_iterator iterU = clusterOverlapTensor.begin(), iterUEnd = clusterOverlapTensor.end(); iterU != iterUEnd; ++iterU)
//    {
//        for (ClusterOverlapMatrix::const_iterator iterV = iterU->second.begin(), iterVEnd = iterU->second.end(); iterV != iterVEnd; ++iterV)
//        {
//            for (ClusterOverlapList::const_iterator iterW = iterV->second.begin(), iterWEnd = iterV->second.end(); iterW != iterWEnd; ++iterW)
//            {
//                const ClusterOverlapInfo &clusterOverlapInfo(iterW->second);
//                std::cout << " Clusters " << iterU->first << ", " << iterV->first << ", " << iterW->first
//                          << ", xOverlap "<< clusterOverlapInfo.GetXOverlap()
//                          << ", chi2UVW " << clusterOverlapInfo.GetChi2UVW()
//                          << ", chi2UWV " << clusterOverlapInfo.GetChi2UWV()
//                          << ", chi2VWU " << clusterOverlapInfo.GetChi2VWU()
//                          << ", chi2Sum " << clusterOverlapInfo.GetChi2Sum() << std::endl;
//                ClusterList clusterListU; clusterListU.insert(iterU->first);
//                ClusterList clusterListV; clusterListV.insert(iterV->first);
//                ClusterList clusterListW; clusterListW.insert(iterW->first);
//                PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
//                PandoraMonitoringApi::VisualizeClusters(&clusterListU, "ClusterListU", RED);
//                PandoraMonitoringApi::VisualizeClusters(&clusterListV, "ClusterListV", BLUE);
//                PandoraMonitoringApi::VisualizeClusters(&clusterListW, "ClusterListW", GREEN);
//                PandoraMonitoringApi::ViewEvent();
//            }
//        }
//    }

    // Form 3D particles
    const PfoList *pPfoList = NULL; std::string pfoListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryPfoListAndSetCurrent(*this, pPfoList, pfoListName));

    ClusterList usedClusters;

    while(true)
    {
        float bestChi2Sum(std::numeric_limits<float>::max());
        Cluster *pBestClusterU(NULL), *pBestClusterV(NULL), *pBestClusterW(NULL);

        for (ClusterOverlapTensor::const_iterator iterU = clusterOverlapTensor.begin(), iterUEnd = clusterOverlapTensor.end(); iterU != iterUEnd; ++iterU)
        {
            for (ClusterOverlapMatrix::const_iterator iterV = iterU->second.begin(), iterVEnd = iterU->second.end(); iterV != iterVEnd; ++iterV)
            {
                for (ClusterOverlapList::const_iterator iterW = iterV->second.begin(), iterWEnd = iterV->second.end(); iterW != iterWEnd; ++iterW)
                {
                    const ClusterOverlapInfo &clusterOverlapInfo(iterW->second);

                    if ((clusterOverlapInfo.GetChi2Sum() < bestChi2Sum) &&
                        (0 == usedClusters.count(iterU->first)) && (0 == usedClusters.count(iterV->first)) && (0 == usedClusters.count(iterW->first)))
                    {
                        bestChi2Sum = clusterOverlapInfo.GetChi2Sum();
                        pBestClusterU = iterU->first;
                        pBestClusterV = iterV->first;
                        pBestClusterW = iterW->first;
                    }
                }
            }
        }

        if (!pBestClusterU || !pBestClusterV || !pBestClusterW)
            break;

        // TODO - correct parameters
        PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
        pfoParameters.m_particleId = 22;
        pfoParameters.m_charge = 0;
        pfoParameters.m_mass = 0.f;
        pfoParameters.m_energy = 0.f;
        pfoParameters.m_momentum = CartesianVector(0., 0., 0.);
        pfoParameters.m_clusterList.insert(pBestClusterU);
        pfoParameters.m_clusterList.insert(pBestClusterV);
        pfoParameters.m_clusterList.insert(pBestClusterW);
        usedClusters.insert(pfoParameters.m_clusterList.begin(), pfoParameters.m_clusterList.end());
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters));

        // DEBUG
std::cout << " Make particle, best Chi2Sum " << bestChi2Sum << std::endl;
ClusterList clusterListU; clusterListU.insert(pBestClusterU);
ClusterList clusterListV; clusterListV.insert(pBestClusterV);
ClusterList clusterListW; clusterListW.insert(pBestClusterW);
PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
PandoraMonitoringApi::VisualizeClusters(&clusterListU, "ClusterListU", RED);
PandoraMonitoringApi::VisualizeClusters(&clusterListV, "ClusterListV", GREEN);
PandoraMonitoringApi::VisualizeClusters(&clusterListW, "ClusterListW", BLUE);
PandoraMonitoringApi::ViewEvent();
    }

    if (!pPfoList->empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SavePfoList(*this, m_outputPfoListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentPfoList(*this, m_outputPfoListName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDParticleCreationAlgorithm::ClusterOverlapInfo::ClusterOverlapInfo(Cluster *pClusterU, Cluster *pClusterV, Cluster *pClusterW) :
    m_isInitialized(false),
    m_xOverlap(0.f),
    m_chi2UVW(-1.f),
    m_chi2UWV(-1.f),
    m_chi2VWU(-1.f)
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
const CartesianVector expU(x, 0., vw2u); PANDORA_MONITORING_API(AddMarkerToVisualization(&expU, "expU", RED, 1.));
const CartesianVector expV(x, 0., uw2v); PANDORA_MONITORING_API(AddMarkerToVisualization(&expV, "expV", GREEN, 1.));
const CartesianVector expW(x, 0., uv2w); PANDORA_MONITORING_API(AddMarkerToVisualization(&expW, "expW", BLUE, 1.));
        }

        if (nSamplingPoints > 0)
        {
            m_isInitialized = true;
            m_xOverlap = maxX - minX;
            m_chi2UVW /= static_cast<float>(nSamplingPoints);
            m_chi2UWV /= static_cast<float>(nSamplingPoints);
            m_chi2VWU /= static_cast<float>(nSamplingPoints);
        }

std::cout << " xOverlap " << m_xOverlap << ", chi2UVW " << m_chi2UVW << ", chi2UWV " << m_chi2UWV << ", chi2VWU " << m_chi2VWU
          << ", sum " << (m_chi2UVW + m_chi2UWV + m_chi2VWU) << std::endl;
//for (OrderedCaloHitList::const_iterator iter = pClusterU->GetOrderedCaloHitList().begin(); iter != pClusterU->GetOrderedCaloHitList().end(); ++iter)
//{
//    const CartesianVector centroid(pClusterU->GetCentroid(iter->first));
//    PANDORA_MONITORING_API(AddMarkerToVisualization(&centroid, "centroidU", DARKRED, 1.));
//}
//for (float x = minXU; x < maxXU; x += xPitch)
//{
//    const CartesianVector fit(centroidFitResultU.GetIntercept() + centroidFitResultU.GetDirection() *
//        (x - centroidFitResultU.GetIntercept().GetX()) * (1.f / centroidFitResultU.GetDirection().GetX()));
//    PANDORA_MONITORING_API(AddMarkerToVisualization(&fit, "fitU", RED, 1.));
//}
//for (OrderedCaloHitList::const_iterator iter = pClusterV->GetOrderedCaloHitList().begin(); iter != pClusterV->GetOrderedCaloHitList().end(); ++iter)
//{
//    const CartesianVector centroid(pClusterV->GetCentroid(iter->first));
//    PANDORA_MONITORING_API(AddMarkerToVisualization(&centroid, "centroidU", DARKGREEN, 1.));
//}
//for (float x = minXV; x < maxXV; x += xPitch)
//{
//    const CartesianVector fit(centroidFitResultV.GetIntercept() + centroidFitResultV.GetDirection() *
//        (x - centroidFitResultV.GetIntercept().GetX()) * (1.f / centroidFitResultV.GetDirection().GetX()));
//    PANDORA_MONITORING_API(AddMarkerToVisualization(&fit, "fitV", GREEN, 1.));
//}
//for (OrderedCaloHitList::const_iterator iter = pClusterW->GetOrderedCaloHitList().begin(); iter != pClusterW->GetOrderedCaloHitList().end(); ++iter)
//{
//    const CartesianVector centroid(pClusterW->GetCentroid(iter->first));
//    PANDORA_MONITORING_API(AddMarkerToVisualization(&centroid, "centroidU", DARKBLUE, 1.));
//}
//for (float x = minXW; x < maxXW; x += xPitch)
//{
//    const CartesianVector fit(centroidFitResultW.GetIntercept() + centroidFitResultW.GetDirection() *
//        (x - centroidFitResultW.GetIntercept().GetX()) * (1.f / centroidFitResultW.GetDirection().GetX()));
//    PANDORA_MONITORING_API(AddMarkerToVisualization(&fit, "fitU", BLUE, 1.));
//}
PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
ClusterList clusterListU; clusterListU.insert(pClusterU);
PandoraMonitoringApi::VisualizeClusters(&clusterListU, "ClusterListU", RED);
ClusterList clusterListV; clusterListV.insert(pClusterV);
PandoraMonitoringApi::VisualizeClusters(&clusterListV, "ClusterListV", GREEN);
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
    ClusterOverlapTensor::const_iterator iter = this->find(pClusterU);

    if (this->end() == iter)
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
    ClusterOverlapList &clusterOverlapList = (*this)[pClusterU][pClusterV];
    ClusterOverlapList::const_iterator iter = clusterOverlapList.find(pClusterW);

    if (clusterOverlapList.end() != iter)
        throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);

    const ClusterOverlapInfo clusterOverlapInfo(pClusterU, pClusterV, pClusterW);

    if (!clusterOverlapInfo.IsInitialized())
        return;

    if (!clusterOverlapList.insert(ClusterOverlapList::value_type(pClusterW, clusterOverlapInfo)).second)
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
