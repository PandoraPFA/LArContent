/**
 *  @file   LArContent/src/LArClusterSplitting/KinkSplittingAlgorithm.cc
 * 
 *  @brief  Implementation of the kink splitting algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArClusterSplitting/KinkSplittingAlgorithm.h"

#include "LArHelpers/LArClusterHelper.h"

using namespace pandora;

namespace lar
{

bool KinkSplittingAlgorithm::IsPossibleSplit(const Cluster *const pCluster) const
{
    if ((1 + pCluster->GetOuterPseudoLayer() - pCluster->GetInnerPseudoLayer()) < m_minClusterLayers)
        return false;

    ClusterHelper::ClusterFitResult innerLayerFit, outerLayerFit, fullLayerFit;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, ClusterHelper::FitStart(pCluster, m_minClusterLayers, innerLayerFit));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, ClusterHelper::FitEnd(pCluster, m_minClusterLayers, outerLayerFit));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, ClusterHelper::FitFullCluster(pCluster, fullLayerFit));

    if (fullLayerFit.GetRms() < m_minOverallScatteringRms)
        return false; 

    if ((innerLayerFit.GetRms() > m_minVertexScatteringRms) && (outerLayerFit.GetRms() > m_minVertexScatteringRms))
        return false;

// Cluster* tempCluster = (Cluster*)(pCluster);
// ClusterList tempList;
// tempList.insert(tempCluster);
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&tempList, "PossibleSplit", BLUE);
// PandoraMonitoringApi::ViewEvent();
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KinkSplittingAlgorithm::FindBestSplitPosition(const Cluster *const pCluster, CartesianVector &splitPosition) const
{
    LArClusterHelper::TwoDSlidingFitResult twoDSlidingFitResult;
    LArClusterHelper::LArTwoDSlidingXZFit(pCluster, m_slidingFitLayerHalfWindow, twoDSlidingFitResult);

    return twoDSlidingFitResult.FindLargestScatter(splitPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KinkSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_slidingFitLayerHalfWindow = 20;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitLayerHalfWindow", m_slidingFitLayerHalfWindow));

    m_minClusterLayers = 20;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLayers", m_minClusterLayers));

    m_minVertexScatteringRms = 0.25;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinVertexScatteringRms", m_minVertexScatteringRms));

    m_minOverallScatteringRms = 0.75;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinOverallScatteringRms", m_minOverallScatteringRms));

    return ClusterSplittingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
