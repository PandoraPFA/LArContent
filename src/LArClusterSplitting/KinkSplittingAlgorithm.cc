/**
 *  @file   LArContent/src/LArClusterSplitting/KinkSplittingAlgorithm.cc
 * 
 *  @brief  Implementation of the kink splitting algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArClusterSplitting/KinkSplittingAlgorithm.h"

#include "LArHelpers/LArParticleIdHelper.h"

using namespace pandora;

namespace lar
{

StatusCode KinkSplittingAlgorithm::Run()
{
    return ClusterSplittingAlgorithm::Run();
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool KinkSplittingAlgorithm::IsPossibleSplit(const Cluster *const pCluster) const
{
    if ((1 + pCluster->GetOuterPseudoLayer() - pCluster->GetInnerPseudoLayer()) < m_minClusterLayers)
        return false;

    ClusterHelper::ClusterFitResult innerLayerFit, outerLayerFit;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, ClusterHelper::FitStart(pCluster, m_minClusterLayers, innerLayerFit));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, ClusterHelper::FitEnd(pCluster, m_minClusterLayers, outerLayerFit));

    if (innerLayerFit.GetDirection().GetCosOpeningAngle(outerLayerFit.GetDirection()) > m_maxCosScatteringAngle)
        return false;

    if ((innerLayerFit.GetRms() > m_minScatteringRms) && (outerLayerFit.GetRms() > m_minScatteringRms))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KinkSplittingAlgorithm::FindBestSplitLayer(const Cluster* const pCluster, unsigned int& splitLayer )
{ 
    LArParticleIdHelper::TwoDSlidingXZFitResult twoDSlidingXZFitResult;
    LArParticleIdHelper::LArTwoDSlidingXZFit(pCluster, twoDSlidingXZFitResult);

    return twoDSlidingXZFitResult.FindLargestScatter(splitLayer);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KinkSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_minClusterLayers = 20;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLayers", m_minClusterLayers));

    m_minScatteringRms = 0.15;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinScatteringRms", m_minScatteringRms));

    m_maxCosScatteringAngle = std::cos(M_PI * 15.f / 180.f);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxCosScatteringAngle", m_maxCosScatteringAngle));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
