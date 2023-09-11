/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/CrossGapsAssociationAlgorithm.cc
 *
 *  @brief  Implementation of the cross gaps association algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/CrossGapsAssociationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CrossGapsAssociationAlgorithm::CrossGapsAssociationAlgorithm() :
    m_minClusterHits(10),
    m_minClusterLayers(6),
    m_slidingFitWindow(20),
    m_maxSamplingPoints(1000),
    m_sampleStepSize(0.5f),
    m_maxUnmatchedSampleRun(8),
    m_maxOnClusterDistance(1.5f),
    m_minMatchedSamplingPoints(10),
    m_minMatchedSamplingFraction(0.5f),
    m_gapTolerance(0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CrossGapsAssociationAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    // ATTN May want to opt-out completely if no gap information available
    // if (PandoraContentApi::GetGeometry(*this)->GetDetectorGapList().empty())
    //     return;

    for (const Cluster *const pCluster : *pClusterList)
    {
        if (pCluster->GetNCaloHits() < m_minClusterHits)
            continue;

        if (1 + pCluster->GetOuterPseudoLayer() - pCluster->GetInnerPseudoLayer() < m_minClusterLayers)
            continue;

        clusterVector.push_back(pCluster);
    }

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByInnerLayer);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CrossGapsAssociationAlgorithm::PopulateClusterAssociationMap(const ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const
{
    TwoDSlidingFitResultMap slidingFitResultMap;

    for (const Cluster *const pCluster : clusterVector)
    {
        try
        {
            const float slidingFitPitch(LArGeometryHelper::GetWirePitch(this->GetPandora(), LArClusterHelper::GetClusterHitType(pCluster)));
            slidingFitResultMap.insert(
                TwoDSlidingFitResultMap::value_type(pCluster, TwoDSlidingFitResult(pCluster, m_slidingFitWindow, slidingFitPitch)));
        }
        catch (StatusCodeException &)
        {
        }
    }

    // ATTN This method assumes that clusters have been sorted by layer
    for (ClusterVector::const_iterator iterI = clusterVector.begin(), iterIEnd = clusterVector.end(); iterI != iterIEnd; ++iterI)
    {
        const Cluster *const pInnerCluster = *iterI;
        TwoDSlidingFitResultMap::const_iterator fitIterI = slidingFitResultMap.find(pInnerCluster);

        if (slidingFitResultMap.end() == fitIterI)
            continue;

        for (ClusterVector::const_iterator iterJ = iterI, iterJEnd = clusterVector.end(); iterJ != iterJEnd; ++iterJ)
        {
            const Cluster *const pOuterCluster = *iterJ;

            if (pInnerCluster == pOuterCluster)
                continue;

            TwoDSlidingFitResultMap::const_iterator fitIterJ = slidingFitResultMap.find(pOuterCluster);

            if (slidingFitResultMap.end() == fitIterJ)
                continue;

            if (!this->AreClustersAssociated(fitIterI->second, fitIterJ->second))
                continue;

            clusterAssociationMap[pInnerCluster].m_forwardAssociations.insert(pOuterCluster);
            clusterAssociationMap[pOuterCluster].m_backwardAssociations.insert(pInnerCluster);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CrossGapsAssociationAlgorithm::IsExtremalCluster(const bool isForward, const Cluster *const pCurrentCluster, const Cluster *const pTestCluster) const
{
    const unsigned int currentLayer(isForward ? pCurrentCluster->GetOuterPseudoLayer() : pCurrentCluster->GetInnerPseudoLayer());
    const unsigned int testLayer(isForward ? pTestCluster->GetOuterPseudoLayer() : pTestCluster->GetInnerPseudoLayer());

    if (isForward && ((testLayer > currentLayer) || ((testLayer == currentLayer) && LArClusterHelper::SortByNHits(pTestCluster, pCurrentCluster))))
        return true;

    if (!isForward && ((testLayer < currentLayer) || ((testLayer == currentLayer) && LArClusterHelper::SortByNHits(pTestCluster, pCurrentCluster))))
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CrossGapsAssociationAlgorithm::AreClustersAssociated(const TwoDSlidingFitResult &innerFitResult, const TwoDSlidingFitResult &outerFitResult) const
{
    if (outerFitResult.GetCluster()->GetInnerPseudoLayer() < innerFitResult.GetCluster()->GetInnerPseudoLayer())
        throw pandora::StatusCodeException(STATUS_CODE_NOT_ALLOWED);

    if (outerFitResult.GetCluster()->GetInnerPseudoLayer() < innerFitResult.GetCluster()->GetOuterPseudoLayer())
        return false;

    return (this->IsAssociated(innerFitResult.GetGlobalMaxLayerPosition(), innerFitResult.GetGlobalMaxLayerDirection(), outerFitResult) &&
            this->IsAssociated(outerFitResult.GetGlobalMinLayerPosition(), outerFitResult.GetGlobalMinLayerDirection() * -1.f, innerFitResult));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CrossGapsAssociationAlgorithm::IsAssociated(
    const CartesianVector &startPosition, const CartesianVector &startDirection, const TwoDSlidingFitResult &targetFitResult) const
{
    const HitType hitType(LArClusterHelper::GetClusterHitType(targetFitResult.GetCluster()));
    const float ratio{LArGeometryHelper::GetWirePitchRatio(this->GetPandora(), hitType)};
    const float sampleStepSizeAdjusted{ratio * m_sampleStepSize};
    unsigned int nSamplingPoints(0), nGapSamplingPoints(0), nMatchedSamplingPoints(0), nUnmatchedSampleRun(0);

    for (unsigned int iSample = 0; iSample < m_maxSamplingPoints; ++iSample)
    {
        ++nSamplingPoints;
        const CartesianVector samplingPoint(startPosition + startDirection * static_cast<float>(iSample) * sampleStepSizeAdjusted);

        if (LArGeometryHelper::IsInGap(this->GetPandora(), samplingPoint, hitType, m_gapTolerance))
        {
            ++nGapSamplingPoints;
            nUnmatchedSampleRun = 0; // ATTN Choose to also reset run when entering gap region
            continue;
        }

        if (this->IsNearCluster(samplingPoint, targetFitResult))
        {
            ++nMatchedSamplingPoints;
            nUnmatchedSampleRun = 0;
        }
        else if (++nUnmatchedSampleRun > m_maxUnmatchedSampleRun)
        {
            break;
        }
    }

    const float expectation((targetFitResult.GetGlobalMaxLayerPosition() - targetFitResult.GetGlobalMinLayerPosition()).GetMagnitude() / sampleStepSizeAdjusted);
    const float matchedSamplingFraction(expectation > 0.f ? static_cast<float>(nMatchedSamplingPoints) / expectation : 0.f);

    if ((nMatchedSamplingPoints > m_minMatchedSamplingPoints) || (matchedSamplingFraction > m_minMatchedSamplingFraction))
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CrossGapsAssociationAlgorithm::IsNearCluster(const CartesianVector &samplingPoint, const TwoDSlidingFitResult &targetFitResult) const
{
    const HitType hitType(LArClusterHelper::GetClusterHitType(targetFitResult.GetCluster()));
    const float ratio{LArGeometryHelper::GetWirePitchRatio(this->GetPandora(), hitType)};
    const float maxOnClusterDistanceAdjusted{ratio * m_maxOnClusterDistance};

    float rL(std::numeric_limits<float>::max()), rT(std::numeric_limits<float>::max());
    targetFitResult.GetLocalPosition(samplingPoint, rL, rT);

    CartesianVector fitPosition(0.f, 0.f, 0.f);

    if (STATUS_CODE_SUCCESS == targetFitResult.GetGlobalFitPosition(rL, fitPosition))
    {
        if ((fitPosition - samplingPoint).GetMagnitudeSquared() < maxOnClusterDistanceAdjusted * maxOnClusterDistanceAdjusted)
            return true;
    }

    CartesianVector fitPositionAtX(0.f, 0.f, 0.f);

    if (STATUS_CODE_SUCCESS == targetFitResult.GetGlobalFitPositionAtX(samplingPoint.GetX(), fitPositionAtX))
    {
        if ((fitPositionAtX - samplingPoint).GetMagnitudeSquared() < maxOnClusterDistanceAdjusted * maxOnClusterDistanceAdjusted)
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CrossGapsAssociationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterHits", m_minClusterHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterLayers", m_minClusterLayers));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxSamplingPoints", m_maxSamplingPoints));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SampleStepSize", m_sampleStepSize));

    if (m_sampleStepSize < std::numeric_limits<float>::epsilon())
    {
        std::cout << "CrossGapsAssociationAlgorithm: Invalid value for SampleStepSize " << m_sampleStepSize << std::endl;
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxUnmatchedSampleRun", m_maxUnmatchedSampleRun));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxOnClusterDistance", m_maxOnClusterDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinMatchedSamplingPoints", m_minMatchedSamplingPoints));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinMatchedSamplingFraction", m_minMatchedSamplingFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "GapTolerance", m_gapTolerance));

    return ClusterAssociationAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
