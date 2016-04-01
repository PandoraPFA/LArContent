/**
 *  @file   LArContent/src/LArTwoDReco/LArClusterAssociation/CrossGapsAssociationAlgorithm.cc
 *
 *  @brief  Implementation of the cross gaps association algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"

#include "LArTwoDReco/LArClusterAssociation/CrossGapsAssociationAlgorithm.h"

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
    m_minMatchedSamplingFraction(0.5f)
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
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));

    for (const Cluster *const pCluster : clusterVector)
    {
        try {(void) slidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pCluster, TwoDSlidingFitResult(pCluster, m_slidingFitWindow, slidingFitPitch)));}
        catch (StatusCodeException &) {}
    }
//ClusterList temp;
//for (const auto &map : slidingFitResultMap) temp.insert(map.first);
//PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f);
//PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &temp, "Clusters", BLUE);
//PandoraMonitoringApi::ViewEvent(this->GetPandora());

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

bool CrossGapsAssociationAlgorithm::IsExtremalCluster(const bool isForward, const Cluster *const pCurrentCluster,  const Cluster *const pTestCluster) const
{
    const unsigned int currentLayer(isForward ? pCurrentCluster->GetOuterPseudoLayer() : pCurrentCluster->GetInnerPseudoLayer());
    const unsigned int testLayer(isForward ? pTestCluster->GetOuterPseudoLayer() : pTestCluster->GetInnerPseudoLayer());

    const float currentEnergy(pCurrentCluster->GetHadronicEnergy());
    const float testEnergy(pTestCluster->GetHadronicEnergy());

    if (isForward && ((testLayer > currentLayer) || ((testLayer == currentLayer) && (testEnergy > currentEnergy))))
        return true;

    if (!isForward && ((testLayer < currentLayer) || ((testLayer == currentLayer) && (testEnergy > currentEnergy))))
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

//ClusterList tempList1, tempList2;
//tempList1.insert(innerFitResult.GetCluster());
//tempList2.insert(outerFitResult.GetCluster());
//std::cout << " innerHits " << innerFitResult.GetCluster()->GetNCaloHits() << ", outerHits " << outerFitResult.GetCluster()->GetNCaloHits() << std::endl;
//PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f);
//PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &tempList1, "InnerCluster", BLUE);
//PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &tempList2, "OuterCluster", GREEN);
//PandoraMonitoringApi::ViewEvent(this->GetPandora());
    return (this->IsAssociated(innerFitResult.GetGlobalMaxLayerPosition(), innerFitResult.GetGlobalMaxLayerDirection(), outerFitResult) &&
        this->IsAssociated(outerFitResult.GetGlobalMinLayerPosition(), outerFitResult.GetGlobalMinLayerDirection() * -1.f, innerFitResult));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CrossGapsAssociationAlgorithm::IsAssociated(const CartesianVector &startPosition, const CartesianVector &startDirection,
    const TwoDSlidingFitResult &targetFitResult) const
{
    const HitType hitType(LArClusterHelper::GetClusterHitType(targetFitResult.GetCluster()));
    unsigned int nSamplingPoints(0), nGapSamplingPoints(0), nMatchedSamplingPoints(0), nUnmatchedSampleRun(0);

    for (unsigned int iSample = 0; iSample < m_maxSamplingPoints; ++iSample)
    {
        ++nSamplingPoints;
        const CartesianVector samplingPoint(startPosition + startDirection * static_cast<float>(iSample) * m_sampleStepSize);

        if (this->IsInGap(samplingPoint, hitType))
        {
//PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &samplingPoint, "samplingPoint_" + TypeToString(iSample), GRAY, 1);
            ++nGapSamplingPoints;
            nUnmatchedSampleRun = 0; // ATTN Choose to also reset run when entering gap region 
            continue;
        }

        if (this->IsNearCluster(samplingPoint, targetFitResult))
        {
//PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &samplingPoint, "samplingPoint_" + TypeToString(iSample), RED, 1);
            ++nMatchedSamplingPoints;
            nUnmatchedSampleRun = 0;
        }
        else
        {
//PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &samplingPoint, "samplingPoint_" + TypeToString(iSample), ORANGE, 1);
            if (++nUnmatchedSampleRun > m_maxUnmatchedSampleRun)
                break;
        }
    }

    const float expectation((targetFitResult.GetGlobalMaxLayerPosition() - targetFitResult.GetGlobalMinLayerPosition()).GetMagnitude() / m_sampleStepSize);
    const float matchedSamplingFraction(expectation > 0.f ? static_cast<float>(nMatchedSamplingPoints) / expectation : 0.f);
//std::cout << "nSamplingPoints " << nSamplingPoints << " nGapSamplingPoints " << nGapSamplingPoints << " nMatchedSamplingPoints " << nMatchedSamplingPoints << std::endl;
//std::cout << "non-gap matched fraction " << static_cast<float>(nMatchedSamplingPoints) / static_cast<float>(nSamplingPoints - nGapSamplingPoints) << std::endl;
//std::cout << "expectation " << expectation << " nMatchedSamplingPoints " << nMatchedSamplingPoints << std::endl;
//PandoraMonitoringApi::ViewEvent(this->GetPandora());
    if ((nMatchedSamplingPoints > m_minMatchedSamplingPoints) || (matchedSamplingFraction > m_minMatchedSamplingFraction))
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CrossGapsAssociationAlgorithm::IsInGap(const CartesianVector &samplingPoint, const HitType hitType) const
{
    for (const DetectorGap *const pDetectorGap : PandoraContentApi::GetGeometry(*this)->GetDetectorGapList())
    {
        if (pDetectorGap->IsInGap(samplingPoint, hitType))
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CrossGapsAssociationAlgorithm::IsNearCluster(const CartesianVector &samplingPoint, const TwoDSlidingFitResult &targetFitResult) const
{
    float rL(std::numeric_limits<float>::max()), rT(std::numeric_limits<float>::max());
    targetFitResult.GetLocalPosition(samplingPoint, rL, rT);

    CartesianVector fitPosition(0.f, 0.f, 0.f);

    if (STATUS_CODE_SUCCESS == targetFitResult.GetGlobalFitPosition(rL, fitPosition))
    {
        if ((fitPosition - samplingPoint).GetMagnitudeSquared() < m_maxOnClusterDistance * m_maxOnClusterDistance)
            return true;
    }

    CartesianVector fitPositionAtX(0.f, 0.f, 0.f);

    if (STATUS_CODE_SUCCESS == targetFitResult.GetGlobalFitPositionAtX(samplingPoint.GetX(), fitPositionAtX))
    {
        if ((fitPositionAtX - samplingPoint).GetMagnitudeSquared() < m_maxOnClusterDistance * m_maxOnClusterDistance)
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CrossGapsAssociationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterHits", m_minClusterHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLayers", m_minClusterLayers));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxSamplingPoints", m_maxSamplingPoints));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SampleStepSize", m_sampleStepSize));

    if (m_sampleStepSize < std::numeric_limits<float>::epsilon())
    {
        std::cout << "CrossGapsAssociationAlgorithm: Invalid value for SampleStepSize " << m_sampleStepSize << std::endl;
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxUnmatchedSampleRun", m_maxUnmatchedSampleRun));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxOnClusterDistance", m_maxOnClusterDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedSamplingPoints", m_minMatchedSamplingPoints));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedSamplingFraction", m_minMatchedSamplingFraction));

    return ClusterAssociationAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
