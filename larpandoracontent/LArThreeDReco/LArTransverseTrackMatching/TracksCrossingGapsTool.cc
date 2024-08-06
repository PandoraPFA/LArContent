/**
 *  @file   larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/TracksCrossingGapsTool.cc
 *
 *  @brief  Implementation of the long tracks tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/LongTracksTool.h"
#include "larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/TracksCrossingGapsTool.h"

using namespace pandora;

namespace lar_content
{

TracksCrossingGapsTool::TracksCrossingGapsTool() :
    m_minMatchedFraction(0.5f),
    m_minMatchedSamplingPoints(10),
    m_minXOverlapFraction(0.9f),
    m_minMatchedSamplingPointRatio(2),
    m_maxGapTolerance(2.f),
    m_sampleStepSize(0.5f),
    m_maxAngleRatio(2)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TracksCrossingGapsTool::Run(ThreeViewTransverseTracksAlgorithm *const pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    if (PandoraContentApi::GetGeometry(*pAlgorithm)->GetDetectorGapList().empty())
        return false;

    ProtoParticleVector protoParticleVector;
    this->FindTracks(pAlgorithm, overlapTensor, protoParticleVector);

    const bool particlesMade(pAlgorithm->CreateThreeDParticles(protoParticleVector));
    return particlesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TracksCrossingGapsTool::FindTracks(
    ThreeViewTransverseTracksAlgorithm *const pAlgorithm, const TensorType &overlapTensor, ProtoParticleVector &protoParticleVector) const
{
    ClusterSet usedClusters;
    ClusterVector sortedKeyClusters;
    overlapTensor.GetSortedKeyClusters(sortedKeyClusters);

    for (const Cluster *const pKeyCluster : sortedKeyClusters)
    {
        if (!pKeyCluster->IsAvailable())
            continue;

        unsigned int nU(0), nV(0), nW(0);
        TensorType::ElementList elementList;
        overlapTensor.GetConnectedElements(pKeyCluster, true, elementList, nU, nV, nW);

        IteratorList iteratorList;
        this->SelectElements(pAlgorithm, elementList, usedClusters, iteratorList);

        // Check that elements are not directly connected and are significantly longer than any other directly connected elements
        for (IteratorList::const_iterator iIter = iteratorList.begin(), iIterEnd = iteratorList.end(); iIter != iIterEnd; ++iIter)
        {
            if (LongTracksTool::HasLongDirectConnections(iIter, iteratorList))
                continue;

            if (!LongTracksTool::IsLongerThanDirectConnections(iIter, elementList, m_minMatchedSamplingPointRatio, usedClusters))
                continue;

            ProtoParticle protoParticle;
            protoParticle.m_clusterList.push_back((*iIter)->GetClusterU());
            protoParticle.m_clusterList.push_back((*iIter)->GetClusterV());
            protoParticle.m_clusterList.push_back((*iIter)->GetClusterW());
            protoParticleVector.push_back(protoParticle);

            usedClusters.insert((*iIter)->GetClusterU());
            usedClusters.insert((*iIter)->GetClusterV());
            usedClusters.insert((*iIter)->GetClusterW());
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TracksCrossingGapsTool::SelectElements(ThreeViewTransverseTracksAlgorithm *const pAlgorithm,
    const TensorType::ElementList &elementList, const pandora::ClusterSet &usedClusters, IteratorList &iteratorList) const
{
    for (TensorType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
    {
        if (usedClusters.count(eIter->GetClusterU()) || usedClusters.count(eIter->GetClusterV()) || usedClusters.count(eIter->GetClusterW()))
            continue;

        if (eIter->GetOverlapResult().GetMatchedFraction() < m_minMatchedFraction)
            continue;

        if (eIter->GetOverlapResult().GetNMatchedSamplingPoints() < m_minMatchedSamplingPoints)
            continue;

        const XOverlap &xOverlap(eIter->GetOverlapResult().GetXOverlap());

        if (xOverlap.GetXOverlapSpan() < std::numeric_limits<float>::epsilon())
            continue;

        // Calculate effective overlap fraction, including the information about gaps
        float xOverlapFractionU(0.f), xOverlapFractionV(0.f), xOverlapFractionW(0.f);
        this->CalculateEffectiveOverlapFractions(pAlgorithm, *eIter, xOverlapFractionU, xOverlapFractionV, xOverlapFractionW);

        if ((xOverlap.GetXSpanU() > std::numeric_limits<float>::epsilon()) && (xOverlapFractionU > m_minXOverlapFraction) &&
            (xOverlap.GetXSpanV() > std::numeric_limits<float>::epsilon()) && (xOverlapFractionV > m_minXOverlapFraction) &&
            (xOverlap.GetXSpanW() > std::numeric_limits<float>::epsilon()) && (xOverlapFractionW > m_minXOverlapFraction))
        {
            iteratorList.push_back(eIter);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TracksCrossingGapsTool::CalculateEffectiveOverlapFractions(ThreeViewTransverseTracksAlgorithm *const pAlgorithm,
    const TensorType::Element &element, float &xOverlapFractionU, float &xOverlapFractionV, float &xOverlapFractionW) const
{
    float xMinEffU(element.GetOverlapResult().GetXOverlap().GetUMinX()), xMaxEffU(element.GetOverlapResult().GetXOverlap().GetUMaxX());
    float xMinEffV(element.GetOverlapResult().GetXOverlap().GetVMinX()), xMaxEffV(element.GetOverlapResult().GetXOverlap().GetVMaxX());
    float xMinEffW(element.GetOverlapResult().GetXOverlap().GetWMinX()), xMaxEffW(element.GetOverlapResult().GetXOverlap().GetWMaxX());
    this->CalculateEffectiveOverlapSpan(pAlgorithm, element, xMinEffU, xMaxEffU, xMinEffV, xMaxEffV, xMinEffW, xMaxEffW);

    const float effectiveXSpanU(xMaxEffU - xMinEffU), effectiveXSpanV(xMaxEffV - xMinEffV), effectiveXSpanW(xMaxEffW - xMinEffW);
    const float minCommonX(std::max(xMinEffU, std::max(xMinEffV, xMinEffW)));
    const float maxCommonX(std::min(xMaxEffU, std::min(xMaxEffV, xMaxEffW)));
    const float effectiveXOverlapSpan(maxCommonX - minCommonX);

    // TODO check that this shouldn't be greater than 1 any more
    xOverlapFractionU = effectiveXSpanU > 0.f ? std::min(1.f, (effectiveXOverlapSpan / effectiveXSpanU)) : 0.f;
    xOverlapFractionV = effectiveXSpanV > 0.f ? std::min(1.f, (effectiveXOverlapSpan / effectiveXSpanV)) : 0.f;
    xOverlapFractionW = effectiveXSpanW > 0.f ? std::min(1.f, (effectiveXOverlapSpan / effectiveXSpanW)) : 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TracksCrossingGapsTool::CalculateEffectiveOverlapSpan(ThreeViewTransverseTracksAlgorithm *const pAlgorithm,
    const TensorType::Element &element, float &xMinEffU, float &xMaxEffU, float &xMinEffV, float &xMaxEffV, float &xMinEffW, float &xMaxEffW) const
{
    const float xMinAll(std::min(xMinEffU, std::min(xMinEffV, xMinEffW)));
    const float xMaxAll(std::max(xMaxEffU, std::max(xMaxEffV, xMaxEffW)));
    const float minCommonX(std::max(xMinEffU, std::max(xMinEffV, xMinEffW)));
    const float maxCommonX(std::min(xMaxEffU, std::min(xMaxEffV, xMaxEffW)));

    float dxUmin(0.f), dxVmin(0.f), dxWmin(0.f);
    float dxUmax(0.f), dxVmax(0.f), dxWmax(0.f);

    // ATTN break out of loops to to avoid finding a non-related gap far from  cluster itself
    const int nSamplingPointsLeft(1 + static_cast<int>((minCommonX - xMinAll) / m_sampleStepSize));
    const int nSamplingPointsRight(1 + static_cast<int>((xMaxAll - maxCommonX) / m_sampleStepSize));

    for (int iSample = 1; iSample <= nSamplingPointsLeft; ++iSample)
    {
        bool gapInU(false), gapInV(false), gapInW(false);
        const float xSample(std::max(xMinAll, minCommonX - static_cast<float>(iSample) * m_sampleStepSize));

        if (!this->PassesGapChecks(pAlgorithm, element, xSample, gapInU, gapInV, gapInW))
            break;

        if (gapInU)
            dxUmin = xMinEffU - xSample;
        if (gapInV)
            dxVmin = xMinEffV - xSample;
        if (gapInW)
            dxWmin = xMinEffW - xSample;
    }

    for (int iSample = 1; iSample <= nSamplingPointsRight; ++iSample)
    {
        bool gapInU(false), gapInV(false), gapInW(false);
        const float xSample(std::min(xMaxAll, maxCommonX + static_cast<float>(iSample) * m_sampleStepSize));

        if (!this->PassesGapChecks(pAlgorithm, element, xSample, gapInU, gapInV, gapInW))
            break;

        if (gapInU)
            dxUmax = xSample - xMaxEffU;
        if (gapInV)
            dxVmax = xSample - xMaxEffV;
        if (gapInW)
            dxWmax = xSample - xMaxEffW;
    }

    xMinEffU -= dxUmin;
    xMaxEffU += dxUmax;
    xMinEffV -= dxVmin;
    xMaxEffV += dxVmax;
    xMinEffW -= dxWmin;
    xMaxEffW += dxWmax;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TracksCrossingGapsTool::PassesGapChecks(ThreeViewTransverseTracksAlgorithm *const pAlgorithm, const TensorType::Element &element,
    const float xSample, bool &gapInU, bool &gapInV, bool &gapInW) const
{
    const TwoDSlidingFitResult &slidingFitResultU(pAlgorithm->GetCachedSlidingFitResult(element.GetClusterU()));
    const TwoDSlidingFitResult &slidingFitResultV(pAlgorithm->GetCachedSlidingFitResult(element.GetClusterV()));
    const TwoDSlidingFitResult &slidingFitResultW(pAlgorithm->GetCachedSlidingFitResult(element.GetClusterW()));

    // If we have access to the global x position in all three clusters, there are no gaps involved (or cluster already spans small gaps)
    CartesianVector fitUPosition(0.f, 0.f, 0.f), fitVPosition(0.f, 0.f, 0.f), fitWPosition(0.f, 0.f, 0.f);
    const StatusCode statusCodeU(slidingFitResultU.GetGlobalFitPositionAtX(xSample, fitUPosition));
    const StatusCode statusCodeV(slidingFitResultV.GetGlobalFitPositionAtX(xSample, fitVPosition));
    const StatusCode statusCodeW(slidingFitResultW.GetGlobalFitPositionAtX(xSample, fitWPosition));

    if ((STATUS_CODE_SUCCESS == statusCodeU) && (STATUS_CODE_SUCCESS == statusCodeV) && (STATUS_CODE_SUCCESS == statusCodeW))
        return false;

    try
    {
        // Note: argument order important - initially assume first view has a gap, but inside CheckXPositionInGap do check other two views
        if ((STATUS_CODE_SUCCESS != statusCodeU) && (!this->IsEndOfCluster(xSample, slidingFitResultU)))
            return this->CheckXPositionInGap(xSample, slidingFitResultU, slidingFitResultV, slidingFitResultW, gapInU, gapInV, gapInW);

        if ((STATUS_CODE_SUCCESS != statusCodeV) && (!this->IsEndOfCluster(xSample, slidingFitResultV)))
            return this->CheckXPositionInGap(xSample, slidingFitResultV, slidingFitResultU, slidingFitResultW, gapInV, gapInU, gapInW);

        if ((STATUS_CODE_SUCCESS != statusCodeW) && (!this->IsEndOfCluster(xSample, slidingFitResultW)))
            return this->CheckXPositionInGap(xSample, slidingFitResultW, slidingFitResultU, slidingFitResultV, gapInW, gapInU, gapInV);
    }
    catch (const StatusCodeException &statusCodeException)
    {
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TracksCrossingGapsTool::CheckXPositionInGap(const float xSample, const TwoDSlidingFitResult &slidingFitResult1,
    const TwoDSlidingFitResult &slidingFitResult2, const TwoDSlidingFitResult &slidingFitResult3, bool &gapIn1, bool &gapIn2, bool &gapIn3) const
{
    CartesianVector fitPosition2(0.f, 0.f, 0.f), fitPosition3(0.f, 0.f, 0.f);

    // If we have the global position at X from the two other clusters, calculate projection in the first view and check for gaps
    if ((STATUS_CODE_SUCCESS == slidingFitResult2.GetGlobalFitPositionAtX(xSample, fitPosition2)) &&
        (STATUS_CODE_SUCCESS == slidingFitResult3.GetGlobalFitPositionAtX(xSample, fitPosition3)))
    {
        const HitType hitType1(LArClusterHelper::GetClusterHitType(slidingFitResult1.GetCluster()));
        const HitType hitType2(LArClusterHelper::GetClusterHitType(slidingFitResult2.GetCluster()));
        const HitType hitType3(LArClusterHelper::GetClusterHitType(slidingFitResult3.GetCluster()));

        const float zSample(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), hitType2, hitType3, fitPosition2.GetZ(), fitPosition3.GetZ()));
        const CartesianVector samplingPoint(xSample, 0.f, zSample);
        return LArGeometryHelper::IsInGap(this->GetPandora(), CartesianVector(xSample, 0.f, zSample), hitType1, m_maxGapTolerance);
    }

    // ATTN Only safe to return here (for efficiency) because gapIn2 and gapIn3 values aren't used by calling function if we return false
    gapIn1 = LArGeometryHelper::IsXSamplingPointInGap(this->GetPandora(), xSample, slidingFitResult1, m_sampleStepSize);

    if (!gapIn1)
        return false;

    // If we dont have a projection at x in the other two clusters, check if they are in gaps or at the end of the cluster
    if ((STATUS_CODE_SUCCESS != slidingFitResult2.GetGlobalFitPositionAtX(xSample, fitPosition2)) &&
        (STATUS_CODE_SUCCESS != slidingFitResult3.GetGlobalFitPositionAtX(xSample, fitPosition3)))
    {
        const bool endIn2(this->IsEndOfCluster(xSample, slidingFitResult2));
        const bool endIn3(this->IsEndOfCluster(xSample, slidingFitResult3));

        if (!endIn2)
            gapIn2 = LArGeometryHelper::IsXSamplingPointInGap(this->GetPandora(), xSample, slidingFitResult2, m_sampleStepSize);

        if (!endIn3)
            gapIn3 = LArGeometryHelper::IsXSamplingPointInGap(this->GetPandora(), xSample, slidingFitResult3, m_sampleStepSize);

        return ((gapIn2 && endIn3) || (gapIn3 && endIn2) || (endIn2 && endIn3));
    }

    // Finally, check whether there is a second gap involved
    if (STATUS_CODE_SUCCESS != slidingFitResult2.GetGlobalFitPositionAtX(xSample, fitPosition2))
    {
        gapIn2 = LArGeometryHelper::IsXSamplingPointInGap(this->GetPandora(), xSample, slidingFitResult2, m_sampleStepSize);
        return (gapIn2 || this->IsEndOfCluster(xSample, slidingFitResult2));
    }
    else
    {
        gapIn3 = LArGeometryHelper::IsXSamplingPointInGap(this->GetPandora(), xSample, slidingFitResult3, m_sampleStepSize);
        return (gapIn3 || this->IsEndOfCluster(xSample, slidingFitResult3));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TracksCrossingGapsTool::IsEndOfCluster(const float xSample, const TwoDSlidingFitResult &slidingFitResult) const
{
    return ((std::fabs(slidingFitResult.GetGlobalMinLayerPosition().GetX() - xSample) < slidingFitResult.GetLayerPitch()) ||
        (std::fabs(slidingFitResult.GetGlobalMaxLayerPosition().GetX() - xSample) < slidingFitResult.GetLayerPitch()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TracksCrossingGapsTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinMatchedFraction", m_minMatchedFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinMatchedSamplingPoints", m_minMatchedSamplingPoints));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinXOverlapFraction", m_minXOverlapFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinMatchedSamplingPointRatio", m_minMatchedSamplingPointRatio));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxGapTolerance", m_maxGapTolerance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SampleStepSize", m_sampleStepSize));

    if (m_sampleStepSize < std::numeric_limits<float>::epsilon())
    {
        std::cout << "TracksCrossingGapsTool: Invalid value for SampleStepSize " << m_sampleStepSize << std::endl;
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxAngleRatio", m_maxAngleRatio));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
