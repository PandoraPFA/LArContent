/**
 *  @file   LArContent/src/LArThreeDReco/LArTransverseTrackMatching/MissingTrackSegmentTool.cc
 * 
 *  @brief  Implementation of the missing track segment tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"
#include "LArHelpers/LArPointingClusterHelper.h"

#include "LArObjects/LArPointingCluster.h"

#include "LArThreeDReco/LArTransverseTrackMatching/LongTracksTool.h"
#include "LArThreeDReco/LArTransverseTrackMatching/MissingTrackSegmentTool.h"

using namespace pandora;

namespace lar
{

bool MissingTrackSegmentTool::Run(ThreeDTransverseTracksAlgorithm *pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraSettings::ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << m_algorithmToolType << std::endl;

    ProtoParticleVector protoParticleVector; ClusterMergeMap clusterMergeMap;
    this->FindTracks(pAlgorithm, overlapTensor, protoParticleVector, clusterMergeMap);

    const bool particlesMade(pAlgorithm->CreateThreeDParticles(protoParticleVector));
    const bool mergesMade(pAlgorithm->MakeClusterMerges(clusterMergeMap));

    return (particlesMade || mergesMade);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MissingTrackSegmentTool::FindTracks(ThreeDTransverseTracksAlgorithm *pAlgorithm, const TensorType &overlapTensor,
    ProtoParticleVector &protoParticleVector, ClusterMergeMap &clusterMergeMap) const
{
    ClusterList usedClusters;

    for (TensorType::const_iterator iterU = overlapTensor.begin(), iterUEnd = overlapTensor.end(); iterU != iterUEnd; ++iterU)
    {
        if (!iterU->first->IsAvailable())
            continue;

        unsigned int nU(0), nV(0), nW(0);
        TensorType::ElementList elementList;
        overlapTensor.GetConnectedElements(iterU->first, true, elementList, nU, nV, nW);

        IteratorList iteratorList;
        this->SelectElements(elementList, usedClusters, iteratorList);

        for (IteratorList::const_iterator iIter = iteratorList.begin(), iIterEnd = iteratorList.end(); iIter != iIterEnd; ++iIter)
        {
            if (LongTracksTool::HasLongDirectConnections(iIter, iteratorList))
                continue;

            if (!LongTracksTool::IsLongerThanDirectConnections(iIter, elementList, m_minMatchedSamplingPointRatio, usedClusters))
                continue;

            if (!this->PassesParticleChecks(pAlgorithm, *(*iIter), usedClusters, clusterMergeMap))
                continue;

            ProtoParticle protoParticle;
            protoParticle.m_clusterListU.insert((*iIter)->GetClusterU());
            protoParticle.m_clusterListV.insert((*iIter)->GetClusterV());
            protoParticle.m_clusterListW.insert((*iIter)->GetClusterW());
            protoParticleVector.push_back(protoParticle);

            usedClusters.insert((*iIter)->GetClusterU());
            usedClusters.insert((*iIter)->GetClusterV());
            usedClusters.insert((*iIter)->GetClusterW());
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MissingTrackSegmentTool::SelectElements(const TensorType::ElementList &elementList, const pandora::ClusterList &usedClusters, IteratorList &iteratorList) const
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
        const float shortSpan(std::min(xOverlap.GetXSpanU(), std::min(xOverlap.GetXSpanV(), xOverlap.GetXSpanW())));
        const float longSpan1(std::max(xOverlap.GetXSpanU(), std::max(xOverlap.GetXSpanV(), xOverlap.GetXSpanW())));
        const float longSpan2(((xOverlap.GetXSpanU() > shortSpan) && (xOverlap.GetXSpanU() < longSpan1)) ? xOverlap.GetXSpanU() :
            ((xOverlap.GetXSpanV() > shortSpan) && (xOverlap.GetXSpanV() < longSpan1)) ? xOverlap.GetXSpanV() : xOverlap.GetXSpanW());

        if ((xOverlap.GetXOverlapSpan() < std::numeric_limits<float>::epsilon()) || (longSpan1 < std::numeric_limits<float>::epsilon()))
            continue;

        if (((shortSpan / xOverlap.GetXOverlapSpan()) < m_minInitialXOverlapFraction) || ((longSpan2 / longSpan1) < m_minInitialXOverlapFraction))
            continue;

        iteratorList.push_back(eIter);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MissingTrackSegmentTool::PassesParticleChecks(ThreeDTransverseTracksAlgorithm *pAlgorithm, const TensorType::Element &element,
    ClusterList &usedClusters, ClusterMergeMap &clusterMergeMap) const
{
    const Particle particle(element);

    ClusterList candidateClusters;
    this->GetCandidateClusters(pAlgorithm, particle, candidateClusters);

    if (candidateClusters.empty())
        return false;

    SlidingFitResultMap slidingFitResultMap;
    this->GetSlidingFitResultMap(pAlgorithm, candidateClusters, slidingFitResultMap);

    if (slidingFitResultMap.empty())
        return false;

    SegmentOverlapMap segmentOverlapMap;
    this->GetSegmentOverlapMap(pAlgorithm, particle, slidingFitResultMap, segmentOverlapMap);

    if (segmentOverlapMap.empty())
        return false;

    return this->MakeDecisions(particle, slidingFitResultMap, segmentOverlapMap, usedClusters, clusterMergeMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MissingTrackSegmentTool::GetCandidateClusters(ThreeDTransverseTracksAlgorithm *pAlgorithm, const Particle &particle,
    ClusterList &candidateClusters) const
{
    const ClusterList &clusterList((TPC_VIEW_U == particle.m_shortHitType) ? pAlgorithm->GetInputClusterListU() :
        (TPC_VIEW_V == particle.m_shortHitType) ? pAlgorithm->GetInputClusterListV() : pAlgorithm->GetInputClusterListW());

    for (ClusterList::const_iterator iter = clusterList.begin(), iterEnd = clusterList.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster(*iter);

        if (pCluster == particle.m_pShortCluster)
            continue;

        if (pCluster->GetNCaloHits() < m_minCaloHitsInCandidateCluster)
            continue;

        candidateClusters.insert(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MissingTrackSegmentTool::GetSlidingFitResultMap(ThreeDTransverseTracksAlgorithm *pAlgorithm, const ClusterList &candidateClusterList,
    SlidingFitResultMap &slidingFitResultMap) const
{
    for (ClusterList::const_iterator iter = candidateClusterList.begin(), iterEnd = candidateClusterList.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster(*iter);

        try
        {
            const TwoDSlidingFitResult &slidingFitResult(pAlgorithm->GetCachedSlidingFitResult(pCluster));
            slidingFitResultMap[pCluster] = slidingFitResult;
            continue;
        }
        catch (StatusCodeException &)
        {
        }

        try
        {
            TwoDSlidingFitResult slidingFitResult;
            LArClusterHelper::LArTwoDSlidingFit(pCluster, pAlgorithm->GetSlidingFitWindow(), slidingFitResult);
            slidingFitResultMap[pCluster] = slidingFitResult;
            continue;
        }
        catch (StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MissingTrackSegmentTool::GetSegmentOverlapMap(ThreeDTransverseTracksAlgorithm *pAlgorithm, const Particle &particle,
    const SlidingFitResultMap &slidingFitResultMap, SegmentOverlapMap &segmentOverlapMap) const
{
    const TwoDSlidingFitResult &fitResult1(pAlgorithm->GetCachedSlidingFitResult(particle.m_pCluster1));
    const TwoDSlidingFitResult &fitResult2(pAlgorithm->GetCachedSlidingFitResult(particle.m_pCluster2));

    const float nPoints1(std::fabs(static_cast<float>(fitResult1.GetMaxLayer() - fitResult1.GetMinLayer())));
    const float nPoints2(std::fabs(static_cast<float>(fitResult2.GetMaxLayer() - fitResult2.GetMinLayer())));

    const unsigned int nPoints(static_cast<unsigned int>(1.f + (nPoints1 + nPoints2) / 2.f));

    for (unsigned n = 0; n <= nPoints; ++n)
    {
        const float x(particle.m_longMinX + (particle.m_longMaxX - particle.m_longMinX) * static_cast<float>(n) / static_cast<float>(nPoints));

        if ((x > particle.m_shortMinX) && (x < particle.m_shortMaxX))
            continue;

        try
        {
            CartesianVector fitVector1(0.f, 0.f, 0.f), fitVector2(0.f, 0.f, 0.f);
            fitResult1.GetGlobalFitPositionAtX(x, fitVector1);
            fitResult2.GetGlobalFitPositionAtX(x, fitVector2);
            const float prediction(LArGeometryHelper::MergeTwoPositions(particle.m_hitType1, particle.m_hitType2, fitVector1.GetZ(), fitVector2.GetZ()));

            for (SlidingFitResultMap::const_iterator iter = slidingFitResultMap.begin(), iterEnd = slidingFitResultMap.end(); iter != iterEnd; ++iter)
            {
                try
                {
                    CartesianVector fitVector(0.f, 0.f, 0.f), fitDirection(0.f, 0.f, 0.f);
                    iter->second.GetGlobalFitPositionAtX(x, fitVector);
                    iter->second.GetGlobalFitDirectionAtX(x, fitDirection);

                    const float delta((prediction - fitVector.GetZ()) * fitDirection.GetX());
                    const float pseudoChi2(delta * delta);

                    SegmentOverlap &segmentOverlap(segmentOverlapMap[iter->first]);
                    ++segmentOverlap.m_nSamplingPoints;
                    segmentOverlap.m_pseudoChi2Sum += pseudoChi2;

                    if (pseudoChi2 < m_pseudoChi2Cut)
                    {
                        ++segmentOverlap.m_nMatchedSamplingPoints;
                        segmentOverlap.m_matchedSamplingMinX = std::min(x, segmentOverlap.m_matchedSamplingMinX);
                        segmentOverlap.m_matchedSamplingMaxX = std::max(x, segmentOverlap.m_matchedSamplingMaxX);
                    }
                }
                catch (StatusCodeException &)
                {
                }
            }
        }
        catch (StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MissingTrackSegmentTool::MakeDecisions(const Particle &particle, const SlidingFitResultMap &slidingFitResultMap,
    const SegmentOverlapMap &segmentOverlapMap, ClusterList &usedClusters, ClusterMergeMap &clusterMergeMap) const
{
    ClusterList possibleMerges;
    float shortMinX(particle.m_shortMinX), shortMaxX(particle.m_shortMaxX);
    bool matchesACluster(false);

    for (SegmentOverlapMap::const_iterator iter = segmentOverlapMap.begin(), iterEnd = segmentOverlapMap.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster(iter->first);
        const SegmentOverlap &segmentOverlap(iter->second);

        if (!this->PassesSamplingCuts(segmentOverlap))
            continue;

        shortMinX = std::min(segmentOverlap.m_matchedSamplingMinX, shortMinX);
        shortMaxX = std::max(segmentOverlap.m_matchedSamplingMaxX, shortMaxX);
        matchesACluster = true;

        // Allow pfo construction if find hits in an unavailable cluster, but can't merge unavailable cluster into this pfo
        if (!usedClusters.insert(pCluster).second || !pCluster->IsAvailable())
            continue;

        if (!this->IsPossibleMerge(pCluster, particle, segmentOverlap, slidingFitResultMap))
            continue;

        possibleMerges.insert(pCluster);
    }

    if (std::fabs(particle.m_longMaxX - particle.m_longMinX) < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    if (!matchesACluster || possibleMerges.empty())
        return false;

    if (((shortMaxX - shortMinX) / (particle.m_longMaxX - particle.m_longMinX)) < m_minFinalXOverlapFraction)
        return false;

    clusterMergeMap[particle.m_pShortCluster].insert(possibleMerges.begin(), possibleMerges.end());
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MissingTrackSegmentTool::PassesSamplingCuts(const SegmentOverlap &segmentOverlap) const
{
    if (0 == segmentOverlap.m_nSamplingPoints)
        return false;

    if ((segmentOverlap.m_nSamplingPoints < m_makePfoMinSamplingPoints) || (segmentOverlap.m_nMatchedSamplingPoints < m_makePfoMinMatchedSamplingPoints))
        return false;

    if ((static_cast<float>(segmentOverlap.m_nMatchedSamplingPoints) / static_cast<float>(segmentOverlap.m_nSamplingPoints)) < m_makePfoMinMatchedFraction)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MissingTrackSegmentTool::IsPossibleMerge(Cluster *const pCluster, const Particle &particle, const SegmentOverlap &segmentOverlap,
    const SlidingFitResultMap &slidingFitResultMap) const
{
    if ((segmentOverlap.m_pseudoChi2Sum / static_cast<float>(segmentOverlap.m_nSamplingPoints)) > m_mergeMaxChi2PerSamplingPoint)
        return false;

    SlidingFitResultMap::const_iterator fitIter = slidingFitResultMap.find(pCluster);

    if (slidingFitResultMap.end() == fitIter)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    float mergeMinX(std::numeric_limits<float>::max()), mergeMaxX(-std::numeric_limits<float>::max());
    fitIter->second.GetMinAndMaxX(mergeMinX, mergeMaxX);

    // cluster should not be wider than the longest span
    if ((mergeMinX < particle.m_longMinX - m_mergeXContainmentTolerance) || (mergeMaxX > particle.m_longMaxX + m_mergeXContainmentTolerance))
        return false;

    // cluster should not overlap with the shortest span
    if ((mergeMinX < particle.m_shortMaxX - m_mergeXContainmentTolerance) && (mergeMaxX > particle.m_shortMinX + m_mergeXContainmentTolerance))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

MissingTrackSegmentTool::Particle::Particle(const TensorType::Element &element)
{
    const XOverlap &xOverlap(element.GetOverlapResult().GetXOverlap());

    m_shortHitType = ((xOverlap.GetXSpanU() < xOverlap.GetXSpanV()) && (xOverlap.GetXSpanU() < xOverlap.GetXSpanW())) ? TPC_VIEW_U :
        ((xOverlap.GetXSpanV() < xOverlap.GetXSpanU()) && (xOverlap.GetXSpanV() < xOverlap.GetXSpanW())) ? TPC_VIEW_V :
        ((xOverlap.GetXSpanW() < xOverlap.GetXSpanU()) && (xOverlap.GetXSpanW() < xOverlap.GetXSpanV())) ? TPC_VIEW_W : CUSTOM;

    if (CUSTOM == m_shortHitType)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    m_pShortCluster = (TPC_VIEW_U == m_shortHitType) ? element.GetClusterU() : (TPC_VIEW_V == m_shortHitType) ? element.GetClusterV() : element.GetClusterW();
    m_pCluster1 = (TPC_VIEW_U == m_shortHitType) ? element.GetClusterV() : element.GetClusterU();
    m_pCluster2 = (TPC_VIEW_W == m_shortHitType) ? element.GetClusterV() : element.GetClusterW();
    m_shortMinX = (TPC_VIEW_U == m_shortHitType) ? xOverlap.GetUMinX() : (TPC_VIEW_V == m_shortHitType) ? xOverlap.GetVMinX() : xOverlap.GetWMinX();
    m_shortMaxX = (TPC_VIEW_U == m_shortHitType) ? xOverlap.GetUMaxX() : (TPC_VIEW_V == m_shortHitType) ? xOverlap.GetVMaxX() : xOverlap.GetWMaxX();
    m_longMinX = (TPC_VIEW_U == m_shortHitType) ? std::min(xOverlap.GetVMinX(), xOverlap.GetWMinX()) : (TPC_VIEW_V == m_shortHitType) ? std::min(xOverlap.GetUMinX(), xOverlap.GetWMinX()) : std::min(xOverlap.GetUMinX(), xOverlap.GetVMinX());
    m_longMaxX = (TPC_VIEW_U == m_shortHitType) ? std::max(xOverlap.GetVMaxX(), xOverlap.GetWMaxX()) : (TPC_VIEW_V == m_shortHitType) ? std::max(xOverlap.GetUMaxX(), xOverlap.GetWMaxX()) : std::max(xOverlap.GetUMaxX(), xOverlap.GetVMaxX());

    m_hitType1 = LArClusterHelper::GetClusterHitType(m_pCluster1);
    m_hitType2 = LArClusterHelper::GetClusterHitType(m_pCluster2);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MissingTrackSegmentTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_minMatchedFraction = 0.9f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedFraction", m_minMatchedFraction));

    m_minMatchedSamplingPoints = 10;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedSamplingPoints", m_minMatchedSamplingPoints));

    m_minMatchedSamplingPointRatio = 2;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedSamplingPointRatio", m_minMatchedSamplingPointRatio));

    m_minInitialXOverlapFraction = 0.75f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinInitialXOverlapFraction", m_minInitialXOverlapFraction));

    m_minFinalXOverlapFraction = 0.75f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinFinalXOverlapFraction", m_minFinalXOverlapFraction));

    m_minCaloHitsInCandidateCluster = 5;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCaloHitsInCandidateCluster", m_minCaloHitsInCandidateCluster));

    m_pseudoChi2Cut = 1.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PseudoChi2Cut", m_pseudoChi2Cut));

    m_makePfoMinSamplingPoints = 5;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MakePfoMinSamplingPoints", m_makePfoMinSamplingPoints));

    m_makePfoMinMatchedSamplingPoints = 5;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MakePfoMinMatchedSamplingPoints", m_makePfoMinMatchedSamplingPoints));

    m_makePfoMinMatchedFraction = 0.8f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MakePfoMinMatchedFraction", m_makePfoMinMatchedFraction));

    m_makePfoMaxImpactParameter = 3.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MakePfoMaxImpactParameter", m_makePfoMaxImpactParameter));

    m_mergeMaxChi2PerSamplingPoint = 0.25f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MergeMaxChi2PerSamplingPoint", m_mergeMaxChi2PerSamplingPoint));

    m_mergeXContainmentTolerance = 1.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MergeXContainmentTolerance", m_mergeXContainmentTolerance));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
