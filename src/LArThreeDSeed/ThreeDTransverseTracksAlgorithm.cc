/**
 *  @file   LArContent/src/LArThreeDSeed/ThreeDTransverseTracksAlgorithm.cc
 * 
 *  @brief  Implementation of the three dimensional transverse tracks algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"

#include "LArThreeDSeed/ThreeDTransverseTracksAlgorithm.h"

using namespace pandora;

namespace lar
{

void ThreeDTransverseTracksAlgorithm::PreparationStep()
{
    ClusterList allClustersList;
    allClustersList.insert(m_clusterVectorU.begin(), m_clusterVectorU.end());
    allClustersList.insert(m_clusterVectorV.begin(), m_clusterVectorV.end());
    allClustersList.insert(m_clusterVectorW.begin(), m_clusterVectorW.end());

    for (ClusterList::const_iterator iter = allClustersList.begin(), iterEnd = allClustersList.end(); iter != iterEnd; ++iter)
    {
        LArClusterHelper::TwoDSlidingFitResult slidingFitResult;
        LArClusterHelper::LArTwoDSlidingFit(*iter, 20, slidingFitResult);

        if (!m_slidingFitResultMap.insert(SlidingFitResultMap::value_type(*iter, slidingFitResult)).second)
            throw StatusCodeException(STATUS_CODE_FAILURE);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDTransverseTracksAlgorithm::CalculateOverlapResult(Cluster *pClusterU, Cluster *pClusterV, Cluster *pClusterW)
{
    SlidingFitResultMap::const_iterator iterU = m_slidingFitResultMap.find(pClusterU);
    SlidingFitResultMap::const_iterator iterV = m_slidingFitResultMap.find(pClusterV);
    SlidingFitResultMap::const_iterator iterW = m_slidingFitResultMap.find(pClusterW);

    if ((m_slidingFitResultMap.end() == iterU) || (m_slidingFitResultMap.end() == iterV) || (m_slidingFitResultMap.end() == iterW))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const LArClusterHelper::TwoDSlidingFitResult &slidingFitResultU(iterU->second);
    const LArClusterHelper::TwoDSlidingFitResult &slidingFitResultV(iterV->second);
    const LArClusterHelper::TwoDSlidingFitResult &slidingFitResultW(iterW->second);

    FitSegmentTensor fitSegmentTensor;
    this->GetFitSegmentTensor(slidingFitResultU, slidingFitResultV, slidingFitResultW, fitSegmentTensor);
    const TrackOverlapResult trackOverlapResult(this->GetBestOverlapResult(fitSegmentTensor));

    const int nLayersSpannedU(slidingFitResultU.GetMaxLayer() - slidingFitResultU.GetMinLayer());
    const int nLayersSpannedV(slidingFitResultV.GetMaxLayer() - slidingFitResultV.GetMinLayer());
    const int nLayersSpannedW(slidingFitResultW.GetMaxLayer() - slidingFitResultW.GetMinLayer());
    const unsigned int meanLayersSpanned(static_cast<unsigned int>((1.f / 3.f) * static_cast<float>(nLayersSpannedU + nLayersSpannedV + nLayersSpannedW)));

    if (0 == meanLayersSpanned)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const float nSamplingPointsPerLayer(static_cast<float>(trackOverlapResult.GetNSamplingPoints()) / static_cast<float>(meanLayersSpanned));

    if (nSamplingPointsPerLayer > m_minOverallMatchedFraction)
         m_overlapTensor.SetOverlapResult(pClusterU, pClusterV, pClusterW, trackOverlapResult);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDTransverseTracksAlgorithm::GetFitSegmentTensor(const TwoDSlidingFitResult &slidingFitResultU, const TwoDSlidingFitResult &slidingFitResultV,
    const TwoDSlidingFitResult &slidingFitResultW, FitSegmentTensor &fitSegmentTensor) const
{
    FitSegmentList fitSegmentListU, fitSegmentListV, fitSegmentListW;
    this->GetFitSegmentList(slidingFitResultU, fitSegmentListU);
    this->GetFitSegmentList(slidingFitResultV, fitSegmentListV);
    this->GetFitSegmentList(slidingFitResultW, fitSegmentListW);

    for (unsigned int indexU = 0, indexUEnd = fitSegmentListU.size(); indexU < indexUEnd; ++indexU)
    {
        const FitSegment &fitSegmentU(fitSegmentListU.at(indexU));

        for (unsigned int indexV = 0, indexVEnd = fitSegmentListV.size(); indexV < indexVEnd; ++indexV)
        {
            const FitSegment &fitSegmentV(fitSegmentListV.at(indexV));

            for (unsigned int indexW = 0, indexWEnd = fitSegmentListW.size(); indexW < indexWEnd; ++indexW)
            {
                const FitSegment &fitSegmentW(fitSegmentListW.at(indexW));

                try
                {
                    const TrackOverlapResult segmentOverlap(this->GetSegmentOverlap(fitSegmentU, fitSegmentV, fitSegmentW,
                        slidingFitResultU, slidingFitResultV, slidingFitResultW));

                    if ((segmentOverlap.GetMatchedFraction() < m_minSegmentMatchedFraction) || (segmentOverlap.GetNMatchedSamplingPoints() < m_minSegmentMatchedPoints)) // TODO
                        continue;

                    if (!fitSegmentTensor[indexU][indexV].insert(FitSegmentToOverlapResultMap::value_type(indexW, segmentOverlap)).second)
                        throw StatusCodeException(STATUS_CODE_FAILURE);
                }
                catch (StatusCodeException &statusCodeException)
                {
                    if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                        throw statusCodeException;
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDTransverseTracksAlgorithm::GetFitSegmentList(const LArClusterHelper::TwoDSlidingFitResult &slidingFitResult, FitSegmentList &fitSegmentList) const
{
    unsigned int nSustainedSteps(0);
    CartesianVector previousPosition(0.f, 0.f, 0.f);
    SlidingFitDirection previousDirection(UNKNOWN), sustainedDirection(UNKNOWN);
    LayerFitResultMap::const_iterator sustainedDirectionStartIter, sustainedDirectionEndIter;

    const LayerFitResultMap &layerFitResultMap(slidingFitResult.GetLayerFitResultMap());

    for (LayerFitResultMap::const_iterator iter = layerFitResultMap.begin(), iterEnd = layerFitResultMap.end(); iter != iterEnd; ++iter)
    {
        CartesianVector position(0.f, 0.f, 0.f);
        slidingFitResult.GetGlobalPosition(iter->second.GetL(), iter->second.GetFitT(), position);

        const CartesianVector delta(position - previousPosition);
        const SlidingFitDirection currentDirection((std::fabs(delta.GetX()) < std::fabs(delta.GetZ()) * -1.f) ?
            UNCHANGED_IN_X : (delta.GetX() > 0.f) ? POSITIVE_IN_X : NEGATIVE_IN_X);

        if (previousDirection == currentDirection)
        {
            if (++nSustainedSteps > 10)
            {
                sustainedDirection = currentDirection;
                sustainedDirectionEndIter = iter;
            }
        }
        else
        {
            if ((POSITIVE_IN_X == sustainedDirection) || (NEGATIVE_IN_X == sustainedDirection))
                fitSegmentList.push_back(FitSegment(slidingFitResult, sustainedDirectionStartIter, sustainedDirectionEndIter));

            nSustainedSteps = 0;
            sustainedDirection = UNKNOWN;
            sustainedDirectionStartIter = iter;
        }

        previousPosition = position;
        previousDirection = currentDirection;
    }

    if ((POSITIVE_IN_X == sustainedDirection) || (NEGATIVE_IN_X == sustainedDirection))
        fitSegmentList.push_back(FitSegment(slidingFitResult, sustainedDirectionStartIter, sustainedDirectionEndIter));
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackOverlapResult ThreeDTransverseTracksAlgorithm::GetSegmentOverlap(const FitSegment &fitSegmentU, const FitSegment &fitSegmentV, const FitSegment &fitSegmentW,
    const TwoDSlidingFitResult &slidingFitResultU, const TwoDSlidingFitResult &slidingFitResultV, const TwoDSlidingFitResult &slidingFitResultW) const
{
    // Assess x-overlap
    const float xSpanU(fitSegmentU.GetMaxX() - fitSegmentU.GetMinX());
    const float xSpanV(fitSegmentV.GetMaxX() - fitSegmentV.GetMinX());
    const float xSpanW(fitSegmentW.GetMaxX() - fitSegmentW.GetMinX());

    const float minX(std::max(fitSegmentU.GetMinX(), std::max(fitSegmentV.GetMinX(), fitSegmentW.GetMinX())));
    const float maxX(std::min(fitSegmentU.GetMaxX(), std::min(fitSegmentV.GetMaxX(), fitSegmentW.GetMaxX())));
    const float xOverlap(maxX - minX);

    if (xOverlap < 0.f)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    // Sampling in x
    const float nPointsU(std::fabs((xOverlap / xSpanU) * static_cast<float>(fitSegmentU.GetEndLayer() - fitSegmentU.GetStartLayer())));
    const float nPointsV(std::fabs((xOverlap / xSpanV) * static_cast<float>(fitSegmentV.GetEndLayer() - fitSegmentV.GetStartLayer())));
    const float nPointsW(std::fabs((xOverlap / xSpanW) * static_cast<float>(fitSegmentW.GetEndLayer() - fitSegmentW.GetStartLayer())));
    const float xPitch(3.f * xOverlap / (nPointsU + nPointsV + nPointsW));

    // Chi2 calculations
    float pseudoChi2Sum(0.f);
    unsigned int nSamplingPoints(0), nMatchedSamplingPoints(0);

    for (float x = minX; x < maxX; x += xPitch)
    {
        try
        {
            CartesianVector fitUVector(0.f, 0.f, 0.f), fitVVector(0.f, 0.f, 0.f), fitWVector(0.f, 0.f, 0.f);
            slidingFitResultU.GetGlobalFitPosition(x, true, fitUVector);
            slidingFitResultV.GetGlobalFitPosition(x, true, fitVVector);
            slidingFitResultW.GetGlobalFitPosition(x, true, fitWVector);

            CartesianVector fitUDirection(0.f, 0.f, 0.f), fitVDirection(0.f, 0.f, 0.f), fitWDirection(0.f, 0.f, 0.f);
            slidingFitResultU.GetGlobalFitDirection(x, true, fitUDirection);
            slidingFitResultV.GetGlobalFitDirection(x, true, fitVDirection);
            slidingFitResultW.GetGlobalFitDirection(x, true, fitWDirection);

            const float u(fitUVector.GetZ()), v(fitVVector.GetZ()), w(fitWVector.GetZ());
            const float uv2w(LArGeometryHelper::MergeTwoPositions(VIEW_U, VIEW_V, u, v));
            const float uw2v(LArGeometryHelper::MergeTwoPositions(VIEW_U, VIEW_W, u, w));
            const float vw2u(LArGeometryHelper::MergeTwoPositions(VIEW_V, VIEW_W, v, w));

            ++nSamplingPoints;
            const float deltaU((vw2u - u) * fitUDirection.GetX());
            const float deltaV((uw2v - v) * fitVDirection.GetX());
            const float deltaW((uv2w - w) * fitWDirection.GetX());

            const float pseudoChi2(deltaW * deltaW + deltaV * deltaV + deltaU * deltaU);
            pseudoChi2Sum += pseudoChi2;

            if (pseudoChi2 < m_pseudoChi2Cut)
                ++nMatchedSamplingPoints;
        }
        catch (StatusCodeException &)
        {
        }
    }

    if (0 == nSamplingPoints)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return TrackOverlapResult(nMatchedSamplingPoints, nSamplingPoints, pseudoChi2Sum);
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackOverlapResult ThreeDTransverseTracksAlgorithm::GetBestOverlapResult(const FitSegmentTensor &fitSegmentTensor) const
{
    // TODO, will need to navigate across fit segment tensor in multiple different directions
    const TrackOverlapResult noMatchResult(0, 1, 0.f);

    if (fitSegmentTensor.empty())
        return noMatchResult;

    TrackOverlapResultVector trackOverlapResultVector;
    unsigned int indexU(0), indexV(0), indexW(0);
    TrackOverlapResult trackOverlapResult(0, 1, 0.f);

    this->GetFirstMatch(fitSegmentTensor, indexU, indexV, indexW, trackOverlapResult);
    this->GetNeighbours(fitSegmentTensor, indexU, indexV, indexW, trackOverlapResult, trackOverlapResultVector);

    TrackOverlapResultVector::const_iterator maxElement = std::max_element(trackOverlapResultVector.begin(), trackOverlapResultVector.end());

    if (trackOverlapResultVector.end() == maxElement)
        return noMatchResult;

    return (*maxElement);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDTransverseTracksAlgorithm::GetFirstMatch(const FitSegmentTensor &fitSegmentTensor, unsigned int &indexU, unsigned int &indexV,
    unsigned int &indexW, TrackOverlapResult &trackOverlapResult) const
{
    if (fitSegmentTensor.empty())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    FitSegmentTensor::const_iterator iterU = fitSegmentTensor.begin();
    const FitSegmentMatrix &fitSegmentMatrix(iterU->second);
    FitSegmentMatrix::const_iterator iterV = fitSegmentMatrix.begin();

    if (fitSegmentMatrix.end() == iterV)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const FitSegmentToOverlapResultMap &fitSegmentToOverlapResultMap(iterV->second);
    FitSegmentToOverlapResultMap::const_iterator iterW = fitSegmentToOverlapResultMap.begin();

    if (fitSegmentToOverlapResultMap.end() == iterW)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    indexU = iterU->first;
    indexV = iterV->first;
    indexW = iterW->first;
    trackOverlapResult = iterW->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDTransverseTracksAlgorithm::GetNeighbours(const FitSegmentTensor &fitSegmentTensor, const unsigned int indexU, const unsigned int indexV,
    const unsigned int indexW, const TrackOverlapResult &trackOverlapResult, TrackOverlapResultVector &trackOverlapResultVector) const
{
    if (fitSegmentTensor.empty())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    bool neighbourFound(false);

    // Care with permutations, not interested in case where index is unchanged
    for (unsigned int iPermutation = 1; iPermutation < 8; ++iPermutation)
    {
        const bool incrementU((iPermutation >> 0) & 0x1);
        const bool incrementV((iPermutation >> 1) & 0x1);
        const bool incrementW((iPermutation >> 2) & 0x1);

        TrackOverlapResult additionalOverlapResult(0, 1, 0.f);
        unsigned int newIndexU(0), newIndexV(0), newIndexW(0);

        if (this->IsPresent(fitSegmentTensor, indexU, indexV, indexW, incrementU, incrementV, incrementW, newIndexU, newIndexV, newIndexW, additionalOverlapResult))
            this->GetNeighbours(fitSegmentTensor, newIndexU, newIndexV, newIndexW, trackOverlapResult + additionalOverlapResult, trackOverlapResultVector);

        if (additionalOverlapResult.GetNMatchedSamplingPoints() > 0)
            neighbourFound = true;
    }

    if (!neighbourFound)
        trackOverlapResultVector.push_back(trackOverlapResult);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ThreeDTransverseTracksAlgorithm::IsPresent(const FitSegmentTensor &fitSegmentTensor, const unsigned int indexU, const unsigned int indexV,
    const unsigned int indexW, const bool incrementU, const bool incrementV, const bool incrementW, unsigned int &newIndexU,
    unsigned int &newIndexV, unsigned int &newIndexW, TrackOverlapResult &trackOverlapResult) const
{
    FitSegmentTensor::const_iterator iterU = fitSegmentTensor.find(indexU);

    if ((fitSegmentTensor.end() == iterU) || (incrementU && (fitSegmentTensor.end() == ++iterU)))
        return false;

    const FitSegmentMatrix &fitSegmentMatrix(iterU->second);
    FitSegmentMatrix::const_iterator iterV = fitSegmentMatrix.find(indexV);

    if ((fitSegmentMatrix.end() == iterV) || (incrementV && ((fitSegmentMatrix.end() == ++iterV))))
        return false;

    const FitSegmentToOverlapResultMap &fitSegmentToOverlapResultMap(iterV->second);
    FitSegmentToOverlapResultMap::const_iterator iterW = fitSegmentToOverlapResultMap.find(indexW);

    if ((fitSegmentToOverlapResultMap.end() == iterW) || (incrementW && ((fitSegmentToOverlapResultMap.end() == ++iterW))))
        return false;

    newIndexU = iterU->first;
    newIndexV = iterV->first;
    newIndexW = iterW->first;
    trackOverlapResult = iterW->second;
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ThreeDTransverseTracksAlgorithm::ExamineTensor()
{
    Cluster *pBestClusterU(NULL), *pBestClusterV(NULL), *pBestClusterW(NULL);
    TrackOverlapResult bestTrackOverlapResult(0, 1, 0.f);

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
                    const TrackOverlapResult &trackOverlapResult(m_overlapTensor.GetOverlapResult(*iterU, *iterV, *iterW));

                    if (trackOverlapResult < bestTrackOverlapResult)
                        continue;

                    bestTrackOverlapResult = trackOverlapResult;
                    pBestClusterU = *iterU;
                    pBestClusterV = *iterV;
                    pBestClusterW = *iterW;
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
    this->BuildProtoParticle(ParticleComponent(pBestClusterU, pBestClusterV, pBestClusterW, bestTrackOverlapResult), protoParticle);
    m_protoParticleVector.push_back(protoParticle);

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDTransverseTracksAlgorithm::BuildProtoParticle(const ParticleComponent &particleComponent, ProtoParticle &protoParticle) const
{
    Cluster *pClusterU(particleComponent.GetClusterU()), *pClusterV(particleComponent.GetClusterV()), *pClusterW(particleComponent.GetClusterW());
    protoParticle.m_clusterListU.insert(pClusterU);
    protoParticle.m_clusterListV.insert(pClusterV);
    protoParticle.m_clusterListW.insert(pClusterW);

    ParticleComponentList particleComponentList;
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
                    if ((pClusterU != *iterU) && (pClusterV != *iterV) && (pClusterW != *iterW))
                        continue;

                    const TrackOverlapResult &trackOverlapResult(m_overlapTensor.GetOverlapResult(*iterU, *iterV, *iterW));
                    particleComponentList.push_back(ParticleComponent(*iterU, *iterV, *iterW, trackOverlapResult));
                }
                catch (StatusCodeException &)
                {
                }
            }
        }
    }

    for (ParticleComponentList::const_iterator iter = particleComponentList.begin(), iterEnd = particleComponentList.end(); iter != iterEnd; ++iter)
    {
        if (this->IsParticleMatch(*iter, protoParticle))
            this->BuildProtoParticle(*iter, protoParticle);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ThreeDTransverseTracksAlgorithm::IsParticleMatch(const ParticleComponent &particleComponent, const ProtoParticle &protoParticle) const
{
    // TODO, return true in case where we have a particle match
    std::cout << "CHECK IS PARTICLE MATCH " << std::endl;
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDTransverseTracksAlgorithm::TidyUp()
{
    m_slidingFitResultMap.clear();
    return ThreeDBaseAlgorithm<TrackOverlapResult>::TidyUp();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDTransverseTracksAlgorithm::FitSegment::FitSegment(const LArClusterHelper::TwoDSlidingFitResult &twoDSlidingFitResult,
        LayerFitResultMap::const_iterator startLayerIter, LayerFitResultMap::const_iterator endLayerIter) :
    m_startLayer(startLayerIter->first),
    m_endLayer(endLayerIter->first)
{
    CartesianVector startLayerPosition(0.f, 0.f, 0.f);
    twoDSlidingFitResult.GetGlobalPosition(startLayerIter->second.GetL(), startLayerIter->second.GetFitT(), startLayerPosition);

    CartesianVector endLayerPosition(0.f, 0.f, 0.f);
    twoDSlidingFitResult.GetGlobalPosition(endLayerIter->second.GetL(), endLayerIter->second.GetFitT(), endLayerPosition);

    m_minX = std::min(startLayerPosition.GetX(), endLayerPosition.GetX());
    m_maxX = std::max(startLayerPosition.GetX(), endLayerPosition.GetX());
    m_startValue = startLayerPosition.GetZ();
    m_endValue = endLayerPosition.GetZ();
    m_isIncreasingX = (endLayerPosition.GetX() > startLayerPosition.GetX());
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDTransverseTracksAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_pseudoChi2Cut = 3.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PseudoChi2Cut", m_pseudoChi2Cut));

    m_minOverallMatchedFraction = 0.2f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinOverallMatchedFraction", m_minOverallMatchedFraction));

    m_minSegmentMatchedFraction = 0.1f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinSegmentMatchedFraction", m_minSegmentMatchedFraction));

    m_minSegmentMatchedPoints = 3;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinSegmentMatchedPoints", m_minSegmentMatchedPoints));

    return ThreeDBaseAlgorithm<TrackOverlapResult>::ReadSettings(xmlHandle);
}

} // namespace lar
