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
#include "LArHelpers/LArVertexHelper.h"

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

                    if ((segmentOverlap.GetMatchedFraction() < m_minSegmentMatchedFraction) || (segmentOverlap.GetNMatchedSamplingPoints() < m_minSegmentMatchedPoints))
                        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

                    fitSegmentTensor[indexU][indexV][indexW] = segmentOverlap;
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

TrackOverlapResult ThreeDTransverseTracksAlgorithm::GetBestOverlapResult(FitSegmentTensor &fitSegmentTensor) const
{
    // TODO - try to work with tensor that a) is const and b) contains only the non-zero entries
    const unsigned int maxIndexU(fitSegmentTensor.size());
    const unsigned int maxIndexV(fitSegmentTensor[0].size());
    const unsigned int maxIndexW(fitSegmentTensor[0][0].size());

    if ((0 == maxIndexU) || (0 == maxIndexV) || (0 == maxIndexW))
        return TrackOverlapResult();

    FitSegmentTensor fitSegmentSumTensor;
    for (unsigned int indexU = 0; indexU < maxIndexU; ++indexU)
    {
        for (unsigned int indexV = 0; indexV < maxIndexV; ++indexV)
        {
            for (unsigned int indexW = 0; indexW < maxIndexW; ++indexW)
            {
                TrackOverlapResultVector trackOverlapResultVector;
                this->GetPreviousOverlapResults(indexU, indexV, indexW, maxIndexU, maxIndexV, maxIndexW, fitSegmentSumTensor, trackOverlapResultVector);
                TrackOverlapResultVector::const_iterator maxElement = std::max_element(trackOverlapResultVector.begin(), trackOverlapResultVector.end());
                TrackOverlapResult maxTrackOverlapResult((trackOverlapResultVector.end() != maxElement) ? *maxElement : TrackOverlapResult());
                fitSegmentSumTensor[indexU][indexV][indexW] = maxTrackOverlapResult + fitSegmentTensor[indexU][indexV][indexW];
            }
        }
    }

    TrackOverlapResult bestTrackOverlapResult;
    for (unsigned int indexU = 0; indexU < maxIndexU; ++indexU)
    {
        for (unsigned int indexV = 0; indexV < maxIndexV; ++indexV)
        {
            for (unsigned int indexW = 0; indexW < maxIndexW; ++indexW)
            {
                const TrackOverlapResult &trackOverlapResult(fitSegmentSumTensor[indexU][indexV][indexW]);

                if (trackOverlapResult > bestTrackOverlapResult)
                    bestTrackOverlapResult = trackOverlapResult;
            }
        }
    }

    return bestTrackOverlapResult;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDTransverseTracksAlgorithm::GetPreviousOverlapResults(const unsigned int indexU, const unsigned int indexV, const unsigned int indexW,
    const unsigned int maxIndexU, const unsigned int maxIndexV, const unsigned int maxIndexW, FitSegmentTensor &fitSegmentSumTensor,
    TrackOverlapResultVector &trackOverlapResultVector) const
{
    for (unsigned int iPermutation = 1; iPermutation < 8; ++iPermutation)
    {
        const bool decrementU((iPermutation >> 0) & 0x1);
        const bool decrementV((iPermutation >> 1) & 0x1);
        const bool decrementW((iPermutation >> 2) & 0x1);

        if ((decrementU && (0 == indexU)) || (decrementV && (0 == indexV)) || (decrementW && (0 == indexW)))
            continue;

        const unsigned int newIndexU(decrementU ? indexU - 1 : indexU);
        const unsigned int newIndexV(decrementV ? indexV - 1 : indexV);
        const unsigned int newIndexW(decrementW ? indexW - 1 : indexW);
        trackOverlapResultVector.push_back(fitSegmentSumTensor[newIndexU][newIndexV][newIndexW]);
    }

    if (trackOverlapResultVector.empty())
        trackOverlapResultVector.push_back(TrackOverlapResult());
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

void ThreeDTransverseTracksAlgorithm::BuildProtoParticle(const ParticleComponent &firstComponent, ProtoParticle &protoParticle) const
{
    Cluster *pClusterU(firstComponent.GetClusterU()), *pClusterV(firstComponent.GetClusterV()), *pClusterW(firstComponent.GetClusterW());
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
        if (protoParticle.m_clusterListU.count(iter->GetClusterU()) && protoParticle.m_clusterListV.count(iter->GetClusterV()) &&
            protoParticle.m_clusterListW.count(iter->GetClusterW()))
        {
            continue;
        }

        if (this->IsParticleMatch(firstComponent, *iter))
        {
            this->BuildProtoParticle(*iter, protoParticle);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ThreeDTransverseTracksAlgorithm::IsParticleMatch(const ParticleComponent &firstComponent, const ParticleComponent &secondComponent) const
{
    return ((this->IsPossibleMatch(firstComponent.GetClusterU(), secondComponent.GetClusterU()) &&
             this->IsPossibleMatch(firstComponent.GetClusterV(), secondComponent.GetClusterV()) &&
             this->IsPossibleMatch(firstComponent.GetClusterW(), secondComponent.GetClusterW())) &&
            (this->IsParticleMatch(firstComponent.GetClusterU(), secondComponent.GetClusterU()) ||
             this->IsParticleMatch(firstComponent.GetClusterV(), secondComponent.GetClusterV()) ||
             this->IsParticleMatch(firstComponent.GetClusterW(), secondComponent.GetClusterW())));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ThreeDTransverseTracksAlgorithm::IsPossibleMatch(Cluster *const pFirstCluster, Cluster *const pSecondCluster) const
{
    if (pFirstCluster == pSecondCluster)
        return true;

    SlidingFitResultMap::const_iterator iter1 = m_slidingFitResultMap.find(pFirstCluster);
    SlidingFitResultMap::const_iterator iter2 = m_slidingFitResultMap.find(pSecondCluster);

    if ((m_slidingFitResultMap.end() == iter1) || (m_slidingFitResultMap.end() == iter2))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const LArClusterHelper::TwoDSlidingFitResult &slidingFitResult1(iter1->second);
    const LArClusterHelper::TwoDSlidingFitResult &slidingFitResult2(iter2->second);

    // Check there is no significant overlap between clusters
    const CartesianVector minLayerPosition1(slidingFitResult1.GetGlobalMinLayerPosition());
    const CartesianVector maxLayerPosition1(slidingFitResult1.GetGlobalMaxLayerPosition());
    const CartesianVector minLayerPosition2(slidingFitResult2.GetGlobalMinLayerPosition());
    const CartesianVector maxLayerPosition2(slidingFitResult2.GetGlobalMaxLayerPosition());

    const CartesianVector linearDirection1((maxLayerPosition1 - minLayerPosition1).GetUnitVector());
    const CartesianVector linearDirection2((maxLayerPosition2 - minLayerPosition2).GetUnitVector());

    const float clusterOverlap1_Pos0((maxLayerPosition1 - minLayerPosition1).GetMagnitude());
    const float clusterOverlap1_Pos1(linearDirection1.GetDotProduct(minLayerPosition2 - minLayerPosition1));
    const float clusterOverlap1_Pos2(linearDirection1.GetDotProduct(maxLayerPosition2 - minLayerPosition1));

    const float clusterOverlap2_Pos0((maxLayerPosition2 - minLayerPosition2).GetMagnitude());
    const float clusterOverlap2_Pos1(linearDirection2.GetDotProduct(minLayerPosition1 - minLayerPosition2));
    const float clusterOverlap2_Pos2(linearDirection2.GetDotProduct(maxLayerPosition1 - minLayerPosition2));

    return ((std::min(clusterOverlap1_Pos1, clusterOverlap1_Pos2) > clusterOverlap1_Pos0 - 2.f) ||
            (std::min(clusterOverlap2_Pos1, clusterOverlap2_Pos2) > clusterOverlap2_Pos0 - 2.f) ||
            (std::max(clusterOverlap1_Pos1, clusterOverlap1_Pos2) < 2.f) ||
            (std::max(clusterOverlap2_Pos1, clusterOverlap2_Pos2) < 2.f));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ThreeDTransverseTracksAlgorithm::IsParticleMatch(Cluster *const pFirstCluster, Cluster *const pSecondCluster) const
{
    if (pFirstCluster == pSecondCluster)
        return false;

    if (!LArVertexHelper::DoesCurrentVertexExist())
        return false;

    SlidingFitResultMap::const_iterator iter1 = m_slidingFitResultMap.find(pFirstCluster);
    SlidingFitResultMap::const_iterator iter2 = m_slidingFitResultMap.find(pSecondCluster);

    if ((m_slidingFitResultMap.end() == iter1) || (m_slidingFitResultMap.end() == iter2))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const LArClusterHelper::TwoDSlidingFitResult &slidingFitResult1(iter1->second);
    const LArClusterHelper::TwoDSlidingFitResult &slidingFitResult2(iter2->second);

    const CartesianVector minLayerPosition1(slidingFitResult1.GetGlobalMinLayerPosition());
    const CartesianVector maxLayerPosition1(slidingFitResult1.GetGlobalMaxLayerPosition());
    const CartesianVector minLayerPosition2(slidingFitResult2.GetGlobalMinLayerPosition());
    const CartesianVector maxLayerPosition2(slidingFitResult2.GetGlobalMaxLayerPosition());

    const bool isForward1(LArVertexHelper::IsForwardInZ3D(pFirstCluster));
    const bool isForward2(LArVertexHelper::IsForwardInZ3D(pSecondCluster));
    const bool isBackward1(LArVertexHelper::IsBackwardInZ3D(pFirstCluster));
    const bool isBackward2(LArVertexHelper::IsBackwardInZ3D(pSecondCluster));

    return ((((minLayerPosition1 - minLayerPosition2).GetMagnitudeSquared() < 4.f) && !(isForward1 && isForward2)  && !(isBackward1 && isBackward2)) ||
            (((maxLayerPosition1 - maxLayerPosition2).GetMagnitudeSquared() < 4.f) && !(isForward1 && isForward2)  && !(isBackward1 && isBackward2)) ||
            (((minLayerPosition1 - maxLayerPosition2).GetMagnitudeSquared() < 4.f) && !(isForward1 && isBackward2) && !(isBackward1 && isForward2)) ||
            (((maxLayerPosition1 - minLayerPosition2).GetMagnitudeSquared() < 4.f) && !(isForward1 && isBackward2) && !(isBackward1 && isForward2)) );
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
