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

void ThreeDTransverseTracksAlgorithm::CalculateOverlapResult(Cluster *pClusterU, Cluster *pClusterV, Cluster *pClusterW)
{
std::cout << " CALCULATE OVERLAP RESULT " << std::endl;
PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
ClusterList clusterListU; clusterListU.insert(pClusterU);
PandoraMonitoringApi::VisualizeClusters(&clusterListU, "ClusterListU", RED);
ClusterList clusterListV; clusterListV.insert(pClusterV);
PandoraMonitoringApi::VisualizeClusters(&clusterListV, "ClusterListV", GREEN);
ClusterList clusterListW; clusterListW.insert(pClusterW);
PandoraMonitoringApi::VisualizeClusters(&clusterListW, "ClusterListW", BLUE);
PandoraMonitoringApi::ViewEvent();
    LArClusterHelper::TwoDSlidingFitResult slidingFitResultU, slidingFitResultV, slidingFitResultW;
    LArClusterHelper::LArTwoDSlidingFit(pClusterU, 20, slidingFitResultU);
    LArClusterHelper::LArTwoDSlidingFit(pClusterV, 20, slidingFitResultV);
    LArClusterHelper::LArTwoDSlidingFit(pClusterW, 20, slidingFitResultW);

    FitSegmentTensor fitSegmentTensor;
    this->GetFitSegmentTensor(slidingFitResultU, slidingFitResultV, slidingFitResultW, fitSegmentTensor);
std::cout << " fitSegmentTensor.size() " << fitSegmentTensor.size() << std::endl;
    const TrackOverlapResult trackOverlapResult(this->GetBestOverlapResult(fitSegmentTensor));

    // const int nLayersSpannedU(slidingFitResultU.GetMaxLayer() - slidingFitResultU.GetMinLayer());
    // const int nLayersSpannedV(slidingFitResultV.GetMaxLayer() - slidingFitResultV.GetMinLayer());
    // const int nLayersSpannedW(slidingFitResultW.GetMaxLayer() - slidingFitResultW.GetMinLayer());
    // const unsigned int nSamplingPoints(static_cast<unsigned int>((1.f / 3.f) * static_cast<float>(nLayersSpannedU + nLayersSpannedV + nLayersSpannedW)));

    // This is the overlap result that actually goes into the final tensor
    if (trackOverlapResult.GetNMatchedSamplingPoints() > 0)
         m_overlapTensor.SetOverlapResult(pClusterU, pClusterV, pClusterW, trackOverlapResult);

std::cout << " FINISHED CLUSTER COMBINATION SetOverlapResult, nMatcheds " << trackOverlapResult.GetNMatchedSamplingPoints() << " nSamples " << trackOverlapResult.GetNSamplingPoints() << " fraction " << trackOverlapResult.GetMatchedFraction() << " chi2 " << trackOverlapResult.GetChi2() << " reducedchi2 " << trackOverlapResult.GetReducedChi2() << std::endl;
PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
PandoraMonitoringApi::VisualizeClusters(&clusterListU, "ClusterListU", RED);
PandoraMonitoringApi::VisualizeClusters(&clusterListV, "ClusterListV", GREEN);
PandoraMonitoringApi::VisualizeClusters(&clusterListW, "ClusterListW", BLUE);
PandoraMonitoringApi::ViewEvent();
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
std::cout << " Segment overlap " << segmentOverlap.GetMatchedFraction() << ", indices, indexU " << indexU << " indexV " << indexV << " indexW " << indexW << std::endl;
                    if ((segmentOverlap.GetMatchedFraction() < 0.1f) || (segmentOverlap.GetNMatchedSamplingPoints() < 5)) // TODO
{std::cout << " segmentOverlap discarded matchedFraction " << segmentOverlap.GetMatchedFraction() << " matchedPoints " << segmentOverlap.GetNMatchedSamplingPoints() << std::endl;
PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
LayerFitResultMap::const_iterator startIter = slidingFitResultU.GetLayerFitResultMap().find(fitSegmentU.GetStartLayer());
LayerFitResultMap::const_iterator endIter = slidingFitResultU.GetLayerFitResultMap().find(fitSegmentU.GetEndLayer());
if ((slidingFitResultU.GetLayerFitResultMap().end() == startIter) || (slidingFitResultU.GetLayerFitResultMap().end() == endIter))
    throw;
for (LayerFitResultMap::const_iterator mIter = startIter; !(mIter == endIter); ++mIter)
{
    CartesianVector position(0.f, 0.f, 0.f);
    slidingFitResultU.GetGlobalPosition(mIter->second.GetL(), mIter->second.GetFitT(), position);
    PANDORA_MONITORING_API(AddMarkerToVisualization(&position, "positionU", RED, 1));
}
startIter = slidingFitResultV.GetLayerFitResultMap().find(fitSegmentV.GetStartLayer());
endIter = slidingFitResultV.GetLayerFitResultMap().find(fitSegmentV.GetEndLayer());
if ((slidingFitResultV.GetLayerFitResultMap().end() == startIter) || (slidingFitResultV.GetLayerFitResultMap().end() == endIter))
    throw;
for (LayerFitResultMap::const_iterator mIter = startIter; !(mIter == endIter); ++mIter)
{
    CartesianVector position(0.f, 0.f, 0.f);
    slidingFitResultV.GetGlobalPosition(mIter->second.GetL(), mIter->second.GetFitT(), position);
    PANDORA_MONITORING_API(AddMarkerToVisualization(&position, "positionV", GREEN, 1));
}
startIter = slidingFitResultW.GetLayerFitResultMap().find(fitSegmentW.GetStartLayer());
endIter = slidingFitResultW.GetLayerFitResultMap().find(fitSegmentW.GetEndLayer());
if ((slidingFitResultW.GetLayerFitResultMap().end() == startIter) || (slidingFitResultW.GetLayerFitResultMap().end() == endIter))
    throw;
for (LayerFitResultMap::const_iterator mIter = startIter; !(mIter == endIter); ++mIter)
{
    CartesianVector position(0.f, 0.f, 0.f);
    slidingFitResultW.GetGlobalPosition(mIter->second.GetL(), mIter->second.GetFitT(), position);
    PANDORA_MONITORING_API(AddMarkerToVisualization(&position, "positionW", BLUE, 1));
}
PandoraMonitoringApi::ViewEvent();

continue;
}

                    if (!fitSegmentTensor[indexU][indexV].insert(FitSegmentToOverlapResultMap::value_type(indexW, segmentOverlap)).second)
                        throw StatusCodeException(STATUS_CODE_FAILURE);
// DEBUG
std::cout << " PUT IN FIT SEGMENT TENSOR indexU " << indexU << " indexV " << indexV << " indexW " << indexW << std::endl;
std::cout << " nMatchedSamplingPoints " << segmentOverlap.GetNMatchedSamplingPoints() << " nSamplingPoints " << segmentOverlap.GetNSamplingPoints() << " fraction " << segmentOverlap.GetMatchedFraction() << std::endl;
PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
LayerFitResultMap::const_iterator startIter = slidingFitResultU.GetLayerFitResultMap().find(fitSegmentU.GetStartLayer());
LayerFitResultMap::const_iterator endIter = slidingFitResultU.GetLayerFitResultMap().find(fitSegmentU.GetEndLayer());
if ((slidingFitResultU.GetLayerFitResultMap().end() == startIter) || (slidingFitResultU.GetLayerFitResultMap().end() == endIter))
    throw;
for (LayerFitResultMap::const_iterator mIter = startIter; !(mIter == endIter); ++mIter)
{
    CartesianVector position(0.f, 0.f, 0.f);
    slidingFitResultU.GetGlobalPosition(mIter->second.GetL(), mIter->second.GetFitT(), position);
    PANDORA_MONITORING_API(AddMarkerToVisualization(&position, "positionU", RED, 1));
}
startIter = slidingFitResultV.GetLayerFitResultMap().find(fitSegmentV.GetStartLayer());
endIter = slidingFitResultV.GetLayerFitResultMap().find(fitSegmentV.GetEndLayer());
if ((slidingFitResultV.GetLayerFitResultMap().end() == startIter) || (slidingFitResultV.GetLayerFitResultMap().end() == endIter))
    throw;
for (LayerFitResultMap::const_iterator mIter = startIter; !(mIter == endIter); ++mIter)
{
    CartesianVector position(0.f, 0.f, 0.f);
    slidingFitResultV.GetGlobalPosition(mIter->second.GetL(), mIter->second.GetFitT(), position);
    PANDORA_MONITORING_API(AddMarkerToVisualization(&position, "positionV", GREEN, 1));
}
startIter = slidingFitResultW.GetLayerFitResultMap().find(fitSegmentW.GetStartLayer());
endIter = slidingFitResultW.GetLayerFitResultMap().find(fitSegmentW.GetEndLayer());
if ((slidingFitResultW.GetLayerFitResultMap().end() == startIter) || (slidingFitResultW.GetLayerFitResultMap().end() == endIter))
    throw;
for (LayerFitResultMap::const_iterator mIter = startIter; !(mIter == endIter); ++mIter)
{
    CartesianVector position(0.f, 0.f, 0.f);
    slidingFitResultW.GetGlobalPosition(mIter->second.GetL(), mIter->second.GetFitT(), position);
    PANDORA_MONITORING_API(AddMarkerToVisualization(&position, "positionW", BLUE, 1));
}
PandoraMonitoringApi::ViewEvent();
                }
                catch (StatusCodeException &statusCodeException)
                {
std::cout << " Segment overlap, exception " << statusCodeException.ToString() << ", indices: indexU " << indexU << " indexV " << indexV << " indexW " << indexW << std::endl;
PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
LayerFitResultMap::const_iterator startIter = slidingFitResultU.GetLayerFitResultMap().find(fitSegmentU.GetStartLayer());
LayerFitResultMap::const_iterator endIter = slidingFitResultU.GetLayerFitResultMap().find(fitSegmentU.GetEndLayer());
if ((slidingFitResultU.GetLayerFitResultMap().end() == startIter) || (slidingFitResultU.GetLayerFitResultMap().end() == endIter))
    throw;
for (LayerFitResultMap::const_iterator mIter = startIter; !(mIter == endIter); ++mIter)
{
    CartesianVector position(0.f, 0.f, 0.f);
    slidingFitResultU.GetGlobalPosition(mIter->second.GetL(), mIter->second.GetFitT(), position);
    PANDORA_MONITORING_API(AddMarkerToVisualization(&position, "positionU", RED, 1));
}
startIter = slidingFitResultV.GetLayerFitResultMap().find(fitSegmentV.GetStartLayer());
endIter = slidingFitResultV.GetLayerFitResultMap().find(fitSegmentV.GetEndLayer());
if ((slidingFitResultV.GetLayerFitResultMap().end() == startIter) || (slidingFitResultV.GetLayerFitResultMap().end() == endIter))
    throw;
for (LayerFitResultMap::const_iterator mIter = startIter; !(mIter == endIter); ++mIter)
{
    CartesianVector position(0.f, 0.f, 0.f);
    slidingFitResultV.GetGlobalPosition(mIter->second.GetL(), mIter->second.GetFitT(), position);
    PANDORA_MONITORING_API(AddMarkerToVisualization(&position, "positionV", GREEN, 1));
}
startIter = slidingFitResultW.GetLayerFitResultMap().find(fitSegmentW.GetStartLayer());
endIter = slidingFitResultW.GetLayerFitResultMap().find(fitSegmentW.GetEndLayer());
if ((slidingFitResultW.GetLayerFitResultMap().end() == startIter) || (slidingFitResultW.GetLayerFitResultMap().end() == endIter))
    throw;
for (LayerFitResultMap::const_iterator mIter = startIter; !(mIter == endIter); ++mIter)
{
    CartesianVector position(0.f, 0.f, 0.f);
    slidingFitResultW.GetGlobalPosition(mIter->second.GetL(), mIter->second.GetFitT(), position);
    PANDORA_MONITORING_API(AddMarkerToVisualization(&position, "positionW", BLUE, 1));
}
PandoraMonitoringApi::ViewEvent();

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
{   fitSegmentList.push_back(FitSegment(slidingFitResult, sustainedDirectionStartIter, sustainedDirectionEndIter));
//DEBUG display
sustainedDirectionEndIter = --(layerFitResultMap.end());
    }

std::cout << " GetFitSegmentList, fitSegmentList.size() " << fitSegmentList.size() << std::endl;
PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
int counter(0);

for (FitSegmentList::const_iterator iter = fitSegmentList.begin(), iterEnd = fitSegmentList.end(); iter != iterEnd; ++iter)
{
    const int startLayer(iter->GetStartLayer());
    const int endLayer(iter->GetEndLayer());

    LayerFitResultMap::const_iterator startIter = layerFitResultMap.find(startLayer);
    LayerFitResultMap::const_iterator endIter = layerFitResultMap.find(endLayer);

    if ((layerFitResultMap.end() == startIter) || (layerFitResultMap.end() == endIter))
        continue;

    ++counter;

    for (LayerFitResultMap::const_iterator mIter = startIter; !(mIter == endIter); ++mIter)
    {
        CartesianVector position(0.f, 0.f, 0.f);
        slidingFitResult.GetGlobalPosition(mIter->second.GetL(), mIter->second.GetFitT(), position);
        
        if (counter % 2 == 0)
        {
            PANDORA_MONITORING_API(AddMarkerToVisualization(&position, "position", RED, 2));
        }
        else
        {
            PANDORA_MONITORING_API(AddMarkerToVisualization(&position, "position", BLUE, 2));
        }
    }
}

for (LayerFitResultMap::const_iterator iter = layerFitResultMap.begin(), iterEnd = layerFitResultMap.end(); iter != iterEnd; ++iter)
{
    CartesianVector position(0.f, 0.f, 0.f);
    slidingFitResult.GetGlobalPosition(iter->second.GetL(), iter->second.GetFitT(), position);
    PANDORA_MONITORING_API(AddMarkerToVisualization(&position, "cluster", GREEN, 1));
}

PandoraMonitoringApi::ViewEvent();
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackOverlapResult ThreeDTransverseTracksAlgorithm::GetSegmentOverlap(const FitSegment &fitSegmentU, const FitSegment &fitSegmentV, const FitSegment &fitSegmentW,
    const TwoDSlidingFitResult &slidingFitResultU, const TwoDSlidingFitResult &slidingFitResultV, const TwoDSlidingFitResult &slidingFitResultW) const
{
    // Assess x-overlap
    const float xSpanU(fitSegmentU.GetMaxX() - fitSegmentU.GetMinX());
    const float xSpanV(fitSegmentV.GetMaxX() - fitSegmentV.GetMinX());
    const float xSpanW(fitSegmentW.GetMaxX() - fitSegmentW.GetMinX());
std::cout << " minXU " << fitSegmentU.GetMinX() << " maxXU " << fitSegmentU.GetMaxX() << std::endl;
std::cout << " minXV " << fitSegmentV.GetMinX() << " maxXV " << fitSegmentV.GetMaxX() << std::endl;
std::cout << " minXW " << fitSegmentW.GetMinX() << " maxXW " << fitSegmentW.GetMaxX() << std::endl;
    const float minX(std::max(fitSegmentU.GetMinX(), std::max(fitSegmentV.GetMinX(), fitSegmentW.GetMinX())));
    const float maxX(std::min(fitSegmentU.GetMaxX(), std::min(fitSegmentV.GetMaxX(), fitSegmentW.GetMaxX())));
    const float xOverlap(maxX - minX);
std::cout << " minX " << minX << " maxX " << maxX << " xOverlap " << xOverlap << std::endl;
    if ((xOverlap < 0.f))// TODO! || ((xOverlap / xSpanU) < 0.3f) || ((xOverlap / xSpanV) < 0.3f) || ((xOverlap / xSpanW) < 0.3f))
        {std::cout << "GetSegmentOverlap:: xOverlap " << xOverlap << std::endl; throw StatusCodeException(STATUS_CODE_NOT_FOUND);}

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
std::cout << " pseudoChi2 " << pseudoChi2 << " m_pseudoChi2Cut " << m_pseudoChi2Cut << " pseudoChi2Sum " << pseudoChi2Sum << std::endl;
const CartesianVector expU(x, 0., vw2u); PANDORA_MONITORING_API(AddMarkerToVisualization(&expU, "expU", RED, 2));
const CartesianVector expV(x, 0., uw2v); PANDORA_MONITORING_API(AddMarkerToVisualization(&expV, "expV", GREEN, 2));
const CartesianVector expW(x, 0., uv2w); PANDORA_MONITORING_API(AddMarkerToVisualization(&expW, "expW", BLUE, 2));
        }
        catch (StatusCodeException &)
        {
        }
    }

    if (0 == nSamplingPoints)
        {std::cout << "GetSegmentOverlap:: nSamplingPoints == 0 " << std::endl; throw StatusCodeException(STATUS_CODE_NOT_FOUND);}

    return TrackOverlapResult(nMatchedSamplingPoints, nSamplingPoints, pseudoChi2Sum);
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackOverlapResult ThreeDTransverseTracksAlgorithm::GetBestOverlapResult(const FitSegmentTensor &fitSegmentTensor) const
{
    const TrackOverlapResult noMatchResult(0, 1, 0.f);

    if (fitSegmentTensor.empty())
        return noMatchResult;

    TrackOverlapResultVector trackOverlapResultVector;
    unsigned int indexU(0), indexV(0), indexW(0);
    TrackOverlapResult trackOverlapResult(0, 1, 0.f);

    this->GetFirstMatch(fitSegmentTensor, indexU, indexV, indexW, trackOverlapResult);
std::cout << " First match " << indexU << ", " << indexV << ", " << indexW << ", nMatchedPoints " << trackOverlapResult.GetNMatchedSamplingPoints() << std::endl;
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
static unsigned int instanceCounter(0);
std::cout << " GetNeighbours instance " << ++instanceCounter << " indices: " << indexU << ", " << indexV << ", " << indexW << std::endl;
    bool neighbourFound(false);

    // Care with permutations, not interested in case where index is unchanged
    for (unsigned int iPermutation = 1; iPermutation < 8; ++iPermutation)
    {
        std::cout << " iPermutation " << iPermutation << std::endl;
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
--instanceCounter;
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
std::cout << " IsPresent iU " << indexU << " iV " << indexV << " iW " << indexW << " inc " << incrementU << ", " << incrementV << ", " << incrementW
<< " nIU " << newIndexU << " nIV " << newIndexV << " nIW " << newIndexW << " nMatch " << iterW->second.GetNMatchedSamplingPoints() << " nSample " << iterW->second.GetNSamplingPoints() << " frac " << iterW->second.GetMatchedFraction() << std::endl;
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ThreeDTransverseTracksAlgorithm::ExamineTensor()
{
    float bestOverlapResult(m_minMatchedFraction);
    Cluster *pBestClusterU(NULL), *pBestClusterV(NULL), *pBestClusterW(NULL);

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
                    const TrackOverlapResult &overlapResult(m_overlapTensor.GetOverlapResult(*iterU, *iterV, *iterW));

                    if (overlapResult.GetNMatchedSamplingPoints() < m_minMatchedPoints)
                        continue;

                    if (overlapResult.GetMatchedFraction() > bestOverlapResult)
                    {
                        bestOverlapResult = overlapResult.GetMatchedFraction();
                        pBestClusterU = *iterU;
                        pBestClusterV = *iterV;
                        pBestClusterW = *iterW;
                    }
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
    protoParticle.m_clusterVectorU.push_back(pBestClusterU);
    protoParticle.m_clusterVectorV.push_back(pBestClusterV);
    protoParticle.m_clusterVectorW.push_back(pBestClusterW);
    m_protoParticleVector.push_back(protoParticle);

    return true;
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

    m_minMatchedFraction = 0.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedFraction", m_minMatchedFraction));

    m_minMatchedPoints = 0;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedPoints", m_minMatchedPoints));

    return ThreeDBaseAlgorithm<TrackOverlapResult>::ReadSettings(xmlHandle);
}

} // namespace lar
