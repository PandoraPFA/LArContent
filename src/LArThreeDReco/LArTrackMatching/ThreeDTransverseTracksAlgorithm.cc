/**
 *  @file   LArContent/src/LArThreeDReco/LArTrackMatching/ThreeDTransverseTracksAlgorithm.cc
 * 
 *  @brief  Implementation of the three dimensional transverse tracks algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArGeometryHelper.h"

#include "LArThreeDReco/LArTrackMatching/ThreeDTransverseTracksAlgorithm.h"

using namespace pandora;

namespace lar
{

bool ThreeDTransverseTracksAlgorithm::SortByNMatchedSamplingPoints(const TensorType::Element &lhs, const TensorType::Element &rhs)
{
    if (lhs.GetOverlapResult().GetNMatchedSamplingPoints() != rhs.GetOverlapResult().GetNMatchedSamplingPoints())
        return (lhs.GetOverlapResult().GetNMatchedSamplingPoints() > rhs.GetOverlapResult().GetNMatchedSamplingPoints());

    if (lhs.GetOverlapResult().GetNSamplingPoints() != rhs.GetOverlapResult().GetNSamplingPoints())
        return (lhs.GetOverlapResult().GetNSamplingPoints() < rhs.GetOverlapResult().GetNSamplingPoints());

    return (lhs.GetOverlapResult().GetXOverlap().GetXOverlapSpan() < rhs.GetOverlapResult().GetXOverlap().GetXOverlapSpan());
}

//------------------------------------------------------------------------------------------------------------------------------------------

const TwoDSlidingFitResult &ThreeDTransverseTracksAlgorithm::GetCachedSlidingFitResult(Cluster *const pCluster) const
{
    TwoDSlidingFitResultMap::const_iterator iter = m_slidingFitResultMap.find(pCluster);

    if (m_slidingFitResultMap.end() == iter)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDTransverseTracksAlgorithm::UpdateForNewCluster(Cluster *const pNewCluster)
{
    TwoDSlidingFitResult slidingFitResult;
    LArClusterHelper::LArTwoDSlidingFit(pNewCluster, m_slidingFitWindow, slidingFitResult);

    if (!m_slidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pNewCluster, slidingFitResult)).second)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    ThreeDBaseAlgorithm::UpdateForNewCluster(pNewCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDTransverseTracksAlgorithm::UpdateUponDeletion(Cluster *const pDeletedCluster)
{
    TwoDSlidingFitResultMap::iterator iter = m_slidingFitResultMap.find(pDeletedCluster);

    if (m_slidingFitResultMap.end() != iter)
        m_slidingFitResultMap.erase(iter);

    ThreeDBaseAlgorithm::UpdateUponDeletion(pDeletedCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDTransverseTracksAlgorithm::PreparationStep()
{
    ClusterList allClustersList;
    allClustersList.insert(m_clusterListU.begin(), m_clusterListU.end());
    allClustersList.insert(m_clusterListV.begin(), m_clusterListV.end());
    allClustersList.insert(m_clusterListW.begin(), m_clusterListW.end());

    for (ClusterList::const_iterator iter = allClustersList.begin(), iterEnd = allClustersList.end(); iter != iterEnd; ++iter)
    {
        TwoDSlidingFitResult slidingFitResult;
        LArClusterHelper::LArTwoDSlidingFit(*iter, m_slidingFitWindow, slidingFitResult);

        if (!m_slidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(*iter, slidingFitResult)).second)
            throw StatusCodeException(STATUS_CODE_FAILURE);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDTransverseTracksAlgorithm::CalculateOverlapResult(Cluster *pClusterU, Cluster *pClusterV, Cluster *pClusterW)
{
    try
    {
        TransverseOverlapResult overlapResult;
        this->CalculateOverlapResult(pClusterU, pClusterV, pClusterW, overlapResult);

        if (overlapResult.IsInitialized())
            m_overlapTensor.SetOverlapResult(pClusterU, pClusterV, pClusterW, overlapResult);
    }
    catch (StatusCodeException &statusCodeException)
    {
        if (!(STATUS_CODE_NOT_FOUND == statusCodeException.GetStatusCode() || 
              STATUS_CODE_NOT_INITIALIZED == statusCodeException.GetStatusCode()))
            throw statusCodeException;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDTransverseTracksAlgorithm::CalculateOverlapResult(Cluster *pClusterU, Cluster *pClusterV, Cluster *pClusterW,
    TransverseOverlapResult &overlapResult)
{
    TwoDSlidingFitResultMap::const_iterator iterU = m_slidingFitResultMap.find(pClusterU);
    TwoDSlidingFitResultMap::const_iterator iterV = m_slidingFitResultMap.find(pClusterV);
    TwoDSlidingFitResultMap::const_iterator iterW = m_slidingFitResultMap.find(pClusterW);

    if ((m_slidingFitResultMap.end() == iterU) || (m_slidingFitResultMap.end() == iterV) || (m_slidingFitResultMap.end() == iterW))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const TwoDSlidingFitResult &slidingFitResultU(iterU->second);
    const TwoDSlidingFitResult &slidingFitResultV(iterV->second);
    const TwoDSlidingFitResult &slidingFitResultW(iterW->second);

    FitSegmentTensor fitSegmentTensor;
    this->GetFitSegmentTensor(slidingFitResultU, slidingFitResultV, slidingFitResultW, fitSegmentTensor);
    const TransverseOverlapResult transverseOverlapResult(this->GetBestOverlapResult(fitSegmentTensor));

    if (!transverseOverlapResult.IsInitialized())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    if ((transverseOverlapResult.GetMatchedFraction() < m_minOverallMatchedFraction) || (transverseOverlapResult.GetNMatchedSamplingPoints() < m_minOverallMatchedPoints))
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    const int nLayersSpannedU(slidingFitResultU.GetMaxLayer() - slidingFitResultU.GetMinLayer());
    const int nLayersSpannedV(slidingFitResultV.GetMaxLayer() - slidingFitResultV.GetMinLayer());
    const int nLayersSpannedW(slidingFitResultW.GetMaxLayer() - slidingFitResultW.GetMinLayer());
    const unsigned int meanLayersSpanned(static_cast<unsigned int>((1.f / 3.f) * static_cast<float>(nLayersSpannedU + nLayersSpannedV + nLayersSpannedW)));

    if (0 == meanLayersSpanned)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const float nSamplingPointsPerLayer(static_cast<float>(transverseOverlapResult.GetNSamplingPoints()) / static_cast<float>(meanLayersSpanned));

    if (nSamplingPointsPerLayer < m_minSamplingPointsPerLayer)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    overlapResult = transverseOverlapResult;
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
                    const TransverseOverlapResult segmentOverlap(this->GetSegmentOverlap(fitSegmentU, fitSegmentV, fitSegmentW,
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

void ThreeDTransverseTracksAlgorithm::GetFitSegmentList(const TwoDSlidingFitResult &slidingFitResult, FitSegmentList &fitSegmentList) const
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
            ++nSustainedSteps;

            if (2 * nSustainedSteps > m_slidingFitWindow)
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

TransverseOverlapResult ThreeDTransverseTracksAlgorithm::GetSegmentOverlap(const FitSegment &fitSegmentU, const FitSegment &fitSegmentV, const FitSegment &fitSegmentW,
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
            const float uv2w(LArGeometryHelper::MergeTwoPositions(TPC_VIEW_U, TPC_VIEW_V, u, v));
            const float uw2v(LArGeometryHelper::MergeTwoPositions(TPC_VIEW_U, TPC_VIEW_W, u, w));
            const float vw2u(LArGeometryHelper::MergeTwoPositions(TPC_VIEW_V, TPC_VIEW_W, v, w));

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

    const TransverseOverlapResult::XOverlap xOverlapObject(fitSegmentU.GetMinX(), fitSegmentU.GetMaxX(), fitSegmentV.GetMinX(),
        fitSegmentV.GetMaxX(), fitSegmentW.GetMinX(), fitSegmentW.GetMaxX(), xOverlap);

    return TransverseOverlapResult(nMatchedSamplingPoints, nSamplingPoints, pseudoChi2Sum, xOverlapObject);
}

//------------------------------------------------------------------------------------------------------------------------------------------

TransverseOverlapResult ThreeDTransverseTracksAlgorithm::GetBestOverlapResult(const FitSegmentTensor &fitSegmentTensor) const
{
    if (fitSegmentTensor.empty())
        return TransverseOverlapResult();

    unsigned int maxIndexU(0), maxIndexV(0), maxIndexW(0);
    for (FitSegmentTensor::const_iterator tIter = fitSegmentTensor.begin(), tIterEnd = fitSegmentTensor.end(); tIter != tIterEnd; ++tIter)
    {
        maxIndexU = std::max(tIter->first, maxIndexU);

        for (FitSegmentMatrix::const_iterator mIter = tIter->second.begin(), mIterEnd = tIter->second.end(); mIter != mIterEnd; ++mIter)
        {
            maxIndexV = std::max(mIter->first, maxIndexV);

            for (FitSegmentToOverlapResultMap::const_iterator rIter = mIter->second.begin(), rIterEnd = mIter->second.end(); rIter != rIterEnd; ++rIter)
                maxIndexW = std::max(rIter->first, maxIndexW);
        }
    }

    FitSegmentTensor fitSegmentSumTensor;
    for (unsigned int indexU = 0; indexU <= maxIndexU; ++indexU)
    {
        for (unsigned int indexV = 0; indexV <= maxIndexV; ++indexV)
        {
            for (unsigned int indexW = 0; indexW <= maxIndexW; ++indexW)
            {
                TransverseOverlapResultVector transverseOverlapResultVector;
                this->GetPreviousOverlapResults(indexU, indexV, indexW, fitSegmentSumTensor, transverseOverlapResultVector);

                TransverseOverlapResultVector::const_iterator maxElement = std::max_element(transverseOverlapResultVector.begin(), transverseOverlapResultVector.end());
                TransverseOverlapResult maxTransverseOverlapResult((transverseOverlapResultVector.end() != maxElement) ? *maxElement : TransverseOverlapResult());

                TransverseOverlapResult thisOverlapResult;
                if (fitSegmentTensor.count(indexU) && (fitSegmentTensor.find(indexU))->second.count(indexV) && ((fitSegmentTensor.find(indexU))->second.find(indexV))->second.count(indexW))
                    thisOverlapResult = (((fitSegmentTensor.find(indexU))->second.find(indexV))->second.find(indexW))->second;

                if (!thisOverlapResult.IsInitialized() && !maxTransverseOverlapResult.IsInitialized())
                    continue;

                fitSegmentSumTensor[indexU][indexV][indexW] = thisOverlapResult + maxTransverseOverlapResult;
            }
        }
    }

    TransverseOverlapResult bestTransverseOverlapResult;
    for (unsigned int indexU = 0; indexU <= maxIndexU; ++indexU)
    {
        for (unsigned int indexV = 0; indexV <= maxIndexV; ++indexV)
        {
            for (unsigned int indexW = 0; indexW <= maxIndexW; ++indexW)
            {
                const TransverseOverlapResult &transverseOverlapResult(fitSegmentSumTensor[indexU][indexV][indexW]);

                if (!transverseOverlapResult.IsInitialized())
                    continue;

                if (transverseOverlapResult > bestTransverseOverlapResult)
                    bestTransverseOverlapResult = transverseOverlapResult;
            }
        }
    }

    return bestTransverseOverlapResult;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDTransverseTracksAlgorithm::GetPreviousOverlapResults(const unsigned int indexU, const unsigned int indexV, const unsigned int indexW,
    FitSegmentTensor &fitSegmentSumTensor, TransverseOverlapResultVector &transverseOverlapResultVector) const
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

        const TransverseOverlapResult &transverseOverlapResult(fitSegmentSumTensor[newIndexU][newIndexV][newIndexW]);

        if (transverseOverlapResult.IsInitialized())
            transverseOverlapResultVector.push_back(transverseOverlapResult);
    }

    if (transverseOverlapResultVector.empty())
        transverseOverlapResultVector.push_back(TransverseOverlapResult());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDTransverseTracksAlgorithm::ExamineTensor()
{
    unsigned int repeatCounter(0);

    for (TensorManipulationToolList::const_iterator iter = m_algorithmToolList.begin(), iterEnd = m_algorithmToolList.end(); iter != iterEnd; )
    {
        if ((*iter)->Run(this, m_overlapTensor))
        {
            iter = m_algorithmToolList.begin();

            if (++repeatCounter > m_nMaxTensorToolRepeats)
                break;
        }
        else
        {
            ++iter;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDTransverseTracksAlgorithm::TidyUp()
{
    m_slidingFitResultMap.clear();
    return ThreeDBaseAlgorithm<TransverseOverlapResult>::TidyUp();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDTransverseTracksAlgorithm::FitSegment::FitSegment(const TwoDSlidingFitResult &twoDSlidingFitResult,
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
    AlgorithmToolList algorithmToolList;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle,
        "TrackTools", algorithmToolList));

    for (AlgorithmToolList::const_iterator iter = algorithmToolList.begin(), iterEnd = algorithmToolList.end(); iter != iterEnd; ++iter)
    {
        TensorManipulationTool *pTensorManipulationTool(dynamic_cast<TensorManipulationTool*>(*iter));

        if (NULL == pTensorManipulationTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolList.push_back(pTensorManipulationTool);
    }

    m_nMaxTensorToolRepeats = 5000;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NMaxTensorToolRepeats", m_nMaxTensorToolRepeats));

    m_slidingFitWindow = 20;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));

    m_pseudoChi2Cut = 3.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PseudoChi2Cut", m_pseudoChi2Cut));

    m_minSegmentMatchedFraction = 0.1f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinSegmentMatchedFraction", m_minSegmentMatchedFraction));

    m_minSegmentMatchedPoints = 3;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinSegmentMatchedPoints", m_minSegmentMatchedPoints));

    m_minOverallMatchedFraction = 0.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinOverallMatchedFraction", m_minOverallMatchedFraction));

    m_minOverallMatchedPoints = 10;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinOverallMatchedPoints", m_minOverallMatchedPoints));

    m_minSamplingPointsPerLayer = 0.1f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinSamplingPointsPerLayer", m_minSamplingPointsPerLayer));

    return ThreeDBaseAlgorithm<TransverseOverlapResult>::ReadSettings(xmlHandle);
}

} // namespace lar
