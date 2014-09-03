/**
 *  @file   LArContent/src/LArThreeDReco/LArTransverseTrackMatching/ThreeDTransverseTracksAlgorithm.cc
 *
 *  @brief  Implementation of the three dimensional transverse tracks algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArGeometryHelper.h"
#include "LArHelpers/LArClusterHelper.h"

#include "LArThreeDReco/LArTransverseTrackMatching/ThreeDTransverseTracksAlgorithm.h"

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
        if (!(STATUS_CODE_NOT_FOUND == statusCodeException.GetStatusCode() || STATUS_CODE_NOT_INITIALIZED == statusCodeException.GetStatusCode()))
            throw statusCodeException;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDTransverseTracksAlgorithm::CalculateOverlapResult(Cluster *pClusterU, Cluster *pClusterV, Cluster *pClusterW,
    TransverseOverlapResult &overlapResult)
{
    const TwoDSlidingFitResult &slidingFitResultU(this->GetCachedSlidingFitResult(pClusterU));
    const TwoDSlidingFitResult &slidingFitResultV(this->GetCachedSlidingFitResult(pClusterV));
    const TwoDSlidingFitResult &slidingFitResultW(this->GetCachedSlidingFitResult(pClusterW));

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
    FitSegmentList fitSegmentListU(slidingFitResultU.GetFitSegmentList());
    FitSegmentList fitSegmentListV(slidingFitResultV.GetFitSegmentList());
    FitSegmentList fitSegmentListW(slidingFitResultW.GetFitSegmentList());

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

    if (xOverlap < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    // Sampling in x
    const float nPointsU(std::fabs((xOverlap / xSpanU) * static_cast<float>(fitSegmentU.GetEndLayer() - fitSegmentU.GetStartLayer())));
    const float nPointsV(std::fabs((xOverlap / xSpanV) * static_cast<float>(fitSegmentV.GetEndLayer() - fitSegmentV.GetStartLayer())));
    const float nPointsW(std::fabs((xOverlap / xSpanW) * static_cast<float>(fitSegmentW.GetEndLayer() - fitSegmentW.GetStartLayer())));

    const unsigned int nPoints(1 + static_cast<unsigned int>((nPointsU + nPointsV + nPointsW) / 3.f));

    // Chi2 calculations
    float pseudoChi2Sum(0.f);
    unsigned int nSamplingPoints(0), nMatchedSamplingPoints(0);

    for (unsigned int n = 0; n <= nPoints; ++n)
    {
        const float x(minX + (maxX - minX) * static_cast<float>(n) / static_cast<float>(nPoints));

        try
        {
            CartesianVector fitUVector(0.f, 0.f, 0.f), fitVVector(0.f, 0.f, 0.f), fitWVector(0.f, 0.f, 0.f);
            CartesianVector fitUDirection(0.f, 0.f, 0.f), fitVDirection(0.f, 0.f, 0.f), fitWDirection(0.f, 0.f, 0.f);
            slidingFitResultU.GetTransverseProjection(x, fitSegmentU, fitUVector, fitUDirection);
            slidingFitResultV.GetTransverseProjection(x, fitSegmentV, fitVVector, fitVDirection);
            slidingFitResultW.GetTransverseProjection(x, fitSegmentW, fitWVector, fitWDirection);
          
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

    const XOverlap xOverlapObject(fitSegmentU.GetMinX(), fitSegmentU.GetMaxX(), fitSegmentV.GetMinX(),
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

    for (TensorToolList::const_iterator iter = m_algorithmToolList.begin(), iterEnd = m_algorithmToolList.end(); iter != iterEnd; )
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
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDTransverseTracksAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    AlgorithmToolList algorithmToolList;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle,
        "TrackTools", algorithmToolList));

    for (AlgorithmToolList::const_iterator iter = algorithmToolList.begin(), iterEnd = algorithmToolList.end(); iter != iterEnd; ++iter)
    {
        TransverseTensorTool *pTransverseTensorTool(dynamic_cast<TransverseTensorTool*>(*iter));

        if (NULL == pTransverseTensorTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolList.push_back(pTransverseTensorTool);
    }

    m_nMaxTensorToolRepeats = 5000;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NMaxTensorToolRepeats", m_nMaxTensorToolRepeats));

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

    return ThreeDTracksBaseAlgorithm<TransverseOverlapResult>::ReadSettings(xmlHandle);
}

} // namespace lar
