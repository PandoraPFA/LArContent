/**
 *  @file   larpandoracontent/LArObjects/LArDiscreteProbabilityVector.cc
 *
 *  @brief  Implementation of the lar discrete probability vector class
 *
 *  $Log: $
 */

#include "larpandoracontent/LArObjects/LArDiscreteProbabilityVector.h"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace lar_content
{

template <typename TX, typename TY>
DiscreteProbabilityVector::DiscreteProbabilityVector(InputData<TX, TY> const &inputData, const TX xUpperBound, 
    const bool useWidths) :
    m_xUpperBound(static_cast<float>(xUpperBound)),
    m_useWidths(useWidths),
    m_discreteProbabilityData(this->InitialiseDiscreteProbabilityData(inputData))
{
    this->VerifyCompleteData();
}

//------------------------------------------------------------------------------------------------------------------------------------------

DiscreteProbabilityVector::DiscreteProbabilityVector(DiscreteProbabilityVector const &discreteProbabilityVector,
    std::mt19937 &randomNumberGenerator) :
    m_xUpperBound(discreteProbabilityVector.m_xUpperBound),
    m_useWidths(discreteProbabilityVector.m_useWidths),
    m_discreteProbabilityData(this->RandomiseDiscreteProbabilityData(discreteProbabilityVector, randomNumberGenerator))
{
    this->VerifyCompleteData();
}

//------------------------------------------------------------------------------------------------------------------------------------------

DiscreteProbabilityVector::DiscreteProbabilityVector(DiscreteProbabilityVector const &discreteProbabilityVector,
    ResamplingPoints const &resamplingPoints) :
    m_xUpperBound(discreteProbabilityVector.m_xUpperBound),
    m_useWidths(discreteProbabilityVector.m_useWidths),
    m_discreteProbabilityData(this->ResampleDiscreteProbabilityData(discreteProbabilityVector, resamplingPoints))
{
    this->VerifyCompleteData();
}


//------------------------------------------------------------------------------------------------------------------------------------------

float DiscreteProbabilityVector::EvaluateCumulativeProbability(const float x) const
{
    if (x - m_discreteProbabilityData.back().GetX() > std::numeric_limits<float>::epsilon())
        return 1.f;

    if (x - m_discreteProbabilityData.front().GetX() < std::numeric_limits<float>::epsilon())
        return 0.f;

    for (size_t iDatum = 1; iDatum < m_discreteProbabilityData.size(); ++iDatum)
    {
        if (x - m_discreteProbabilityData.at(iDatum).GetX() > std::numeric_limits<float>::epsilon())
            continue;

        const float xLow(m_discreteProbabilityData.at(iDatum-1).GetX());
        const float yLow(m_discreteProbabilityData.at(iDatum-1).GetCumulativeDatum());
        const float xHigh(m_discreteProbabilityData.at(iDatum).GetX());
        const float yHigh(m_discreteProbabilityData.at(iDatum).GetCumulativeDatum());

        if (std::fabs(xHigh-xLow) < std::numeric_limits<float>::epsilon())
            throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

        const float m((yHigh-yLow)/(xHigh-xLow));
        const float c(yLow-m*xLow);

        return m*x+c;
    }

    throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    return 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename TX, typename TY>
DiscreteProbabilityVector::DiscreteProbabilityData DiscreteProbabilityVector::InitialiseDiscreteProbabilityData(
    InputData<TX, TY> inputData) const
{
    if (2 > inputData.size())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    std::sort(inputData.begin(), inputData.end(), DiscreteProbabilityVector::SortInputDataByX<TX,TY>);

    if (inputData.back().first - m_xUpperBound > std::numeric_limits<float>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    const float normalisation(this->CalculateNormalisation(inputData));
    float accumulationDatum(0.f);

    DiscreteProbabilityData data;
    for (size_t iDatum = 0; iDatum < inputData.size()-1; ++iDatum)
    {
        const float x(static_cast<float>(inputData.at(iDatum).first));
        const float deltaX(static_cast<float>(inputData.at(iDatum+1).first)-x);
        const float densityDatum(static_cast<float>(inputData.at(iDatum).second) / normalisation);
        accumulationDatum+=densityDatum*(m_useWidths ? deltaX : 1.f);
        data.emplace_back(DiscreteProbabilityVector::DiscreteProbabilityDatum(x, densityDatum, accumulationDatum, deltaX));
    }
    const float x(static_cast<float>(inputData.back().first));
    const float deltaX(m_xUpperBound-x);
    const float densityDatum(static_cast<float>(inputData.back().second) / normalisation);
    accumulationDatum+=densityDatum*(m_useWidths ? deltaX : 1.f);
    data.emplace_back(DiscreteProbabilityVector::DiscreteProbabilityDatum(x, densityDatum, accumulationDatum, deltaX));

    return data;
}

//------------------------------------------------------------------------------------------------------------------------------------------

DiscreteProbabilityVector::DiscreteProbabilityData DiscreteProbabilityVector::ResampleDiscreteProbabilityData(
    DiscreteProbabilityVector const & discreteProbabilityVector, ResamplingPoints const & resamplingPoints) const
{
    if (2 > resamplingPoints.size())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    DiscreteProbabilityData resampledProbabilityData;

    float prevCumulativeData(0.f);
    for (size_t iSample = 0; iSample < resamplingPoints.size()-1; ++iSample)
    {
        const float xResampled(resamplingPoints.at(iSample));
        const float deltaX(resamplingPoints.at(iSample+1)-xResampled);
        const float cumulativeDatumResampled(discreteProbabilityVector.EvaluateCumulativeProbability(xResampled));
        const float densityDatumResampled((cumulativeDatumResampled-prevCumulativeData)/(m_useWidths ? deltaX : 1.f));
        resampledProbabilityData.emplace_back(DiscreteProbabilityVector::DiscreteProbabilityDatum(xResampled, 
            densityDatumResampled, cumulativeDatumResampled, deltaX));
        prevCumulativeData = cumulativeDatumResampled;
    }

    const float xResampled(resamplingPoints.back());
    const float deltaX(m_xUpperBound-xResampled);
    const float cumulativeDatumResampled(discreteProbabilityVector.EvaluateCumulativeProbability(xResampled));
    const float densityDatumResampled((cumulativeDatumResampled-prevCumulativeData)/(m_useWidths ? deltaX : 1.f));
    resampledProbabilityData.emplace_back(DiscreteProbabilityVector::DiscreteProbabilityDatum(xResampled, densityDatumResampled, 
        cumulativeDatumResampled, deltaX));

    return resampledProbabilityData;
}

//------------------------------------------------------------------------------------------------------------------------------------------

DiscreteProbabilityVector::DiscreteProbabilityData DiscreteProbabilityVector::RandomiseDiscreteProbabilityData(
   DiscreteProbabilityVector const &discreteProbabilityVector, std::mt19937 &randomNumberGenerator) const
{
    DiscreteProbabilityData randomisedProbabilityData;

    std::vector<size_t> randomisedElements(discreteProbabilityVector.GetSize());
    std::iota (std::begin(randomisedElements), std::end(randomisedElements), 0);
    std::shuffle(std::begin(randomisedElements), std::end(randomisedElements), randomNumberGenerator);

    float xPos(discreteProbabilityVector.GetX(0));
    float cumulativeProbability(0.f);
    for (size_t iElement = 0; iElement < discreteProbabilityVector.GetSize(); ++iElement)
    {
        const size_t randomElementIndex(randomisedElements.at(iElement));
        const float deltaX(discreteProbabilityVector.GetWidth(randomElementIndex));
        const float probabilityDensity(discreteProbabilityVector.GetProbabilityDensity(randomElementIndex));
        cumulativeProbability+=probabilityDensity*(m_useWidths ? deltaX : 1.f);
        randomisedProbabilityData.emplace_back(DiscreteProbabilityVector::DiscreteProbabilityDatum(xPos, probabilityDensity, 
            cumulativeProbability, deltaX));
        xPos+=deltaX;
    }

    return randomisedProbabilityData;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename TX, typename TY>
bool DiscreteProbabilityVector::SortInputDataByX(InputDatum<TX, TY> const &lhs, InputDatum<TX, TY> const &rhs)
{
    float deltaX(static_cast<float>(rhs.first) - static_cast<float>(lhs.first));
    if (std::fabs(deltaX) < std::numeric_limits<float>::epsilon())
        return (lhs.second < rhs.second);

    return (lhs.first < rhs.first);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename TX, typename TY>
float DiscreteProbabilityVector::CalculateNormalisation(InputData<TX, TY> const &inputData) const
{
    if (2 > inputData.size())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    float normalisation(0.f);

    for (size_t iDatum = 0; iDatum < inputData.size()-1; ++iDatum)
    {
        const float y(static_cast<float>(inputData.at(iDatum).second));
        normalisation += y*(m_useWidths ? 
            static_cast<float>(inputData.at(iDatum+1).first) - static_cast<float>(inputData.at(iDatum).first) : 1.f);
    }
    const float y(static_cast<float>(inputData.back().second));
    normalisation += y*(m_useWidths ? m_xUpperBound - static_cast<float>(inputData.back().first) : 1.f);

    return normalisation;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template DiscreteProbabilityVector::DiscreteProbabilityVector(InputData<int, float> const&, int const, bool const);
template DiscreteProbabilityVector::DiscreteProbabilityVector(InputData<float, int> const&, float const, bool const);
template DiscreteProbabilityVector::DiscreteProbabilityVector(InputData<float, float> const&, float const, bool const);
template DiscreteProbabilityVector::DiscreteProbabilityVector(InputData<int, int> const&, int const, bool const);

template DiscreteProbabilityVector::DiscreteProbabilityData DiscreteProbabilityVector::InitialiseDiscreteProbabilityData(
    InputData<int, float>) const;
template DiscreteProbabilityVector::DiscreteProbabilityData DiscreteProbabilityVector::InitialiseDiscreteProbabilityData(
    InputData<float, int>) const;
template DiscreteProbabilityVector::DiscreteProbabilityData DiscreteProbabilityVector::InitialiseDiscreteProbabilityData(
    InputData<float, float>) const;
template DiscreteProbabilityVector::DiscreteProbabilityData DiscreteProbabilityVector::InitialiseDiscreteProbabilityData(
    InputData<int, int>) const;

template bool DiscreteProbabilityVector::SortInputDataByX(InputDatum<int, float> const&, InputDatum<int, float> const&);
template bool DiscreteProbabilityVector::SortInputDataByX(InputDatum<float, int> const&, InputDatum<float, int> const&);
template bool DiscreteProbabilityVector::SortInputDataByX(InputDatum<float, float> const&, InputDatum<float, float> const&);
template bool DiscreteProbabilityVector::SortInputDataByX(InputDatum<int, int> const&, InputDatum<int, int> const&);

template float DiscreteProbabilityVector::CalculateNormalisation(InputData<int, float> const&) const;
template float DiscreteProbabilityVector::CalculateNormalisation(InputData<float, int> const&) const;
template float DiscreteProbabilityVector::CalculateNormalisation(InputData<float, float> const&) const;
template float DiscreteProbabilityVector::CalculateNormalisation(InputData<int, int> const&) const;

} // namespace lar_content
