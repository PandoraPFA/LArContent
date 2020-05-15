/**
 *  @file   larpandoracontent/LArObjects/LArDiscreteProbabilityVector.cc
 *
 *  @brief  Implementation of the lar discrete probability vetor class
 *
 *  $Log: $
 */

//#include "Pandora/PandoraInputTypes.h"

#include "larpandoracontent/LArObjects/LArDiscreteProbabilityVector.h"

namespace lar_content
{

template <typename TX, typename TY>
DiscreteProbabilityVector::DiscreteProbabilityVector(InputData<TX, TY> const &inputData, const TX xUpperBound) :
    m_xUpperBound(static_cast<float>(xUpperBound)),
    m_discreteProbabilityData(InitialiseDiscreteProbabilityData(inputData))
{
    VerifyCompleteData();
}

//------------------------------------------------------------------------------------------------------------------------------------------

DiscreteProbabilityVector::DiscreteProbabilityVector(DiscreteProbabilityVector const &discreteProbabilityVector, ResamplingPoints const &resamplingPoints) :
    m_xUpperBound(discreteProbabilityVector.m_xUpperBound),
    m_discreteProbabilityData(ResampleDiscreteProbabilityData(discreteProbabilityVector, resamplingPoints))
{
    VerifyCompleteData();
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DiscreteProbabilityVector::EvaluateCumulativeProbability(const float x) const
{
    if (x > m_discreteProbabilityData.back().GetX())
        return 1.f;

    if (x < m_discreteProbabilityData.front().GetX())
        return 0.f;

    for (size_t iDatum = 1; iDatum < m_discreteProbabilityData.size(); ++iDatum)
    {
        if (x > m_discreteProbabilityData.at(iDatum).GetX())
            continue;

        float xLow(m_discreteProbabilityData.at(iDatum-1).GetX());
        float yLow(m_discreteProbabilityData.at(iDatum-1).GetCumulativeDatum());
        float xHigh(m_discreteProbabilityData.at(iDatum).GetX());
        float yHigh(m_discreteProbabilityData.at(iDatum).GetCumulativeDatum());

        if (std::fabs(xHigh-xLow) < std::numeric_limits<float>::epsilon())
            throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

        float m((yHigh-yLow)/(xHigh-xLow));
        float c(yLow-m*xLow);

        return m*x+c;
    }

    throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    return 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename TX, typename TY>
DiscreteProbabilityVector::DiscreteProbabilityData DiscreteProbabilityVector::InitialiseDiscreteProbabilityData(DiscreteProbabilityVector::InputData<TX, TY> inputData)
{
    if (2 > inputData.size())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);


    std::sort(inputData.begin(), inputData.end(), DiscreteProbabilityVector::SortInputDataByX<TX,TY>);

    if (inputData.back().first > m_xUpperBound)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    float normalisation(CalculateNormalisation(inputData));
    float accumulationDatum(0.f);

    DiscreteProbabilityData data;
    for (size_t iDatum = 0; iDatum < inputData.size(); ++iDatum)
    {
        float x(inputData.at(iDatum).first);
        float densityDatum(static_cast<float>(inputData.at(iDatum).second) / normalisation);
        accumulationDatum+=densityDatum;
        data.emplace_back(DiscreteProbabilityVector::DiscreteProbabilityDatum(x, densityDatum, accumulationDatum));
    }
    return data;
}

//------------------------------------------------------------------------------------------------------------------------------------------


DiscreteProbabilityVector::DiscreteProbabilityData DiscreteProbabilityVector::ResampleDiscreteProbabilityData(DiscreteProbabilityVector const & discreteProbabilityVector, ResamplingPoints const & resamplingPoints)
{
    DiscreteProbabilityData resampledProbabilityData;

    float prevCumulativeData(0.f);
    for (size_t iSample = 0; iSample < resamplingPoints.size(); iSample++)
    {
        float xResampled(resamplingPoints.at(iSample));
        float cumulativeDatumResampled(discreteProbabilityVector.EvaluateCumulativeProbability(xResampled));
        float densityDatumResampled(cumulativeDatumResampled-prevCumulativeData);
        resampledProbabilityData.emplace_back(DiscreteProbabilityVector::DiscreteProbabilityDatum(xResampled, densityDatumResampled, cumulativeDatumResampled));
        prevCumulativeData = cumulativeDatumResampled;
    }

    return resampledProbabilityData;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename TX, typename TY>
bool DiscreteProbabilityVector::SortInputDataByX(InputDatum<TX, TY> lhs, InputDatum<TX, TY> rhs)
{
    float deltaX(static_cast<float>(rhs.first) - static_cast<float>(lhs.first));
    if (std::fabs(deltaX) < std::numeric_limits<float>::epsilon())
        return (lhs.second < rhs.second);

    return (lhs.first < rhs.first);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename TX, typename TY>
float DiscreteProbabilityVector::CalculateNormalisation(InputData<TX, TY> const &inputData)
{
    if (2 > inputData.size())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    float normalisation(0.f);

    for (size_t iDatum = 0; iDatum < inputData.size()-1; ++iDatum)
    {
        float deltaX(static_cast<float>(inputData.at(iDatum+1).first) - static_cast<float>(inputData.at(iDatum).first));
        float y(static_cast<float>(inputData.at(iDatum).second));
        normalisation += y*deltaX;
    }
    float deltaX(m_xUpperBound - static_cast<float>(inputData.back().first));
    float y(static_cast<float>(inputData.back().second));
    normalisation += y*deltaX;

    return normalisation;
}

template DiscreteProbabilityVector::DiscreteProbabilityVector(DiscreteProbabilityVector::InputData<int, float> const&, int const);
template DiscreteProbabilityVector::DiscreteProbabilityVector(DiscreteProbabilityVector::InputData<float, int> const&, float const);
template DiscreteProbabilityVector::DiscreteProbabilityVector(DiscreteProbabilityVector::InputData<float, float> const&, float const);
template DiscreteProbabilityVector::DiscreteProbabilityVector(DiscreteProbabilityVector::InputData<int, int> const&, int const);

template DiscreteProbabilityVector::DiscreteProbabilityData DiscreteProbabilityVector::InitialiseDiscreteProbabilityData(DiscreteProbabilityVector::InputData<int, float>);
template DiscreteProbabilityVector::DiscreteProbabilityData DiscreteProbabilityVector::InitialiseDiscreteProbabilityData(DiscreteProbabilityVector::InputData<float, int>);
template DiscreteProbabilityVector::DiscreteProbabilityData DiscreteProbabilityVector::InitialiseDiscreteProbabilityData(DiscreteProbabilityVector::InputData<float, float>);
template DiscreteProbabilityVector::DiscreteProbabilityData DiscreteProbabilityVector::InitialiseDiscreteProbabilityData(DiscreteProbabilityVector::InputData<int, int>);

template bool DiscreteProbabilityVector::SortInputDataByX(DiscreteProbabilityVector::InputDatum<int, float>, DiscreteProbabilityVector::InputDatum<int, float>);
template bool DiscreteProbabilityVector::SortInputDataByX(DiscreteProbabilityVector::InputDatum<float, int>, DiscreteProbabilityVector::InputDatum<float, int>);
template bool DiscreteProbabilityVector::SortInputDataByX(DiscreteProbabilityVector::InputDatum<float, float>, DiscreteProbabilityVector::InputDatum<float, float>);
template bool DiscreteProbabilityVector::SortInputDataByX(DiscreteProbabilityVector::InputDatum<int, int>, DiscreteProbabilityVector::InputDatum<int, int>);

template float DiscreteProbabilityVector::CalculateNormalisation(InputData<int, float> const&);
template float DiscreteProbabilityVector::CalculateNormalisation(InputData<float, int> const&);
template float DiscreteProbabilityVector::CalculateNormalisation(InputData<float, float> const&);
template float DiscreteProbabilityVector::CalculateNormalisation(InputData<int, int> const&);






} // namespace lar_content
