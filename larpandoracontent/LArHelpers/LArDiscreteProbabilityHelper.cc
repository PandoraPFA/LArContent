/**
 *  @file   larpandoracontent/LArHelpers/LArDiscreteProbabilityHelper.cc
 *
 *  @brief  Implementation of the discrete probability helper helper class.
 *
 *  $Log: $
 */
#include "Pandora/PandoraInternal.h"

#include "larpandoracontent/LArHelpers/LArDiscreteProbabilityHelper.h"

namespace lar_content
{

template <typename T>
float LArDiscreteProbabilityHelper::CalculateCorrelationCoefficientPValueFromPermutationTest(
    const T &t1, const T &t2, std::mt19937 &randomNumberGenerator, const unsigned int nPermutations)
{
    if (1 > nPermutations)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    const float rNominal(LArDiscreteProbabilityHelper::CalculateCorrelationCoefficient(t1, t2));

    unsigned int nExtreme(0);
    for (unsigned int iPermutation = 0; iPermutation < nPermutations; ++iPermutation)
    {
        const float rRandomised(LArDiscreteProbabilityHelper::CalculateCorrelationCoefficient(
            LArDiscreteProbabilityHelper::MakeRandomisedSample(t1, randomNumberGenerator),
            LArDiscreteProbabilityHelper::MakeRandomisedSample(t2, randomNumberGenerator)));

        if ((rRandomised - rNominal) > std::numeric_limits<float>::epsilon())
            nExtreme++;
    }

    return static_cast<float>(nExtreme) / static_cast<float>(nPermutations);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
float LArDiscreteProbabilityHelper::CalculateCorrelationCoefficientPValueFromStudentTDistribution(
    const T &t1, const T &t2, const unsigned int nIntegrationSteps, const float upperLimit)
{
    const float correlation(LArDiscreteProbabilityHelper::CalculateCorrelationCoefficient(t1, t2));
    const float dof(static_cast<float>(LArDiscreteProbabilityHelper::GetSize(t1)) - 2.f);

    if (0 > dof)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    const float tTestStatisticDenominator(1.f - correlation - correlation);

    if (tTestStatisticDenominator < std::numeric_limits<float>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

    const float tTestStatistic(correlation * std::sqrt(dof) / (std::sqrt(tTestStatisticDenominator)));
    const float tDistCoeff(std::tgamma(0.5f * (dof + 1.f)) / std::tgamma(0.5f * dof) / (std::sqrt(dof * M_PI)));

    const float dx((upperLimit - tTestStatistic) / static_cast<float>(nIntegrationSteps));
    float integral(tDistCoeff * std::pow(1.f + tTestStatistic * tTestStatistic / dof, -0.5f * (dof + 1.f)) +
        tDistCoeff * std::pow(1.f + upperLimit * upperLimit / dof, -0.5f * (dof + 1.f)));
    for (unsigned int iStep = 1; iStep < nIntegrationSteps; ++iStep)
    {
        integral += 2.f * tDistCoeff *
            std::pow(1.f + (tTestStatistic + static_cast<float>(iStep) * dx) * (tTestStatistic + static_cast<float>(iStep) * dx) / dof,
                -0.5f * (dof + 1.f));
    }

    return integral * dx / 2.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
float LArDiscreteProbabilityHelper::CalculateCorrelationCoefficient(const T &t1, const T &t2)
{
    const unsigned int size1(LArDiscreteProbabilityHelper::GetSize(t1));
    const unsigned int size2(LArDiscreteProbabilityHelper::GetSize(t2));
    if (size1 != size2)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    if (2 > size1)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    const float mean1(LArDiscreteProbabilityHelper::CalculateMean(t1));
    const float mean2(LArDiscreteProbabilityHelper::CalculateMean(t2));

    float variance1(0.f), variance2(0.f), covariance(0.f);

    for (unsigned int iElement = 0; iElement < size1; ++iElement)
    {
        const float diff1(LArDiscreteProbabilityHelper::GetElement(t1, iElement) - mean1);
        const float diff2(LArDiscreteProbabilityHelper::GetElement(t2, iElement) - mean2);

        variance1 += diff1 * diff1;
        variance2 += diff2 * diff2;
        covariance += diff1 * diff2;
    }

    if (variance1 < std::numeric_limits<float>::epsilon() || variance2 < std::numeric_limits<float>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

    const float sqrtVars(std::sqrt(variance1 * variance2));
    if (sqrtVars < std::numeric_limits<float>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

    return covariance / sqrtVars;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
float LArDiscreteProbabilityHelper::CalculateMean(const T &t)
{
    const unsigned int size(LArDiscreteProbabilityHelper::GetSize(t));
    if (1 > size)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

    float mean(0.f);
    for (unsigned int iElement = 0; iElement < size; ++iElement)
        mean += LArDiscreteProbabilityHelper::GetElement(t, iElement);

    return mean / static_cast<float>(size);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template float LArDiscreteProbabilityHelper::CalculateCorrelationCoefficientPValueFromPermutationTest(
    const DiscreteProbabilityVector &, const DiscreteProbabilityVector &, std::mt19937 &, const unsigned int);
template float LArDiscreteProbabilityHelper::CalculateCorrelationCoefficientPValueFromPermutationTest(
    const pandora::FloatVector &, const pandora::FloatVector &, std::mt19937 &, const unsigned int);

template float LArDiscreteProbabilityHelper::CalculateCorrelationCoefficientPValueFromStudentTDistribution(
    const DiscreteProbabilityVector &, const DiscreteProbabilityVector &, const unsigned int, const float);
template float LArDiscreteProbabilityHelper::CalculateCorrelationCoefficientPValueFromStudentTDistribution(
    const pandora::FloatVector &, const pandora::FloatVector &, const unsigned int, const float);

template float LArDiscreteProbabilityHelper::CalculateCorrelationCoefficient(const DiscreteProbabilityVector &, const DiscreteProbabilityVector &);
template float LArDiscreteProbabilityHelper::CalculateCorrelationCoefficient(const pandora::FloatVector &, const pandora::FloatVector &);

template float LArDiscreteProbabilityHelper::CalculateMean(const DiscreteProbabilityVector &);
template float LArDiscreteProbabilityHelper::CalculateMean(const pandora::FloatVector &);

} // namespace lar_content
