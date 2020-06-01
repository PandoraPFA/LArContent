/**
 *  @file   larpandoracontent/LArHelpers/LArDiscreteProbabilityHelper.h
 *
 *  @brief  Header file for the discrete probability helper class.
 *
 *  $Log: $
 */
#ifndef LAR_DISCRETE_PROBABILITY_HELPER_HELPER_H
#define LAR_DISCRETE_PROBABILITY_HELPER_HELPER_H 1

#include "larpandoracontent/LArObjects/LArDiscreteProbabilityVector.h"

#include <random>

namespace lar_content
{

/**
 *  @brief  LArDiscreteProbabilityHelper class
 */
class LArDiscreteProbabilityHelper
{
public:
    /**
     *  @brief  Calculate P value for measured correlation coefficient between two datasets via a permutation test 
     *
     *  @param  t1 the first input dataset
     *  @param  t2 the second input dataset
     *  @param  randomNumberGenerator the random number generator to shuffle the datasets
     *  @param  nPermutations the number of permutations to run
     *
     *  @return the p-value
     */
    template <typename T>
    static float CalculateCorrelationCoefficientPValueFromPermutationTest(const T &t1, const T &t2, 
        std::mt19937 &randomNumberGenerator, const size_t nPermutations);

    /**
     *  @brief  Calculate P value for measured correlation coefficient between two datasets via a integrating the student T distribution
     *
     *  @param  t1 the first input dataset
     *  @param  t2 the second input dataset
     *  @param  nIntegrationSteps how many steps to use in the trapezium integration
     *
     *  @return the p-value
     */
    template <typename T>
    static float CalculateCorrelationCoefficientPValueFromStudentTDistribution(const T &t1, const T &t2, 
        const size_t nIntegrationSteps);

    /**
     *  @brief  Calculate the correlation coefficient between two datasets 
     *
     *  @param  t1 the first input dataset
     *  @param  t2 the second input dataset
     *
     *  @return the correlation coefficient
     */
    template <typename T>
    static float CalculateCorrelationCoefficient(const T &t1, const T &t2);

    /**
     *  @brief  Calculate the mean of a dataset 
     *
     *  @param  t the dataset
     *
     *  @return the mean
     */
    template <typename T>
    static float CalculateMean(const T &t);

private:
    /**
     *  @brief  Creates a randomised copy of a dataset 
     *
     *  @param  t the dataset to be shuffled
     *  @param  randomNumberGenerator the random number generator
     *
     *  @return the reshuffled dataset
     */
    template <typename T>
    static T MakeRandomisedSample(const T &t, std::mt19937 &randomNumberGenerator);

    /**
     *  @brief  An implementation for making a randomised copy of a dataset
     *
     *  @param  t the dataset to be shuffled
     *  @param  randomNumberGenerator the random number generator
     *
     *  @return the reshuffled dataset
     */
    template <typename T>
    static T MakeRandomisedSampleImpl(const T &t, std::mt19937 &randomNumberGenerator);

    /**
     *  @brief  An implementation for making a randomised copy of dataset (dataset is an std::vector)
     *
     *  @param  t the std::vector-based dataset to be shuffled
     *  @param  randomNumberGenerator the random number generator
     *
     *  @return the reshuffled std::vector
     */
    template <typename T>
    static std::vector<T> MakeRandomisedSampleImpl(const std::vector<T> &t, std::mt19937 &randomNumberGenerator);

    /**
     *  @brief  Gets the size of a dataset
     *
     *  @param  t the dataset
     *
     *  @return the dataset size
     */
    template <typename T>
    static size_t GetSize(const T &t);

    /**
     *  @brief  An implementation of a getter for the size of a dataset
     *
     *  @param  t the dataset
     *
     *  @return the dataset size
     */
    template <typename T>
    static size_t GetSizeImpl(const T &t);

    /**
     *  @brief  An implementation of a getter for the size of a dataset (dataset is an std::vector)
     *
     *  @param  t the std::vector dataset
     *
     *  @return the std::vector-based dataset size
     */
    template <typename T>
    static size_t GetSizeImpl(const std::vector<T> &t);

    /**
     *  @brief  Get an element in a dataset 
     *
     *  @param  t the dataset
     *  @param  index the index of the element
     *
     *  @return the dataset element
     */
    template <typename T>
    static float GetElement(const T &t, const size_t index);

    /**
     *  @brief  An implementation of a getter for an element in a dataset 
     *
     *  @param  t the dataset
     *  @param  index the index of the element
     *
     *  @return the dataset element
     */
    template <typename T>
    static float GetElementImpl(const T &t, const size_t index);

    /**
     *  @brief  An implementation of a getter for an element in a dataset (dataset is an std::vector)
     *
     *  @param  t the std::vector dataset
     *  @param  index the index of the element
     *
     *  @return the std::vector-based dataset element
     */
    template <typename T>
    static float GetElementImpl(const std::vector<T> &t, const size_t index);
};

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
float LArDiscreteProbabilityHelper::CalculateCorrelationCoefficientPValueFromPermutationTest(const T &t1, const T &t2, 
    std::mt19937 &randomNumberGenerator, const size_t nPermutations)
{
    if (1 > nPermutations)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    float rNominal(CalculateCorrelationCoefficient(t1,t2));

    int nExtreme(0);
    for (size_t iPermutation = 0; iPermutation < nPermutations; ++iPermutation)
    {
        float rRandomised(CalculateCorrelationCoefficient(MakeRandomisedSample(t1,randomNumberGenerator),MakeRandomisedSample(t2,randomNumberGenerator)));
        if ((rRandomised-rNominal) > std::numeric_limits<float>::epsilon())
            nExtreme++;
    }

    return static_cast<float>(nExtreme)/static_cast<float>(nPermutations);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
float LArDiscreteProbabilityHelper::CalculateCorrelationCoefficientPValueFromStudentTDistribution(const T &t1, 
    const T &t2, const size_t nIntegrationSteps)
{
    float correlation(LArDiscreteProbabilityHelper::CalculateCorrelationCoefficient(t1,t2));
    float dof(static_cast<float>(GetSize(t1)) - 2.f);
    float tTestStatistic(correlation*sqrt(dof)/(sqrt(1.-correlation*correlation)));
    float tDistCoeff(std::tgamma(0.5 * (dof+1.)) / std::tgamma(0.5*dof) / (std::sqrt(dof*M_PI)));

    int nSteps(nIntegrationSteps);
    float upperLimit(15.f);
    float dx((upperLimit-tTestStatistic)/static_cast<float>(nSteps));
    float integral(tDistCoeff*std::pow( 1.0 + tTestStatistic*tTestStatistic/dof, -0.5 *(dof + 1.0)) + 
            tDistCoeff*std::pow( 1.0 + upperLimit*upperLimit/dof, -0.5 *(dof + 1.0)));
    for (int iStep = 1; iStep < nSteps; iStep++)
        integral+=2. * tDistCoeff*std::pow( 
            1.0 + (tTestStatistic + static_cast<float>(iStep)*dx)*(tTestStatistic + static_cast<float>(iStep)*dx)/dof, -0.5 *(dof + 1.0));
    integral *= dx/2.0;

    return integral;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
float LArDiscreteProbabilityHelper::CalculateCorrelationCoefficient(const T &t1, const T &t2)
{
    if (GetSize(t1) != (GetSize(t2)))
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    if (2 > GetSize(t1))
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    float mean1(CalculateMean(t1)); 
    float mean2(CalculateMean(t2));

    float variance1(0.f);
    float variance2(0.f);
    float covariance(0.f);

    for (size_t iElement = 0; iElement < GetSize(t1); iElement++)
    {
        float element1(GetElement(t1,iElement));
        float element2(GetElement(t2,iElement));

        variance1 += (element1-mean1)*(element1-mean1);
        variance2 += (element2-mean2)*(element2-mean2);
        covariance += (element1-mean1)*(element2-mean2);
    }

    float sqrtVars(std::sqrt(variance1*variance2));
    if(sqrtVars < std::numeric_limits<float>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

    return covariance /= sqrtVars;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
float LArDiscreteProbabilityHelper::CalculateMean(const T &t)
{
    if (0 == GetSize(t))
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

    float mean(0.f);
    for (size_t iElement = 0; iElement < GetSize(t); ++iElement)
    {
        mean+=GetElement(t,iElement);
    }
    mean /= static_cast<float>(GetSize(t));

    return mean;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline T LArDiscreteProbabilityHelper::MakeRandomisedSample(const T &t, std::mt19937 &randomNumberGenerator)
{
    return MakeRandomisedSampleImpl(t, randomNumberGenerator);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <>
inline DiscreteProbabilityVector LArDiscreteProbabilityHelper::MakeRandomisedSampleImpl(const DiscreteProbabilityVector &t, std::mt19937 &randomNumberGenerator)
{
    return DiscreteProbabilityVector(t,randomNumberGenerator);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline std::vector<T> LArDiscreteProbabilityHelper::MakeRandomisedSampleImpl(const std::vector<T> &t, std::mt19937 &randomNumberGenerator)
{
    std::vector<T> randomisedVector(t);
    std::shuffle(randomisedVector.begin(), randomisedVector.end(), randomNumberGenerator);

    return randomisedVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline size_t LArDiscreteProbabilityHelper::GetSize(const T &t)
{
    return GetSizeImpl(t);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <>
inline size_t LArDiscreteProbabilityHelper::GetSizeImpl(const DiscreteProbabilityVector& t)
{
    return t.GetSize();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline size_t LArDiscreteProbabilityHelper::GetSizeImpl(const std::vector<T> &t)
{
    return t.size();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline float LArDiscreteProbabilityHelper::GetElement(const T &t, const size_t index)
{
    return GetElementImpl(t, index);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <>
inline float LArDiscreteProbabilityHelper::GetElementImpl(const DiscreteProbabilityVector& t, const size_t index)
{
    return static_cast<float>(t.GetProbability(index));
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
inline float LArDiscreteProbabilityHelper::GetElementImpl(const std::vector<T> &t, const size_t index)
{
    return static_cast<float>(t.at(index));
}

} // namespace lar_content
#endif // #ifndef LAR_DISCRETE_PROBABILITY_HELPER_H
