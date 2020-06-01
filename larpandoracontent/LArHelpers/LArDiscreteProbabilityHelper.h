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

#include <algorithm>
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
     *  @brief  Calculate P value for measured correlation coefficient between two datasets via a integrating the student T dist.
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
inline T LArDiscreteProbabilityHelper::MakeRandomisedSample(const T &t, std::mt19937 &randomNumberGenerator)
{
    return LArDiscreteProbabilityHelper::MakeRandomisedSampleImpl(t, randomNumberGenerator);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <>
inline DiscreteProbabilityVector LArDiscreteProbabilityHelper::MakeRandomisedSampleImpl(const DiscreteProbabilityVector &t, 
    std::mt19937 &randomNumberGenerator)
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
    return LArDiscreteProbabilityHelper::GetSizeImpl(t);
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
    return LArDiscreteProbabilityHelper::GetElementImpl(t, index);
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
