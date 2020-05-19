/**
 *  @file   larpandoracontent/LArHelpers/LArDiscreteProbabilityHelper.h
 *
 *  @brief  Header file for the discrete probability helper helper class.
 *
 *  $Log: $
 */
#ifndef LAR_DISCRETE_PROBABILITY_HELPER_HELPER_H
#define LAR_DISCRETE_PROBABILITY_HELPER_HELPER_H 1

#include "larpandoracontent/LArObjects/LArDiscreteProbabilityVector.h"

/**
 *  @brief  DiscreteCumulativeDistributionHelper class
 */

namespace lar_content
{

class LArDiscreteProbabilityHelper
{
public:

    template <typename T>
    static float CalculateCorrelationCoefficient(const T &t1, const T &t2);

    template <typename T>
    static float CalculateMean(const T &t);

private:



    template <typename T>
    static size_t GetSize(const T &t);

    template <typename T>
    static size_t GetSizeImpl(const T &t);

    template <typename T>
    static size_t GetSizeImpl(const std::vector<T> &t);

    template <typename T>
    static float GetElement(const T &t, const size_t element);

    template <typename T>
    static float GetElementImpl(const T &t, const size_t element);

    template <typename T>
    static float GetElementImpl(const std::vector<T> &t, const size_t element);
   /**
    *  @brief  CumulDistLinearInterpolation
    *
    *  @param  x position for which we want the interpolated y
    *  @param  the discrete cumulatie distribution one wishes to interpolate
    *
    *  @return 
    */
   //static float CumulDistLinearInterpolation(const float &xPos, const DiscreteCumulativeDistribution &distribution);
};


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

    if(variance1*variance2 < std::numeric_limits<float>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

    return covariance /= std::sqrt(variance1*variance2);
}

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

template <typename T>
inline size_t LArDiscreteProbabilityHelper::GetSize(const T &t)
{
    return GetSizeImpl(t);
}

template <>
inline size_t LArDiscreteProbabilityHelper::GetSizeImpl(const DiscreteProbabilityVector& t)
{
    return t.GetSize();
}

template <typename T>
inline size_t LArDiscreteProbabilityHelper::GetSizeImpl(const std::vector<T> &t)
{
    return t.size();
}

template <typename T>
inline float LArDiscreteProbabilityHelper::GetElement(const T &t, const size_t element)
{
    return GetElementImpl(t, element);
}

template <>
inline float LArDiscreteProbabilityHelper::GetElementImpl(const DiscreteProbabilityVector& t, const size_t element)
{
    return t.GetProbabilityDensity(element);
}

template<typename T>
inline float LArDiscreteProbabilityHelper::GetElementImpl(const std::vector<T> &t, const size_t element)
{
    return t.at(element);
}

} // namespace lar_content
#endif // #ifndef LAR_DISCRETE_PROBABILITY_HELPER_HELPER_H
