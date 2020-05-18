/**
 *  @file   larpandoracontent/LArHelpers/LArFileHelper.cc
 *
 *  @brief  Implementation of the file helper class.
 *
 *  $Log: $
 */

#include "Pandora/StatusCodes.h"
//
#include "larpandoracontent/LArHelpers/LArDiscreteProbabilityHelper.h"
//
#include <limits>
#include <cmath>
//#include <sys/stat.h>


namespace lar_content
{

template <typename T>
float LArDiscreteProbabilityHelper::CalculateCorrelationCoefficient(const T &t1, const T &t2)
{
    if (GetSize(t1) != (GetSize(t2)))
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    if (3 > GetSize(t1))
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

    variance1 /= (static_cast<float>(GetSize(t1))-1.f);
    variance2 /= (static_cast<float>(GetSize(t2))-1.f);
    covariance /= (static_cast<float>(GetSize(t1))-1.f);

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

template float LArDiscreteProbabilityHelper::CalculateCorrelationCoefficient(DiscreteProbabilityVector const&, DiscreteProbabilityVector const&);
template float LArDiscreteProbabilityHelper::CalculateCorrelationCoefficient(std::vector<float> const&, std::vector<float> const&);

template float LArDiscreteProbabilityHelper::CalculateMean(DiscreteProbabilityVector const&);
template float LArDiscreteProbabilityHelper::CalculateMean(std::vector<float> const&);




} // namespace lar_content
