/**
 *  @file   larpandoracontent/LArObjects/LArDiscreteProbabilityVector.h
 *
 *  @brief  Header file for the lar discrete cumulative distribution class.
 *
 *  $Log: $
 */
#ifndef LAR_DISCRETE_PROBABILITY_VECTOR_H
#define LAR_DISCRETE_PROBABILITY_VECTOR_H 1

#include <vector>
#include <iostream>

#include "Pandora/PandoraInternal.h"
#include "Pandora/StatusCodes.h"

namespace lar_content
{

/**
 *  @brief  DiscreteProbabilityVector class
 */
class DiscreteProbabilityVector
{
public:

    template <typename TX, typename TY>
    using InputDatum = std::pair<TX, TY>;

    template <typename TX, typename TY>
    using InputData = std::vector<InputDatum<TX, TY> >;

    /**
     *  @brief  Constructor
     *
     *  @param  param description
     */
    template <typename TX, typename TY>
    DiscreteProbabilityVector(InputData<TX, TY> inputData);

    /**
     *  @brief  Store data to later convert to a cumulative distribution
     *
     *  @param  x independent variable
     *  @param  y dependent variable
     */
    //template <typename T1, typename T2>
    //void CollectInputData(T1 x, T2 y);


private:

    class DiscreteProbabilityDatum
    {
        public:
            DiscreteProbabilityDatum(const float &x, const float &densityDatum, const float &cumulativeDatum);

            float GetX();
            float GetDensityDatum();
            float GetCumulativeDatum();

        private:
            float m_x;                     ///< The x coordinate
            float m_densityDatum;          ///< The probability density value
            float m_cumulativeDatum;       ///< The cumulative probability value
    };
    typedef std::vector<DiscreteProbabilityDatum> DiscreteProbabilityData;

    template<typename TX, typename TY>
    DiscreteProbabilityData InitialiseDiscreteProbabilityData(DiscreteProbabilityVector::InputData<TX, TY> &inputData);

    template <typename TX, typename TY>
    static bool SortInputDataByX(InputDatum<TX, TY> lhs, InputDatum<TX, TY> rhs);

    template <typename TX, typename TY>
    float CalculateNormalisation(InputData<TX, TY> const &inputData);

    const DiscreteProbabilityData m_discreteProbabilityData;
};

template <typename TX, typename TY>
inline DiscreteProbabilityVector::DiscreteProbabilityVector(DiscreteProbabilityVector::InputData<TX, TY> inputData) :
    m_discreteProbabilityData(InitialiseDiscreteProbabilityData(inputData))
{
}

inline DiscreteProbabilityVector::DiscreteProbabilityDatum::DiscreteProbabilityDatum(const float &x, const float &densityDatum, const float &cumulativeDatum) :
    m_x(x),
    m_densityDatum(densityDatum),
    m_cumulativeDatum(cumulativeDatum)
{
}

inline float DiscreteProbabilityVector::DiscreteProbabilityDatum::GetX()
{
    return m_x;
}

inline float DiscreteProbabilityVector::DiscreteProbabilityDatum::GetDensityDatum()
{
    return m_densityDatum;
}

inline float DiscreteProbabilityVector::DiscreteProbabilityDatum::GetCumulativeDatum()
{
    return m_cumulativeDatum;
}

template <typename TX, typename TY>
inline DiscreteProbabilityVector::DiscreteProbabilityData DiscreteProbabilityVector::InitialiseDiscreteProbabilityData(DiscreteProbabilityVector::InputData<TX, TY> &inputData)
{
    std::sort(inputData.begin(), inputData.end(), DiscreteProbabilityVector::SortInputDataByX<TX,TY>);

    float normalisation(CalculateNormalisation(inputData));
    float accumulationDatum(0.f);

    DiscreteProbabilityData data;
    for (size_t iDatum = 0; iDatum < inputData.size()-1; ++iDatum)
    {
        float x(inputData.at(iDatum).first);
        float densityDatum(static_cast<float>(inputData.at(iDatum).second) / normalisation);
        accumulationDatum+=densityDatum;
        data.emplace_back(DiscreteProbabilityVector::DiscreteProbabilityDatum(x, densityDatum, accumulationDatum));
        std::cout<<iDatum << "  " << x << "  " << densityDatum << "  " << accumulationDatum << std::endl;
    }
    return data;
}

template <typename TX, typename TY>
inline bool DiscreteProbabilityVector::SortInputDataByX(InputDatum<TX, TY> lhs, InputDatum<TX, TY> rhs)
{
    float deltaX(static_cast<float>(rhs.first) - static_cast<float>(lhs.first));
    if (std::fabs(deltaX) < std::numeric_limits<float>::epsilon())
        return (lhs.second < rhs.second);

    return (lhs.first < rhs.first);
}

template <typename TX, typename TY>
inline float DiscreteProbabilityVector::CalculateNormalisation(InputData<TX, TY> const &inputData)
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

    return normalisation;
}


} // namespace lar_content

#endif // £ifndef  LAR_DISCRETE_PROBABILITY_VECTOR_H

