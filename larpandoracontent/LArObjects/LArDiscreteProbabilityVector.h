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
#include <cmath>

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

    typedef std::vector<float> ResamplingPoints;

    /**
     *  @brief  Constructor
     *
     *  @param  param description
     */
    template <typename TX, typename TY>
    DiscreteProbabilityVector(InputData<TX, TY> const &inputData, const TX xUpperBound);

    /**
     *  @brief  Constructor
     *
     *  @param  param description
     */
    DiscreteProbabilityVector(DiscreteProbabilityVector const &discreteProbabilityVector, ResamplingPoints const &resamplingPoints);

    /**
     *  @brief  Evaluate cumulative probability at arbritrary x
     *
     *  @param  x the x value
     */
    float EvaluateCumulativeProbability(float x) const;

    void Print();


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

            float GetX() const;
            float GetDensityDatum() const;
            float GetCumulativeDatum() const;

        private:
            const float m_x;                     ///< The x coordinate
            const float m_densityDatum;          ///< The probability density value
            const float m_cumulativeDatum;       ///< The cumulative probability value
    };
    typedef std::vector<DiscreteProbabilityDatum> DiscreteProbabilityData;

    typedef std::vector<std::pair<float, float > > DiscreteCumulativeProbabilityData;

    template<typename TX, typename TY>
    DiscreteProbabilityData InitialiseDiscreteProbabilityData(DiscreteProbabilityVector::InputData<TX, TY> inputData);

    DiscreteProbabilityData ResampleDiscreteProbabilityData(DiscreteProbabilityVector const &discreteProbabilityVector, ResamplingPoints const &resamplingPoints);

    template <typename TX, typename TY>
    static bool SortInputDataByX(InputDatum<TX, TY> lhs, InputDatum<TX, TY> rhs);

    template <typename TX, typename TY>
    float CalculateNormalisation(InputData<TX, TY> const &inputData);

    void VerifyCompleteData();

    const float m_xUpperBound;
    const DiscreteProbabilityData m_discreteProbabilityData;
};

inline void DiscreteProbabilityVector::Print()
{
    for (size_t iDatum = 0; iDatum < m_discreteProbabilityData.size(); ++iDatum)
    {
        std::cout<<"Datum: " << iDatum << "  x: " << m_discreteProbabilityData.at(iDatum).GetX() << "  density: " << m_discreteProbabilityData.at(iDatum).GetDensityDatum() << "  cumulative: " << m_discreteProbabilityData.at(iDatum).GetCumulativeDatum() << std::endl;
    }
}

inline DiscreteProbabilityVector::DiscreteProbabilityDatum::DiscreteProbabilityDatum(const float &x, const float &densityDatum, const float &cumulativeDatum) :
    m_x(x),
    m_densityDatum(densityDatum),
    m_cumulativeDatum(cumulativeDatum)
{
}

inline float DiscreteProbabilityVector::DiscreteProbabilityDatum::GetX() const
{
    return m_x;
}

inline float DiscreteProbabilityVector::DiscreteProbabilityDatum::GetDensityDatum() const
{
    return m_densityDatum;
}

inline float DiscreteProbabilityVector::DiscreteProbabilityDatum::GetCumulativeDatum() const
{
    return m_cumulativeDatum;
}




inline void DiscreteProbabilityVector::VerifyCompleteData()
{
    if (2 > m_discreteProbabilityData.size())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    if (m_discreteProbabilityData.back().GetX() > m_xUpperBound)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    return;
}


} // namespace lar_content

#endif // £ifndef  LAR_DISCRETE_PROBABILITY_VECTOR_H

