/**
 *  @file   larpandoracontent/LArObjects/LArDiscreteProbabilityVector.h
 *
 *  @brief  Header file for the lar discrete probability vector class
 *
 *  $Log: $
 */
#ifndef LAR_DISCRETE_PROBABILITY_VECTOR_H
#define LAR_DISCRETE_PROBABILITY_VECTOR_H 1

#include "Pandora/StatusCodes.h"

#include <utility>
#include <vector>


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
     *  @param  inputData the data used to construt the probability vector
     *  @param  the upper bound of the probability vector
     */
    template <typename TX, typename TY>
    DiscreteProbabilityVector(InputData<TX, TY> const &inputData, const TX xUpperBound);

    /**
     *  @brief  Constructor
     *
     *  @param  discreteProbabilityVector a discrete probability vector to resample from
     *  @param  The points to resample the discrete probability vector with 
     */
    DiscreteProbabilityVector(DiscreteProbabilityVector const &discreteProbabilityVector, ResamplingPoints const &resamplingPoints);

    /**
     *  @brief  Evaluate the cumulative probability at arbitrary x
     *
     *  @param  x the x value
     *
     *  @return the cumulative probability
     */
    float EvaluateCumulativeProbability(float x) const;

    /**
     *  @brief  Get the size of the probability vector
     *
     *  @return the probability vector size
     */
    size_t GetSize() const;

    /**
     *  @brief  Get the x value of the element in the vector
     *
     *  @param  index the index in the vector
     *
     *  @return the x value
     */
    float GetX(const size_t index) const;

    /**
     *  @brief  Get the probability density value of the element in the vector
     *
     *  @param  index the index in the vector
     *
     *  @return the probablity density
     */
    float GetProbabilityDensity(const size_t index) const;

    /**
     *  @brief  Get the cumulative probability value of the element in the vector
     *
     *  @param  index the index in the vector
     *
     *  @return the cumulative probability
     */
    float GetCumulativeProbability(const size_t index) const;

    /**
     *  @brief  Get all information stored at a particular index
     *
     *  @param  index the index in the vector
     *  @param  x the x value
     *  @param  probabilityDensity the probability density value
     *  @param  cumulativeProbability the cumulative probability value
     */
    void GetAllAtIndex(const size_t index, float &x, float &probabilityDensity, float &cumulativeProbability) const;

private:

    /**
     *  @brief  DiscreteProbabilityData class
     */
    class DiscreteProbabilityDatum
    {
        public:
            /**
             *  @brief  Constructor
             *
             *  @param  x the x value
             *  @param  densityDatum the probability density for the corresponding x
             *  @param  cumulativeDatum the cumulative probability for the corresponding x
             */
            DiscreteProbabilityDatum(const float &x, const float &densityDatum, const float &cumulativeDatum);

            /**
             *  @brief  Get the x value for the datum
             *
             *  @return the x value
             */
            float GetX() const;

            /**
             *  @brief  Get the probability density for the datum
             *
             *  @return the probability density
             */
            float GetDensityDatum() const;

            /**
             *  @brief  Get the cumulative probability for the datum
             *
             *  @return the cumulative probability
             */
            float GetCumulativeDatum() const;

        private:
            const float m_x;                     ///< The x coordinate
            const float m_densityDatum;          ///< The probability density value
            const float m_cumulativeDatum;       ///< The cumulative probability value
    };

    typedef std::vector<DiscreteProbabilityDatum> DiscreteProbabilityData;

    typedef std::vector<std::pair<float, float > > DiscreteCumulativeProbabilityData;

    /**
     *  @brief  Get a initialised probability data vector from the input data
     *
     *  @param  inputData the input data
     *
     *  @return a fully-initialised discrete probability data vector
     */
    template<typename TX, typename TY>
    DiscreteProbabilityData InitialiseDiscreteProbabilityData(InputData<TX, TY> inputData) const;

    /**
     *  @brief  Get a resampled probability data vector by resampling another probability data vector
     *
     *  @param  discreteProbabilityVector another discrete probability vector
     *  @param  resamplingPoints the points to resample the discrete probability vector with
     *
     *  @return a resampled probability data vector
     */
    DiscreteProbabilityData ResampleDiscreteProbabilityData(DiscreteProbabilityVector const &discreteProbabilityVector, 
        ResamplingPoints const &resamplingPoints) const;

    /**
     *  @brief  Sort the input data according to their x value
     *
     *  @param  lhs the first InputDatum
     *  @param  lhs the second InputDatum
     *
     *  @return a bool dictacting swapping the two InputDatums
     */
    template <typename TX, typename TY>
    static bool SortInputDataByX(InputDatum<TX, TY> lhs, InputDatum<TX, TY> rhs);

    /**
     *  @brief  Calculate the probability normalisation
     *
     *  @param  inputData the input data
     *
     *  @return the probability normalisation
     */
    template <typename TX, typename TY>
    float CalculateNormalisation(InputData<TX, TY> const &inputData) const;

    /**
     *  @brief  Verify the integrity of the complete probability vector
     */
    void VerifyCompleteData() const;

    /**
     *  @brief  Verify the integrity of the element request
     *
     *  @param  index the index in the probability vector
     */
    void VerifyElementRequest(const size_t index) const;

    const float m_xUpperBound;                                              ///< the upper bound of the probability vector
    const DiscreteProbabilityData m_discreteProbabilityData;                ///< the probability data
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline size_t DiscreteProbabilityVector::GetSize() const
{
    return m_discreteProbabilityData.size();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DiscreteProbabilityVector::GetX(const size_t index) const
{
    VerifyElementRequest(index);

    return m_discreteProbabilityData.at(index).GetX();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DiscreteProbabilityVector::GetProbabilityDensity(const size_t index) const
{
    VerifyElementRequest(index);

    return m_discreteProbabilityData.at(index).GetDensityDatum();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DiscreteProbabilityVector::GetCumulativeProbability(const size_t index) const
{
    VerifyElementRequest(index);

    return m_discreteProbabilityData.at(index).GetCumulativeDatum();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void DiscreteProbabilityVector::GetAllAtIndex(const size_t index, float &x, float &probabilityDensity,
    float &cumulativeProbability) const
{
    VerifyElementRequest(index);

    const DiscreteProbabilityDatum &theDatum(m_discreteProbabilityData.at(index));
    x = theDatum.GetX();
    probabilityDensity = theDatum.GetDensityDatum();
    cumulativeProbability = theDatum.GetCumulativeDatum();

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline DiscreteProbabilityVector::DiscreteProbabilityDatum::DiscreteProbabilityDatum(const float &x,
    const float &densityDatum, const float &cumulativeDatum) :
    m_x(x),
    m_densityDatum(densityDatum),
    m_cumulativeDatum(cumulativeDatum)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DiscreteProbabilityVector::DiscreteProbabilityDatum::GetX() const
{
    return m_x;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DiscreteProbabilityVector::DiscreteProbabilityDatum::GetDensityDatum() const
{
    return m_densityDatum;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DiscreteProbabilityVector::DiscreteProbabilityDatum::GetCumulativeDatum() const
{
    return m_cumulativeDatum;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void DiscreteProbabilityVector::VerifyCompleteData() const
{
    if (2 > m_discreteProbabilityData.size())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

    if (m_discreteProbabilityData.back().GetX() > m_xUpperBound)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void DiscreteProbabilityVector::VerifyElementRequest(const size_t index) const
{
    if (GetSize() < index || 0 > index)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_OUT_OF_RANGE);

    return;
}

} // namespace lar_content

#endif // £ifndef  LAR_DISCRETE_PROBABILITY_VECTOR_H

