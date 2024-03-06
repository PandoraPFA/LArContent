/**
 *  @file   larpandoracontent/LArObjects/LArDiscreteProbabilityVector.h
 *
 *  @brief  Header file for the lar discrete probability vector class
 *
 *  $Log: $
 */
#ifndef LAR_DISCRETE_PROBABILITY_VECTOR_H
#define LAR_DISCRETE_PROBABILITY_VECTOR_H 1

#include "Pandora/PandoraInternal.h"
#include "Pandora/StatusCodes.h"

#include <limits>
#include <random>
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
    using InputData = std::vector<InputDatum<TX, TY>>;

    typedef InputData<float, float> AllFloatInputData;

    typedef pandora::FloatVector ResamplingPoints;

    /**
     *  @brief  Constructor
     *
     *  @param  inputData the data used to construct the probability vector
     *  @param  xUpperBound the upper bound of the probability vector
     *  @param  useWidths bool controlling whether the bin widths are used in calculations
     */
    template <typename TX, typename TY>
    DiscreteProbabilityVector(const InputData<TX, TY> &inputData, const TX xUpperBound, const bool useWidths);

    /**
     *  @brief  Constructor
     *
     *  @param  discreteProbabilityVector a discrete probability vector to resample from
     *  @param  resamplingPoints the points to resample the discrete probability vector with
     */
    DiscreteProbabilityVector(const DiscreteProbabilityVector &discreteProbabilityVector, const ResamplingPoints &resamplingPoints);

    /**
     *  @brief  Constructor
     *
     *  @param  discreteProbabilityVector a discrete probability vector to randomly rearrange
     *  @param  randomNumberGenerator the random number generator for the random reshuffling
     */
    DiscreteProbabilityVector(const DiscreteProbabilityVector &discreteProbabilityVector, std::mt19937 &randomNumberGenerator);

    /**
     *  @brief  Evaluate the cumulative probability at arbitrary x
     *
     *  @param  x the x value
     *
     *  @return the cumulative probability
     */
    float EvaluateCumulativeProbability(const float x) const;

    /**
     *  @brief  Get the size of the probability vector
     *
     *  @return the probability vector size
     */
    unsigned int GetSize() const;

    /**
     *  @brief  Get the x value of the element in the vector
     *
     *  @param  index the index in the vector
     *
     *  @return the x value
     */
    float GetX(const unsigned int index) const;

    /**
     *  @brief  Get the probability value of the element in the vector
     *
     *  @param  index the index in the vector
     *
     *  @return the probability
     */
    float GetProbability(const unsigned int index) const;

    /**
     *  @brief  Get the probability density value of the element in the vector
     *
     *  @param  index the index in the vector
     *
     *  @return the probability density
     */
    float GetProbabilityDensity(const unsigned int index) const;

    /**
     *  @brief  Get the cumulative probability value of the element in the vector
     *
     *  @param  index the index in the vector
     *
     *  @return the cumulative probability
     */
    float GetCumulativeProbability(const unsigned int index) const;

    /**
     *  @brief  Get the width of the element in the vectorr
     *
     *  @param  index the index in the vector
     *
     *  @return the width of the probability bin
     */
    float GetWidth(const unsigned int index) const;

    /**
     *  @brief  Get all information stored at a particular index
     *
     *  @param  index the index in the vector
     *  @param  x the x value
     *  @param  probabilityDensity the probability density value
     *  @param  cumulativeProbability the cumulative probability value
     *  @param  width the width of the probability bin
     */
    void GetAllAtIndex(const unsigned int index, float &x, float &probabilityDensity, float &cumulativeProbability, float &width) const;

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
         *  @param  width the width of the bin
         */
        DiscreteProbabilityDatum(const float x, const float densityDatum, const float cumulativeDatum, const float width);

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

        /**
         *  @brief  Get the width of the datum
         *
         *  @return the width
         */
        float GetWidth() const;

    private:
        float m_x;               ///< The x coordinate
        float m_densityDatum;    ///< The probability density value
        float m_cumulativeDatum; ///< The cumulative probability value
        float m_width;           ///< The width of the probability bin
    };

    typedef std::vector<DiscreteProbabilityDatum> DiscreteProbabilityData;

    /**
     *  @brief  Get a initialised probability data vector from the input data
     *
     *  @param  inputData the input data
     *
     *  @return a fully-initialised discrete probability data vector
     */
    template <typename TX, typename TY>
    DiscreteProbabilityData InitialiseDiscreteProbabilityData(InputData<TX, TY> inputData) const;

    /**
     *  @brief  Get a resampled probability data vector by resampling another probability data vector
     *
     *  @param  discreteProbabilityVector another discrete probability vector
     *  @param  resamplingPoints the points to resample the discrete probability vector with
     *
     *  @return a resampled probability data vector
     */
    DiscreteProbabilityData ResampleDiscreteProbabilityData(
        const DiscreteProbabilityVector &discreteProbabilityVector, const ResamplingPoints &resamplingPoints) const;

    /**
     *  @brief  Get a randomised probability data vector in which the x values are unchanged, the probability density is
     *          randomised and the cumulative probability is recalculated
     *
     *  @param  discreteProbabilityVector another discrete probability vector
     *  @param  randomNumberGenerator the random number generator for the random reshuffling
     *
     *  @return a resampled probability data vector
     */
    DiscreteProbabilityData RandomiseDiscreteProbabilityData(
        const DiscreteProbabilityVector &discreteProbabilityVector, std::mt19937 &randomNumberGenerator) const;

    /**
     *  @brief  Sort the input data according to their x value
     *
     *  @param  lhs the first InputDatum
     *  @param  lhs the second InputDatum
     *
     *  @return a bool dictacting swapping the two InputDatums
     */
    template <typename TX, typename TY>
    static bool SortInputDataByX(const InputDatum<TX, TY> &lhs, const InputDatum<TX, TY> &rhs);

    /**
     *  @brief  Calculate the probability normalisation
     *
     *  @param  inputData the input data
     *
     *  @return the probability normalisation
     */
    template <typename TX, typename TY>
    float CalculateNormalisation(const InputData<TX, TY> &inputData) const;

    /**
     *  @brief  Verify the integrity of the complete probability vector
     */
    void VerifyCompleteData() const;

    /**
     *  @brief  Verify the integrity of the element request
     *
     *  @param  index the index in the probability vector
     */
    void VerifyElementRequest(const unsigned int index) const;

    float m_xUpperBound;                               ///< the upper bound of the probability vector
    bool m_useWidths;                                  ///< controls whether bin widths are used in calculations
    DiscreteProbabilityData m_discreteProbabilityData; ///< the probability data
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int DiscreteProbabilityVector::GetSize() const
{
    return m_discreteProbabilityData.size();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DiscreteProbabilityVector::GetX(const unsigned int index) const
{
    this->VerifyElementRequest(index);

    return m_discreteProbabilityData.at(index).GetX();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DiscreteProbabilityVector::GetProbability(const unsigned int index) const
{
    this->VerifyElementRequest(index);

    return m_discreteProbabilityData.at(index).GetDensityDatum() * (m_useWidths ? m_discreteProbabilityData.at(index).GetWidth() : 1.f);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DiscreteProbabilityVector::GetProbabilityDensity(const unsigned int index) const
{
    this->VerifyElementRequest(index);

    return m_discreteProbabilityData.at(index).GetDensityDatum();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DiscreteProbabilityVector::GetCumulativeProbability(const unsigned int index) const
{
    this->VerifyElementRequest(index);

    return m_discreteProbabilityData.at(index).GetCumulativeDatum();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DiscreteProbabilityVector::GetWidth(const unsigned int index) const
{
    this->VerifyElementRequest(index);

    return m_discreteProbabilityData.at(index).GetWidth();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void DiscreteProbabilityVector::GetAllAtIndex(
    const unsigned int index, float &x, float &probabilityDensity, float &cumulativeProbability, float &width) const
{
    this->VerifyElementRequest(index);

    const DiscreteProbabilityDatum &theDatum(m_discreteProbabilityData.at(index));
    x = theDatum.GetX();
    probabilityDensity = theDatum.GetDensityDatum();
    cumulativeProbability = theDatum.GetCumulativeDatum();
    width = theDatum.GetWidth();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline DiscreteProbabilityVector::DiscreteProbabilityDatum::DiscreteProbabilityDatum(
    const float x, const float densityDatum, const float cumulativeDatum, const float width) :
    m_x(x), m_densityDatum(densityDatum), m_cumulativeDatum(cumulativeDatum), m_width(width)
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

inline float DiscreteProbabilityVector::DiscreteProbabilityDatum::GetWidth() const
{
    return m_width;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void DiscreteProbabilityVector::VerifyCompleteData() const
{
    if (2 > m_discreteProbabilityData.size())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

    if (m_discreteProbabilityData.back().GetX() - m_xUpperBound > std::numeric_limits<float>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void DiscreteProbabilityVector::VerifyElementRequest(const unsigned int index) const
{
    if (this->GetSize() < index)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_OUT_OF_RANGE);

    return;
}

} // namespace lar_content

#endif // #ifndef  LAR_DISCRETE_PROBABILITY_VECTOR_H
