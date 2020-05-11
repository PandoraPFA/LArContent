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

namespace lar_content
{

/**
 *  @brief  DiscreteProbabilityVector class
 */
class DiscreteProbabilityVector
{
public:

    template<typename T1, typename T2>
    using InputData = std::vector<std::pair<T1, T2> >;

    /**
     *  @brief  Constructor
     *
     *  @param  param description
     */
    template <typename T1, typename T2>
    DiscreteProbabilityVector(InputData<T1, T2> inputData);

    /**
     *  @brief  Store data to later convert to a cumulative distribution
     *
     *  @param  x independent variable
     *  @param  y dependent variable
     */
    template <typename T1, typename T2>
    void CollectInputData(T1 x, T2 y);

    template <typename T1, typename T2>
    void CollectCumulativeData(T1 x, T2 y);


    /**
     *  @brief  Get size of input data vector
     *
     *  @return input data vector size
     */
    size_t GetInputVectorSize() const;

    /**
     *  @brief  Get X value for specific element from input vector
     *
     *  @param index the index in the DCD
     *
     *  @return X value for element in DCD
     */
    float GetXFromInputVector(int index) const;

    /**
     *  @brief  Get size of the DCD
     *
     *  @return size of DCD
     */
    size_t GetSize() const;

    /**
     *  @brief  Get X value for specific element from input vector
     *
     *  @param index the index in the input vector
     *  @param x the x-value for the element in DCD
     *  @param y the y-value for the element in DCD
     */
    void GetXandY(int index, float &x, float &y) const;

    /**
     *  @brief  Create the cumulative distribution from the collected data
     */
    void CreateCumulativeDistribution();


private:

    typedef std::vector<std::pair<float, float> > InputDataVector;
    typedef InputDataVector CumulativeDataVector;

    InputDataVector m_InputDataHolder;                 ///< Input data storage which gets convered into a cumulative distribution
    CumulativeDataVector m_CumulativeDataVector;       ///< Hold the cumulative distribution
    //float       m_uMinX;                        ///< The min x value in the u view
};

template <typename T1, typename T2>
inline DiscreteProbabilityVector::DiscreteProbabilityVector(InputData<T1, T2> inputData)
{
    std::cout<<"Input: 0 " << inputData[0].first << "  " << inputData[1].second << std::endl;
}

template <typename T1, typename T2>
inline void DiscreteProbabilityVector::CollectInputData(T1 x, T2 y)
{
    float X = static_cast<float>(x);
    float Y = static_cast<float>(y);

    m_InputDataHolder.emplace_back(X,Y);
}

template <typename T1, typename T2>
inline void DiscreteProbabilityVector::CollectCumulativeData(T1 x, T2 y)
{
    float X = static_cast<float>(x);
    float Y = static_cast<float>(y);

    m_CumulativeDataVector.emplace_back(X,Y);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline size_t DiscreteProbabilityVector::GetInputVectorSize() const
{
    return m_InputDataHolder.size();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DiscreteProbabilityVector::GetXFromInputVector(int index) const
{
    return m_InputDataHolder.at(index).first;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline size_t DiscreteProbabilityVector::GetSize() const
{
    return m_CumulativeDataVector.size();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void DiscreteProbabilityVector::GetXandY(int index, float &x, float &y) const
{
    x = m_CumulativeDataVector.at(index).first;
    y = m_CumulativeDataVector.at(index).second;

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void DiscreteProbabilityVector::CreateCumulativeDistribution()
{
    //float yLast(0.);
    //if (1 < m_InputDataHolder.size())
    //{
    //    float xStart(m_InputDataHolder[0].first);
    //    float xSecond(m_InputDataHolder[1].first);
    //    m_InputDataHolder.insert(m_InputDataHolder.begin(),std::pair<float,float>(xStart-(xSecond-xStart), 0.f));

    //    float xLast(m_InputDataHolder[m_InputDataHolder.size()-1].first);
    //    float xSecondLast(m_InputDataHolder[m_InputDataHolder.size()-2].first);
    //    CollectInputData(xLast + (xLast-xSecondLast), 0.f);

    //}
    for (size_t iData = 0; iData < m_InputDataHolder.size(); iData++)
    {
        float x = m_InputDataHolder[iData].first;
        float y = m_InputDataHolder[iData].second;
        if (0 != m_CumulativeDataVector.size()) y += m_CumulativeDataVector.back().second;
        m_CumulativeDataVector.emplace_back(x,y);
        //yLast+=y;
    }
    if (0 != m_CumulativeDataVector.size())
    {
        float yLast = m_CumulativeDataVector.back().second;
        for (size_t iData = 0; iData < m_CumulativeDataVector.size(); ++iData)
        {
            m_CumulativeDataVector[iData].second /= yLast;
        }
    }
    //std::cout<<"yLast: " << yLast << std::endl;
    //std::cout<<"In size: " << m_InputDataHolder.size() << "  out size: " << m_CumulativeDataVector.size() << std::endl;
}

/**
 *  @brief  x overlap result + operator
 *
 *  @param  lhs the first x overlap result to add
 *  @param  rhs the second x overlap result to add
 */
//DiscreteProbabilityVector operator+(const DiscreteProbabilityVector &lhs, const DiscreteProbabilityVector &rhs);

/*


inline DiscreteProbabilityVector::DiscreteProbabilityVector(const float uMinX, const float uMaxX, const float vMinX, const float vMaxX,
        const float wMinX, const float wMaxX, const float xOverlapSpan) :
    m_uMinX(uMinX),
    m_uMaxX(uMaxX),
    m_vMinX(vMinX),
    m_vMaxX(vMaxX),
    m_wMinX(wMinX),
    m_wMaxX(wMaxX),
    m_xOverlapSpan(xOverlapSpan)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DiscreteProbabilityVector::GetUMinX() const
{
    return m_uMinX;
}
*/



} // namespace lar_content

#endif // £ifndef  LAR_DISCRETE_CUMULATIVE_DISTRIBUTION_H

