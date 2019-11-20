/**
 *  @file   larpandoracontent/LArObjects/LArDiscreteCumulativeDistribution.h
 *
 *  @brief  Header file for the lar discrete cumulative distribution class.
 *
 *  $Log: $
 */
#ifndef LAR_DISCRETE_CUMULATIVE_DISTRIBUTION_H
#define LAR_DISCRETE_CUMULATIVE_DISTRIBUTION_H 1

namespace lar_content
{

/**
 *  @brief  DiscreteCumulativeDistribution class
 */
class DiscreteCumulativeDistribution
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  param description
     */
    //DiscreteCumulativeDistribution();

    /**
     *  @brief  Store data to later convert to a cumulative distribution
     *
     *  @param  x independent variable
     *  @param  y dependent variable
     */
    template <typename T1, T2>
    void CollectInputData(T1 x, T2 y) const;

    /**
     *  @brief  Get size of input data vector
     *
     *  @return input data vector size
     */
    size_t GetInputDataSize() const;



private:

    typedef std::vector<std::pair<float, float> > InputDataVector;
    typedef InputDataVector CumulativeDataVector;

    InputDataVector m_InputDataHolder;                 ///< Input data storage which gets convered into a cumulative distribution
    CumulativeDataVector m_CumulativeDataVector;       ///< Hold the cumulative distribution
    //float       m_uMinX;                        ///< The min x value in the u view
};

template <typename T1, T2>
inline void DiscreteCumulativeDistribution::CollectInputData(T1 x, T2 y) const
{
    float X = static_cast<float>(x);
    float Y = static_cast<float>(y);

    m_InputDataHolder.emplace_back(x,y);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline size_t DiscreteCumulativeDistribution::GetInputDataSize() const
{
    return m_InputDataHolder.size();
}

/**
 *  @brief  x overlap result + operator
 *
 *  @param  lhs the first x overlap result to add
 *  @param  rhs the second x overlap result to add
 */
//DiscreteCumulativeDistribution operator+(const DiscreteCumulativeDistribution &lhs, const DiscreteCumulativeDistribution &rhs);

/*


inline DiscreteCumulativeDistribution::DiscreteCumulativeDistribution(const float uMinX, const float uMaxX, const float vMinX, const float vMaxX,
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

inline float DiscreteCumulativeDistribution::GetUMinX() const
{
    return m_uMinX;
}
*/



} // namespace lar_content

#endif // £ifndef  LAR_DISCRETE_CUMULATIVE_DISTRIBUTION_H

