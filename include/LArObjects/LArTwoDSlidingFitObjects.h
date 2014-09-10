/**
 *  @file   LArContent/include/LArObjects/LArTwoDSlidingFitObjects.h
 *
 *  @brief  Header file for the lar two dimensional sliding fit objects.
 *
 *  $Log: $
 */
#ifndef LAR_TWO_D_SLIDING_FIT_OBJECTS_H
#define LAR_TWO_D_SLIDING_FIT_OBJECTS_H 1

#include <cmath>
#include <map>
#include <vector>

namespace lar_content
{

/**
 *  @brief  TransverseDirection enum
 */
enum TransverseDirection
{
    POSITIVE_IN_X,
    NEGATIVE_IN_X,
    UNCHANGED_IN_X,
    UNKNOWN
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  class LayerFitResult
 */
class LayerFitResult
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  l the l coordinate
     *  @param  fitT the fitted t coordinate
     *  @param  gradient the fitted gradient dt/dl
     *  @param  rms the rms of the fit residuals
     */
    LayerFitResult(const double l, const double fitT, const double gradient, const double rms);

    /**
     *  @brief  Get the l coordinate
     *
     *  @return the l coordinate
     */
    double GetL() const;

    /**
     *  @brief  Get the fitted t coordinate
     *
     *  @return the fitted t coordinate
     */
    double GetFitT() const;

    /**
     *  @brief  Get the fitted gradient dt/dz
     *
     *  @return the fitted gradient dt/dl
     */
    double GetGradient() const;

    /**
     *  @brief  Get the rms of the fit residuals
     *
     *  @return the rms of the fit residuals
     */
    double GetRms() const;

private:
    double          m_l;                                    ///< The l coordinate
    double          m_fitT;                                 ///< The fitted t coordinate
    double          m_gradient;                             ///< The fitted gradient dt/dl
    double          m_rms;                                  ///< The rms of the fit residuals
};

typedef std::map<int, LayerFitResult> LayerFitResultMap;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  LayerFitContribution class
 */
class LayerFitContribution
{
public:
    /**
     *  @brief  Default constructor
     */
    LayerFitContribution();

    /**
     *  @brief  Add point to layer fit
     *
     *  @param  l the longitudinal coordinate
     *  @param  t the transverse coordinate
     */
    void AddPoint(const float l, const float t);

    /**
     *  @brief  Get the sum t
     *
     *  @return the sum t
     */
    double GetSumT() const;

    /**
     *  @brief  Get the sum l
     *
     *  @return the sum l
     */
    double GetSumL() const;

    /**
     *  @brief  Get the sum t * t
     *
     *  @return the sum t * t
     */
    double GetSumTT() const;

    /**
     *  @brief  Get the sum l * t
     *
     *  @return the sum l * t
     */
    double GetSumLT() const;

    /**
     *  @brief  Get the sum l * l
     *
     *  @return the sum z * z
     */
    double GetSumLL() const;

    /**
     *  @brief  Get the number of points used
     *
     *  @return the number of points used
     */
    unsigned int GetNPoints() const;

private:
    double          m_sumT;                                 ///< The sum t
    double          m_sumL;                                 ///< The sum l
    double          m_sumTT;                                ///< The sum t * t
    double          m_sumLT;                                ///< The sum l * t
    double          m_sumLL;                                ///< The sum l * l
    unsigned int    m_nPoints;                              ///< The number of points used
};

typedef std::map<int, LayerFitContribution> LayerFitContributionMap;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  LayerInterpolation class
 */
class LayerInterpolation
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param firstayerIter  the iterator for the upstream layer
     *  @param secondLayerIter  the iterator for the downstream layer
     *  @param firstWeight  the weight to be applied to the upstream layer
     *  @param secondWeight  the weight to be applied to the downstream layer
     */
    LayerInterpolation(const LayerFitResultMap::const_iterator &firstLayerIter, const LayerFitResultMap::const_iterator &secondLayerIter,
        const float firstWeight, const float secondWeight);

    /**
     *  @brief  Get the start layer iterator
     *
     *  @return the iterator for the start layer
     */
    LayerFitResultMap::const_iterator GetStartLayerIter() const;

    /**
     *  @brief  Get the end layer iterator
     *
     *  @return the iterator for the end layer
     */
    LayerFitResultMap::const_iterator GetEndLayerIter() const;

    /**
     *  @brief  Get the start layer weight
    *
     *  @return the weight for the start layer
     */
    float GetStartLayerWeight() const;

    /**
     *  @brief  Get the end layer weight
     *
     *  @return the weight for the end layer
     */
    float GetEndLayerWeight() const;

private:
    LayerFitResultMap::const_iterator m_startLayerIter;     ///< The start layer iterator
    LayerFitResultMap::const_iterator m_endLayerIter;       ///< The end layer iterator
    float                             m_startLayerWeight;   ///< The start layer weight
    float                             m_endLayerWeight;     ///< The end layer weight
};

typedef std::vector<LayerInterpolation> LayerInterpolationList;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  FitSegment class
 */
class FitSegment
{
    public:
    /**
     *  @brief Constructor
     *
     *  @param startLayer the start layer
     *  @param endLayer the end layer
     *  @param startX the x position at the start layer
     *  @param endX the x position at the end layer
     */
    FitSegment(const int startLayer, const int endLayer, const float startX, const float endX);

    /**
     *  @brief  Get start layer
     *
     *  @return the start layer
     */
    int GetStartLayer() const;

    /**
     *  @brief  Get end layer
     *
     *  @return the end layer
     */
    int GetEndLayer() const;

    /**
     *  @brief  Get the minimum x value
     *
     *  @return the minimum x value
     */
    float GetMinX() const;

    /**
     *  @brief  Get the maximum x value
     *
     *  @return the maximum x value
     */
    float GetMaxX() const;

    /**
     *  @brief  Whether the x coordinate increases between the start and end layers
     *
     *  @return boolean
     */
    bool IsIncreasingX() const;

private:
    int             m_startLayer;                           ///< The start layer
    int             m_endLayer;                             ///< The end layer
    float           m_minX;                                 ///< The minimum x value
    float           m_maxX;                                 ///< The maximum x value
    bool            m_isIncreasingX;                        ///< Whether the x coordinate increases between the start and end layers
};

typedef std::vector<FitSegment> FitSegmentList;

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline LayerFitResult::LayerFitResult(const double l, const double fitT, const double gradient, const double rms) :
    m_l(l),
    m_fitT(fitT),
    m_gradient(gradient),
    m_rms(rms)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LayerFitResult::GetL() const
{
    return m_l;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LayerFitResult::GetFitT() const
{
    return m_fitT;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LayerFitResult::GetGradient() const
{
    return m_gradient;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LayerFitResult::GetRms() const
{
    return m_rms;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline LayerFitContribution::LayerFitContribution() :
    m_sumT(0.),
    m_sumL(0.),
    m_sumTT(0.),
    m_sumLT(0.),
    m_sumLL(0.),
    m_nPoints(0)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LayerFitContribution::AddPoint(const float l, const float t)
{
    const double T = static_cast<double>(t);
    const double L = static_cast<double>(l);

    m_sumT += T;
    m_sumL += L;
    m_sumTT += T * T;
    m_sumLT += L * T;
    m_sumLL += L * L;
    ++m_nPoints;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LayerFitContribution::GetSumT() const
{
    return m_sumT;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LayerFitContribution::GetSumL() const
{
    return m_sumL;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LayerFitContribution::GetSumLT() const
{
    return m_sumLT;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LayerFitContribution::GetSumLL() const
{
    return m_sumLL;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LayerFitContribution::GetSumTT() const
{
    return m_sumTT;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int LayerFitContribution::GetNPoints() const
{
    return m_nPoints;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline LayerInterpolation::LayerInterpolation(const LayerFitResultMap::const_iterator &startLayerIter,
        const LayerFitResultMap::const_iterator &endLayerIter, const float startLayerWeight, const float endLayerWeight) :
    m_startLayerIter(startLayerIter),
    m_endLayerIter(endLayerIter),
    m_startLayerWeight(startLayerWeight),
    m_endLayerWeight(endLayerWeight)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline LayerFitResultMap::const_iterator LayerInterpolation::GetStartLayerIter() const
{
    return m_startLayerIter;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline LayerFitResultMap::const_iterator LayerInterpolation::GetEndLayerIter() const
{
    return m_endLayerIter;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LayerInterpolation::GetStartLayerWeight() const
{
    return m_startLayerWeight;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LayerInterpolation::GetEndLayerWeight() const
{
    return m_endLayerWeight;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline FitSegment::FitSegment(const int startLayer, const int endLayer, const float startX, const float endX) :
    m_startLayer(startLayer),
    m_endLayer(endLayer)
{
    m_minX = std::min(startX, endX);
    m_maxX = std::max(startX, endX);
    m_isIncreasingX = (endX > startX);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int FitSegment::GetStartLayer() const
{
    return m_startLayer;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int FitSegment::GetEndLayer() const
{
    return m_endLayer;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float FitSegment::GetMinX() const
{
    return m_minX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float FitSegment::GetMaxX() const
{
    return m_maxX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool FitSegment::IsIncreasingX() const
{
    return m_isIncreasingX;
}

} // namespace lar_content

#endif // #ifndef LAR_TWO_D_SLIDING_FIT_OBJECTS_H
