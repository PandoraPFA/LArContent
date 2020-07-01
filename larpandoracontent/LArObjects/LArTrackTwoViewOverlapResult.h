/**
 *  @file   larpandoracontent/LArObjects/LArTrackTwoViewOverlapResult.h
 *
 *  @brief  Header file for the lar track two view overlap result class.
 *
 *  $Log: $
 */
#ifndef LAR_TRACK_TWO_VIEW_OVERLAP_RESULT_H
#define LAR_TRACK_TWO_VIEW_OVERLAP_RESULT_H 1

#include "Pandora/StatusCodes.h"

#include "larpandoracontent/LArObjects/LArTwoViewXOverlap.h"

#include <vector>

namespace lar_content
{

/**
 *  @brief  TrackTwoViewOverlapResult class
 */
class TrackTwoViewOverlapResult
{
public:
    /**
     *  @brief  Default constructor
     */
    TrackTwoViewOverlapResult();

    /**
     *  @brief  constructor
     *
     *  @param matchingScore
     */
    TrackTwoViewOverlapResult(const float matchingScore);

    /**
     *  @brief  Copy constructor
     *
     *  @param  rhs
     */
    TrackTwoViewOverlapResult(const TrackTwoViewOverlapResult &rhs);

    /**
     *  @brief  Destructor
     */
    virtual ~TrackTwoViewOverlapResult();

    /**
     *  @brief  Whether the track overlap result has been initialized
     *
     *  @return boolean
     */
    bool IsInitialized() const;

    /**
     *  @brief  Get the matching score of the overlap result
     *
     *  @return matching score
     */
    float GetMatchingScore() const;

    /**
     *  @brief  Track two view overlap result less than operator
     *
     *  @param  rhs the track two view overlap result for comparison
     */
    bool operator<(const TrackTwoViewOverlapResult &rhs) const;

    /**
     *  @brief  Track two view overlap result greater than operator
     *
     *  @param  rhs the track two view overlap result for comparison
     */
    bool operator>(const TrackTwoViewOverlapResult &rhs) const;

    /**
     *  @brief  Track overlap result assigment operator
     *
     *  @param  rhs the track overlap result to assign
     */
    TrackTwoViewOverlapResult &operator=(const TrackTwoViewOverlapResult &rhs);

protected:
    bool            m_isInitialized;                ///< Whether the track overlap result has been initialized
    float           m_matchingScore;                ///< The compatability score for the two objects associated with the overlap result
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  TwoViewTransverseOverlapResult class
 */
class TwoViewTransverseOverlapResult : public TrackTwoViewOverlapResult
{
public:
    /**
     *  @brief  Default constructor
     */
    TwoViewTransverseOverlapResult();

    /**
     *  @brief  Constructor
     *
     *  @param  matchingScore the matching candidate matching score
     *  @param  nSamplingPoints the number of sampling points used in the matching
     *  @param  correlationCoefficient the corerlation coefficient for the matching candidate
     *  @param locallyMatchedFraction the fraction of locally matching regions in the matching candidate
     *  @param  twoViewXOverlap the description of the geometrical overlap for the matching candidate
     */
    TwoViewTransverseOverlapResult(const float matchingScore, const unsigned int nSamplingPoints, 
        const float correlationCoefficient, const float locallyMatchedFraction, const TwoViewXOverlap &twoViewXOverlap);

    /**
     *  @brief  Copy constructor
     *
     *  @param  rhs the rhs
     */
    TwoViewTransverseOverlapResult(const TwoViewTransverseOverlapResult &rhs);

    /**
     *  @brief  Destructor
     */
    ~TwoViewTransverseOverlapResult();

    /**
     *  @brief  Get the number of sampling points
     *
     *  @return the number of sampling points
     */
    unsigned int GetNSamplingPoints() const;

    /**
     *  @brief  Get the correlation coefficient
     *
     *  @return the correlation coefficient
     */
    float GetCorrelationCoefficient() const;

    /**
     *  @brief  Get the locally matched fraction
     *
     *  @return the locally matched fraction
     */
    float GetLocallyMatchedFraction() const;

    /**
     *  @brief  Get the two view x overlap object
     *
     *  @return the two view x overlap object
     */
    const TwoViewXOverlap &GetTwoViewXOverlap() const;

    /**
     *  @brief  Track overlap result assigment operator
     *
     *  @param  rhs the track overlap result to assign
     */
    TwoViewTransverseOverlapResult &operator=(const TwoViewTransverseOverlapResult &rhs);

private:
    unsigned int           m_nSamplingPoints;                     ///< The number of sampling points
    float                  m_correlationCoefficient;              ///< The correlation coefficient
    float                  m_locallyMatchedFraction;              ///< The locally matched fraction
    TwoViewXOverlap        m_twoViewXOverlap;                     ///< The two view x overlap object
};

typedef std::vector<TwoViewTransverseOverlapResult> TwoViewTransverseOverlapResultVector;

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline bool TrackTwoViewOverlapResult::IsInitialized() const
{
    return m_isInitialized;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackTwoViewOverlapResult::GetMatchingScore() const
{
    if (m_isInitialized)
        return m_matchingScore;

    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int TwoViewTransverseOverlapResult::GetNSamplingPoints() const
{
    if (m_isInitialized)
        return m_nSamplingPoints;
   
    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TwoViewTransverseOverlapResult::GetCorrelationCoefficient() const
{
    if (m_isInitialized)
        return m_correlationCoefficient;
   
    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TwoViewTransverseOverlapResult::GetLocallyMatchedFraction() const
{
    if (m_isInitialized)
        return m_locallyMatchedFraction;
   
    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const TwoViewXOverlap &TwoViewTransverseOverlapResult::GetTwoViewXOverlap() const
{
    if (m_isInitialized)
        return m_twoViewXOverlap;

    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);
}

} // namespace lar_content

#endif // #ifndef LAR_TRACK_TWO_VIEW_OVERLAP_RESULT_H
