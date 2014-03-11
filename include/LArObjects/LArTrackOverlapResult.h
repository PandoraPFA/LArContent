/**
 *  @file   LArContent/include/LArObjects/LArTrackOverlapResult.h
 * 
 *  @brief  Header file for the lar track overlap result class.
 * 
 *  $Log: $
 */
#ifndef LAR_TRACK_OVERLAP_RESULT_H
#define LAR_TRACK_OVERLAP_RESULT_H 1

#include <vector>

namespace lar
{

/**
 *  @brief  TrackOverlapResult class
 */
class TrackOverlapResult
{
public:
    /**
     *  @brief  Default constructor
     */
    TrackOverlapResult();

    /**
     *  @brief  Constructor
     * 
     *  @param  nMatchedSamplingPoints
     *  @param  nSamplingPoints
     *  @param  chi2
     */
    TrackOverlapResult(const unsigned int nMatchedSamplingPoints, const unsigned int nSamplingPoints, const float chi2);

    /**
     *  @brief  Copy constructor
     * 
     *  @param  rhs
     */
    TrackOverlapResult(const TrackOverlapResult &rhs);

    /**
     *  @brief  Whether the track overlap result has been initialized
     *
     *  @return boolean
     */
    bool IsInitialized() const;

    /**
     *  @brief  Get the number of matched sampling points
     *
     *  @return the number of matched sampling points
     */
    unsigned int GetNMatchedSamplingPoints() const;

    /**
     *  @brief  Get the number of sampling points
     *
     *  @return the number of sampling points
     */
    unsigned int GetNSamplingPoints() const;

    /**
     *  @brief  Get the fraction of sampling points resulting in a match
     *
     *  @return the fraction of sampling points resulting in a match
     */
    float GetMatchedFraction() const;

    /**
     *  @brief  Get the absolute chi2 value
     *
     *  @return the absolute chi2 value
     */
    float GetChi2() const;

    /**
     *  @brief  Get the chi2 per samping point value
     *
     *  @return the chi2 per samping point value
     */
    float GetReducedChi2() const;

    /**
     *  @brief  Track overlap result less than operator
     * 
     *  @param  rhs the track overlap result for comparison
     */
    bool operator<(const TrackOverlapResult &rhs) const;

    /**
     *  @brief  Track overlap result greater than operator
     * 
     *  @param  rhs the track overlap result for comparison
     */
    bool operator>(const TrackOverlapResult &rhs) const;

    /**
     *  @brief  Track overlap result assigment operator
     * 
     *  @param  rhs the track overlap result to assign
     */
    TrackOverlapResult &operator=(const TrackOverlapResult &rhs);

private:
    bool            m_isInitialized;                ///< Whether the track overlap result has been initialized
    unsigned int    m_nMatchedSamplingPoints;       ///< The number of matched sampling points
    unsigned int    m_nSamplingPoints;              ///< The number of sampling points
    float           m_matchedFraction;              ///< The fraction of sampling points resulting in a match
    float           m_chi2;                         ///< The absolute chi2 value
    float           m_reducedChi2;                  ///< The chi2 per samping point value
};

typedef std::vector<TrackOverlapResult> TrackOverlapResultVector;

/**
 *  @brief  Track overlap result + operator
 * 
 *  @param  lhs the first track overlap result to add
 *  @param  rhs the second track overlap result to add
 */
TrackOverlapResult operator+(const TrackOverlapResult &lhs, const TrackOverlapResult &rhs);

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool TrackOverlapResult::IsInitialized() const
{
    return m_isInitialized;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int TrackOverlapResult::GetNMatchedSamplingPoints() const
{
    return m_nMatchedSamplingPoints;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int TrackOverlapResult::GetNSamplingPoints() const
{
    return m_nSamplingPoints;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackOverlapResult::GetMatchedFraction() const
{
    return m_matchedFraction;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackOverlapResult::GetChi2() const
{
    return m_chi2;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackOverlapResult::GetReducedChi2() const
{
    return m_reducedChi2;
}

} // namespace lar

#endif // #ifndef LAR_TRACK_OVERLAP_RESULT_H
