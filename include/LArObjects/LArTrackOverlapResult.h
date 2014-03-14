/**
 *  @file   LArContent/include/LArObjects/LArTrackOverlapResult.h
 * 
 *  @brief  Header file for the lar track overlap result class.
 * 
 *  $Log: $
 */
#ifndef LAR_TRACK_OVERLAP_RESULT_H
#define LAR_TRACK_OVERLAP_RESULT_H 1

#include "Pandora/StatusCodes.h"

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
     *  @brief  Destructor
     */
    ~TrackOverlapResult();

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
     *  @brief  Track overlap result assigment operator
     * 
     *  @param  rhs the track overlap result to assign
     */
    TrackOverlapResult &operator=(const TrackOverlapResult &rhs);

protected:
    bool            m_isInitialized;                ///< Whether the track overlap result has been initialized
    unsigned int    m_nMatchedSamplingPoints;       ///< The number of matched sampling points
    unsigned int    m_nSamplingPoints;              ///< The number of sampling points
    float           m_matchedFraction;              ///< The fraction of sampling points resulting in a match
    float           m_chi2;                         ///< The absolute chi2 value
    float           m_reducedChi2;                  ///< The chi2 per samping point value
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  TransverseOverlapResult class
 */
class TransverseOverlapResult : public TrackOverlapResult
{
public:
    /**
     *  @brief  XOverlap class
     */
    class XOverlap
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  xSpanU the x span in the u view
         *  @param  xSpanV the x span in the v view
         *  @param  xSpanW the x span in the w view
         *  @param  xOverlapSpan the x overlap span
         */
        XOverlap(const float xSpanU, const float xSpanV, const float xSpanW, const float xOverlapSpan);

        /**
         *  @brief  Get the x span in the u view
         * 
         *  @return the x span in the u view
         */
        float GetXSpanU() const;

        /**
         *  @brief  Get the x span in the v view
         * 
         *  @return the x span in the v view
         */
        float GetXSpanV() const;

        /**
         *  @brief  Get the x span in the w view
         * 
         *  @return the x span in the w view
         */
        float GetXSpanW() const;

        /**
         *  @brief  Get the x overlap span
         * 
         *  @return the x overlap span
         */
        float GetXOverlapSpan() const;

    private:
        float       m_xSpanU;                       ///< The x span in the u view
        float       m_xSpanV;                       ///< The x span in the v view
        float       m_xSpanW;                       ///< The x span in the w view
        float       m_xOverlapSpan;                 ///< The x overlap span
    };

    /**
     *  @brief  Default constructor
     */
    TransverseOverlapResult();

    /**
     *  @brief  Constructor
     * 
     *  @param  nMatchedSamplingPoints
     *  @param  nSamplingPoints
     *  @param  chi2
     *  @param  xOverlap
     */
    TransverseOverlapResult(const unsigned int nMatchedSamplingPoints, const unsigned int nSamplingPoints, const float chi2,
        const XOverlap &xOverlap);

    /**
     *  @brief  Copy constructor
     * 
     *  @param  rhs
     */
    TransverseOverlapResult(const TransverseOverlapResult &rhs);

    /**
     *  @brief  Destructor
     */
    ~TransverseOverlapResult();

    /**
     *  @brief  Get the x overlap object
     * 
     *  @return the x overlap object
     */
    const XOverlap &GetXOverlap() const;

    /**
     *  @brief  Track overlap result less than operator
     * 
     *  @param  rhs the track overlap result for comparison
     */
    bool operator<(const TransverseOverlapResult &rhs) const;

    /**
     *  @brief  Track overlap result greater than operator
     * 
     *  @param  rhs the track overlap result for comparison
     */
    bool operator>(const TransverseOverlapResult &rhs) const;

    /**
     *  @brief  Track overlap result assigment operator
     * 
     *  @param  rhs the track overlap result to assign
     */
    TransverseOverlapResult &operator=(const TransverseOverlapResult &rhs);

private:
    XOverlap        m_xOverlap;                     ///< The x overlap object
};

typedef std::vector<TransverseOverlapResult> TransverseOverlapResultVector;

/**
 *  @brief  Transverse overlap result + operator
 * 
 *  @param  lhs the first transverse overlap result to add
 *  @param  rhs the second transverse overlap result to add
 */
TransverseOverlapResult operator+(const TransverseOverlapResult &lhs, const TransverseOverlapResult &rhs);

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline bool TrackOverlapResult::IsInitialized() const
{
    return m_isInitialized;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int TrackOverlapResult::GetNMatchedSamplingPoints() const
{
    if (m_isInitialized)
        return m_nMatchedSamplingPoints;

    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int TrackOverlapResult::GetNSamplingPoints() const
{
    if (m_isInitialized)
        return m_nSamplingPoints;

    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackOverlapResult::GetMatchedFraction() const
{
    if (m_isInitialized)
        return m_matchedFraction;

    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackOverlapResult::GetChi2() const
{
    if (m_isInitialized)
        return m_chi2;

    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackOverlapResult::GetReducedChi2() const
{
    if (m_isInitialized)
        return m_reducedChi2;

    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const TransverseOverlapResult::XOverlap &TransverseOverlapResult::GetXOverlap() const
{
    if (m_isInitialized)
        return m_xOverlap;

    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline TransverseOverlapResult::XOverlap::XOverlap(const float xSpanU, const float xSpanV, const float xSpanW, const float xOverlapSpan) :
    m_xSpanU(xSpanU),
    m_xSpanV(xSpanV),
    m_xSpanW(xSpanW),
    m_xOverlapSpan(xOverlapSpan)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TransverseOverlapResult::XOverlap::GetXSpanU() const
{
    return m_xSpanU;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TransverseOverlapResult::XOverlap::GetXSpanV() const
{
    return m_xSpanV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TransverseOverlapResult::XOverlap::GetXSpanW() const
{
    return m_xSpanW;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TransverseOverlapResult::XOverlap::GetXOverlapSpan() const
{
    return m_xOverlapSpan;
}

} // namespace lar

#endif // #ifndef LAR_TRACK_OVERLAP_RESULT_H
