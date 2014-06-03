/**
 *  @file   LArContent/include/LArObjects/LArTrackOverlapResult.h
 * 
 *  @brief  Header file for the lar track overlap result class.
 * 
 *  $Log: $
 */
#ifndef LAR_TRACK_OVERLAP_RESULT_H
#define LAR_TRACK_OVERLAP_RESULT_H 1

#include "Pandora/PandoraInputTypes.h"
#include "Pandora/StatusCodes.h"

#include <cmath>
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
         *  @param  uMinX min x value in the u view
         *  @param  uMaxX max x value in the u view
         *  @param  vMinX min x value in the v view
         *  @param  vMaxX max x value in the v view
         *  @param  wMinX min x value in the w view
         *  @param  wMaxX max x value in the w view
         *  @param  xOverlapSpan the x overlap span
         */
        XOverlap(const float uMinX, const float uMaxX, const float vMinX, const float vMaxX, const float wMinX, const float wMaxX, const float xOverlapSpan);

        /**
         *  @brief  Get the min x value in the u view
         * 
         *  @return the min x value in the u view
         */
        float GetUMinX() const;

        /**
         *  @brief  Get the max x value in the u view
         * 
         *  @return the max x value in the u view
         */
        float GetUMaxX() const;

        /**
         *  @brief  Get the min x value in the v view
         * 
         *  @return the min x value in the v view
         */
        float GetVMinX() const;

        /**
         *  @brief  Get the max x value in the v view
         * 
         *  @return the max x value in the v view
         */
        float GetVMaxX() const;

        /**
         *  @brief  Get the min x value in the w view
         * 
         *  @return the min x value in the w view
         */
        float GetWMinX() const;

        /**
         *  @brief  Get the max x value in the w view
         * 
         *  @return the max x value in the w view
         */
        float GetWMaxX() const;

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
        float       m_uMinX;                        ///< The min x value in the u view
        float       m_uMaxX;                        ///< The max x value in the u view
        float       m_vMinX;                        ///< The min x value in the v view
        float       m_vMaxX;                        ///< The max x value in the v view
        float       m_wMinX;                        ///< The min x value in the w view
        float       m_wMaxX;                        ///< The max x value in the w view
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

/**
 *  @brief  LongitudinalOverlapResult class
 */
class LongitudinalOverlapResult : public TrackOverlapResult
{
public:
    /**
     *  @brief  Default constructor
     */
    LongitudinalOverlapResult();

    /**
     *  @brief  Constructor
     * 
     *  @param  trackOverlapResult
     *  @param  innerChi2
     *  @param  outerChi2
     */
    LongitudinalOverlapResult(const TrackOverlapResult trackOverlapResult, const float innerChi2, const float outerChi2);

    /**
     *  @brief  Constructor
     * 
     *  @param  nMatchedSamplingPoints
     *  @param  nSamplingPoints
     *  @param  chi2
     *  @param  innerChi2
     *  @param  outerChi2
     */
    LongitudinalOverlapResult(const unsigned int nMatchedSamplingPoints, const unsigned int nSamplingPoints, const float chi2,
        const float innerChi2, const float outerChi2);

    /**
     *  @brief  Copy constructor
     * 
     *  @param  rhs
     */
    LongitudinalOverlapResult(const LongitudinalOverlapResult &rhs);

    /**
     *  @brief  Destructor
     */
    ~LongitudinalOverlapResult();

    /**
     *  @brief  
     * 
     *  @return 
     */
    float GetInnerChi2() const;

    /**
     *  @brief  
     * 
     *  @return 
     */
    float GetOuterChi2() const;

    /**
     *  @brief  Track overlap result assigment operator
     * 
     *  @param  rhs the track overlap result to assign
     */
    LongitudinalOverlapResult &operator=(const LongitudinalOverlapResult &rhs);

private:
    float       m_innerChi2;                        ///< The inner chi squared
    float       m_outerChi2;                        ///< The outer chi squared
};

typedef std::vector<LongitudinalOverlapResult> LongitudinalOverlapResultVector;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  FragmentOverlapResult class
 */
class FragmentOverlapResult : public TrackOverlapResult
{
public:
    /**
     *  @brief  Default constructor
     */
    FragmentOverlapResult();

    /**
     *  @brief  Constructor
     * 
     *  @param  trackOverlapResult
     *  @param  caloHitList
     *  @param  clusterList
     */
    FragmentOverlapResult(const TrackOverlapResult trackOverlapResult, const pandora::CaloHitList &caloHitList,
        const pandora::ClusterList &clusterList);

    /**
     *  @brief  Constructor
     * 
     *  @param  nMatchedSamplingPoints
     *  @param  nSamplingPoints
     *  @param  chi2
     *  @param  caloHitList
     *  @param  clusterList
     */
    FragmentOverlapResult(const unsigned int nMatchedSamplingPoints, const unsigned int nSamplingPoints, const float chi2,
        const pandora::CaloHitList &caloHitList, const pandora::ClusterList &clusterList);

    /**
     *  @brief  Copy constructor
     * 
     *  @param  rhs
     */
    FragmentOverlapResult(const FragmentOverlapResult &rhs);

    /**
     *  @brief  Destructor
     */
    ~FragmentOverlapResult();

    /**
     *  @brief  Get the list of fragment-associated hits
     * 
     *  @return the list of fragment-associated hits
     */
    const pandora::CaloHitList &GetFragmentCaloHitList() const;

    /**
     *  @brief  Get the list of fragment-associated clusters
     * 
     *  @return the list of fragment-associated clusters
     */
    const pandora::ClusterList &GetFragmentClusterList() const;

    /**
     *  @brief  Get the fragment hit type
     * 
     *  @return the fragment hit type
     */
    pandora::HitType GetFragmentHitType() const;

    /**
     *  @brief  Fragments overlap result assigment operator
     * 
     *  @param  rhs the track overlap result to assign
     */
    FragmentOverlapResult &operator=(const FragmentOverlapResult &rhs);

private:
    pandora::CaloHitList    m_caloHitList;      ///< The list of fragment-associated hits
    pandora::ClusterList    m_clusterList;      ///< The list of fragment-associated clusters
};

typedef std::vector<FragmentOverlapResult> FragmentOverlapResultVector;

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

inline TransverseOverlapResult::XOverlap::XOverlap(const float uMinX, const float uMaxX, const float vMinX, const float vMaxX,
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

inline float TransverseOverlapResult::XOverlap::GetUMinX() const
{
    return m_uMinX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TransverseOverlapResult::XOverlap::GetUMaxX() const
{
    return m_uMaxX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TransverseOverlapResult::XOverlap::GetVMinX() const
{
    return m_vMinX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TransverseOverlapResult::XOverlap::GetVMaxX() const
{
    return m_vMaxX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TransverseOverlapResult::XOverlap::GetWMinX() const
{
    return m_wMinX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TransverseOverlapResult::XOverlap::GetWMaxX() const
{
    return m_wMaxX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TransverseOverlapResult::XOverlap::GetXSpanU() const
{
    return std::fabs(m_uMaxX - m_uMinX);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TransverseOverlapResult::XOverlap::GetXSpanV() const
{
    return std::fabs(m_vMaxX - m_vMinX);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TransverseOverlapResult::XOverlap::GetXSpanW() const
{
    return std::fabs(m_wMaxX - m_wMinX);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TransverseOverlapResult::XOverlap::GetXOverlapSpan() const
{
    return m_xOverlapSpan;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline float LongitudinalOverlapResult::GetInnerChi2() const
{
    return m_innerChi2;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LongitudinalOverlapResult::GetOuterChi2() const
{
    return m_outerChi2;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CaloHitList &FragmentOverlapResult::GetFragmentCaloHitList() const
{
    return m_caloHitList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::ClusterList &FragmentOverlapResult::GetFragmentClusterList() const
{
    return m_clusterList;
}

} // namespace lar

#endif // #ifndef LAR_TRACK_OVERLAP_RESULT_H
