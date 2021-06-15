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

#include "larpandoracontent/LArObjects/LArTrackOverlapResult.h"
#include "larpandoracontent/LArObjects/LArTwoViewXOverlap.h"

#include <vector>

namespace lar_content
{
/**
 *  @brief  TwoViewDeltaRayOverlapResult class
 */
class TwoViewDeltaRayOverlapResult
{
public:
    /**
     *  @brief Default constructor
     */
    TwoViewDeltaRayOverlapResult();

    /**
     *  @brief  Constructor
     *
     *  @param  xOverlap the xOverlap of the two clusters
     *  @param  commonMuonPfoList the list of muon pfos whose clusters lie close to delta ray clusters in all views
     *  @param  pBestMatchedCluster the largest cluster that lies on the projected hits
     *  @param  matchedClusterList the list of clusters that lie on the projected hits
     *  @param  reducedChiSquared the reduced chi squared
     */
    TwoViewDeltaRayOverlapResult(const TwoViewXOverlap &xOverlap, const pandora::PfoList &commonMuonPfoList,
        const pandora::Cluster *const pBestMatchedCluster, const pandora::ClusterList &matchedClusterList, const float reducedChiSquared);

    /**
     *  @brief  Copy constructor
     *
     *  @param  rhs
     */
    TwoViewDeltaRayOverlapResult(const TwoViewDeltaRayOverlapResult &rhs);

    /**
     *  @brief  Destructor
     */
    virtual ~TwoViewDeltaRayOverlapResult();

    /**
     *  @brief  Whether the track overlap result has been initialized
     *
     *  @return boolean
     */
    bool IsInitialized() const;

    /**
     *  @brief  Get the x overlap object
     *
     *  @return the x overlap object
     */
    const TwoViewXOverlap &GetXOverlap() const;

    /**
     *  @brief  Get the best matched cluster
     *
     *  @return the best matched cluster
     */
    const pandora::Cluster *GetBestMatchedCluster() const;

    /**
     *  @brief  Get the common muon pfo list
     *
     *  @return the common muon pfo list
     */
    const pandora::PfoList &GetCommonMuonPfoList() const;

    /**
     *  @brief  Get the matched cluster list
     *
     *  @return the matched cluster list
     */
    const pandora::ClusterList &GetMatchedClusterList() const;

    /**
     *  @brief  Get the reduced chi squared value
     *
     *  @return the reduced chi squared value
     */
    float GetReducedChiSquared() const;

    /**
     *  @brief  Get the best matched available cluster
     *
     *  @return the best matched available cluster
     */
    const pandora::Cluster *GetBestMatchedAvailableCluster() const;

    /**
     *  @brief  Track overlap result assigment operator
     *
     *  @param  rhs the track overlap result to assign
     */
    TwoViewDeltaRayOverlapResult &operator=(const TwoViewDeltaRayOverlapResult &rhs);

    /**
     *  @brief  Track two view overlap result less than operator
     *
     *  @param  rhs the track two view overlap result for comparison
     */
    bool operator<(const TwoViewDeltaRayOverlapResult &rhs) const;

private:
    bool m_isInitialized;                          ///< Whether the two view delta ray overlap result has been initialized
    TwoViewXOverlap m_xOverlap;                    ///< The x overlap object
    pandora::PfoList m_commonMuonPfoList;          ///< The list of muon pfos whose clusters lie close to delta ray clusters in all views
    const pandora::Cluster *m_pBestMatchedCluster; ///< The largest cluster that lies on the projected hits
    pandora::ClusterList m_matchedClusterList;     ///< The list of clusters that lie on the projected hits
    float m_reducedChiSquared;                     ///< The reduced chi squared of the best matched cluster
};

//------------------------------------------------------------------------------------------------------------------------------------------

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
    bool m_isInitialized;  ///< Whether the track overlap result has been initialized
    float m_matchingScore; ///< The compatability score for the two objects associated with the overlap result
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
     *  @param  downsamplingFactor the downsampling factor
     *  @param  nSamplingPoints the number of sampling points used in the matching
     *  @param  nMatchedSamplingPoints the number of matched sampling points
     *  @param  correlationCoefficient the corerlation coefficient for the matching candidate
     *  @param  twoViewXOverlap the description of the geometrical overlap for the matching candidate
     */
    TwoViewTransverseOverlapResult(const float matchingScore, const float downsamplingFactor, const unsigned int nSamplingPoints,
        const unsigned int nMatchedSamplingPoints, const float correlationCoefficient, const TwoViewXOverlap &twoViewXOverlap);

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
     *  @brief  Get the number of matched sampling points
     *
     *  @return the number of matched sampling points
     */
    unsigned int GetNMatchedSamplingPoints() const;

    /**
     *  @brief  Get the number of re-upsampled sampling points
     *
     *  @return the number of re-upsampled sampling points
     */
    unsigned int GetNReUpsampledSamplingPoints() const;

    /**
     *  @brief  Get the number of matched re-upsampled sampling points
     *
     *  @return the number of matched re-upsampled sampling points
     */
    unsigned int GetNMatchedReUpsampledSamplingPoints() const;

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
     *  @brief  Track two view overlap result less than operator
     *
     *  @param  rhs the track two view overlap result for comparison
     */
    bool operator<(const TwoViewTransverseOverlapResult &rhs) const;

    /**
     *  @brief  Track overlap result assigment operator
     *
     *  @param  rhs the track overlap result to assign
     */
    TwoViewTransverseOverlapResult &operator=(const TwoViewTransverseOverlapResult &rhs);

private:
    float m_downsamplingFactor;            ///< The downsampling factor
    unsigned int m_nSamplingPoints;        ///< The number of sampling points
    unsigned int m_nMatchedSamplingPoints; ///< The number of matched sampling points
    float m_correlationCoefficient;        ///< The correlation coefficient
    TwoViewXOverlap m_twoViewXOverlap;     ///< The two view x overlap object
};

typedef std::vector<TwoViewTransverseOverlapResult> TwoViewTransverseOverlapResultVector;

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline bool TwoViewDeltaRayOverlapResult::IsInitialized() const
{
    return m_isInitialized;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const TwoViewXOverlap &TwoViewDeltaRayOverlapResult::GetXOverlap() const
{
    return m_xOverlap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Cluster *TwoViewDeltaRayOverlapResult::GetBestMatchedCluster() const
{
    return m_pBestMatchedCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::PfoList &TwoViewDeltaRayOverlapResult::GetCommonMuonPfoList() const
{
    return m_commonMuonPfoList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::ClusterList &TwoViewDeltaRayOverlapResult::GetMatchedClusterList() const
{
    return m_matchedClusterList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TwoViewDeltaRayOverlapResult::GetReducedChiSquared() const
{
    return m_reducedChiSquared;
}

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

inline unsigned int TwoViewTransverseOverlapResult::GetNMatchedSamplingPoints() const
{
    if (m_isInitialized)
        return m_nMatchedSamplingPoints;

    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int TwoViewTransverseOverlapResult::GetNReUpsampledSamplingPoints() const
{
    if (m_isInitialized)
        return static_cast<unsigned int>(m_downsamplingFactor * static_cast<float>(m_nSamplingPoints));

    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int TwoViewTransverseOverlapResult::GetNMatchedReUpsampledSamplingPoints() const
{
    if (m_isInitialized)
        return static_cast<unsigned int>(m_downsamplingFactor * static_cast<float>(m_nMatchedSamplingPoints));

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
        return (m_nSamplingPoints > 0 ? static_cast<float>(m_nMatchedSamplingPoints) / static_cast<float>(m_nSamplingPoints) : 0.f);

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
