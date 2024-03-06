/**
 *  @file   larpandoracontent/LArObjects/LArTrackTwoViewOverlapResult.cc
 *
 *  @brief  Implementation of the lar track two view overlap result class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArObjects/LArTrackTwoViewOverlapResult.h"
#include "Objects/Cluster.h"

using namespace pandora;

namespace lar_content
{

TwoViewDeltaRayOverlapResult::TwoViewDeltaRayOverlapResult() :
    m_isInitialized(false),
    m_xOverlap(TwoViewXOverlap(0.f, 0.f, 0.f, 0.f)),
    m_commonMuonPfoList(),
    m_pBestMatchedCluster(nullptr),
    m_matchedClusterList(),
    m_reducedChiSquared(std::numeric_limits<float>::max())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TwoViewDeltaRayOverlapResult::TwoViewDeltaRayOverlapResult(const TwoViewXOverlap &xOverlap, const PfoList &commonMuonPfoList,
    const Cluster *const pBestMatchedCluster, const ClusterList &matchedClusterList, const float reducedChiSquared) :
    m_isInitialized(true),
    m_xOverlap(xOverlap),
    m_commonMuonPfoList(commonMuonPfoList),
    m_pBestMatchedCluster(pBestMatchedCluster),
    m_matchedClusterList(matchedClusterList),
    m_reducedChiSquared(reducedChiSquared)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TwoViewDeltaRayOverlapResult::TwoViewDeltaRayOverlapResult(const TwoViewDeltaRayOverlapResult &rhs) :
    m_isInitialized(rhs.m_isInitialized),
    m_xOverlap(rhs.GetXOverlap()),
    m_commonMuonPfoList(rhs.GetCommonMuonPfoList()),
    m_pBestMatchedCluster(rhs.GetBestMatchedCluster()),
    m_matchedClusterList(rhs.GetMatchedClusterList()),
    m_reducedChiSquared(rhs.GetReducedChiSquared())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TwoViewDeltaRayOverlapResult::~TwoViewDeltaRayOverlapResult()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster *TwoViewDeltaRayOverlapResult::GetBestMatchedAvailableCluster() const
{
    unsigned int highestNHits(0);
    const Cluster *pBestMatchedCluster(nullptr);

    for (const Cluster *const pMatchedCluster : m_matchedClusterList)
    {
        if (!pMatchedCluster->IsAvailable())
            continue;

        if (pMatchedCluster->GetNCaloHits() > highestNHits)
        {
            highestNHits = pMatchedCluster->GetNCaloHits();
            pBestMatchedCluster = pMatchedCluster;
        }
    }

    return pBestMatchedCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

TwoViewDeltaRayOverlapResult &TwoViewDeltaRayOverlapResult::operator=(const TwoViewDeltaRayOverlapResult &rhs)
{
    m_isInitialized = rhs.m_isInitialized;
    m_xOverlap = rhs.GetXOverlap();
    m_commonMuonPfoList = rhs.GetCommonMuonPfoList();
    m_pBestMatchedCluster = rhs.GetBestMatchedCluster();
    m_matchedClusterList = rhs.GetMatchedClusterList();
    m_reducedChiSquared = rhs.GetReducedChiSquared();

    return *this;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewDeltaRayOverlapResult::operator<(const TwoViewDeltaRayOverlapResult &rhs) const
{
    if (this == &rhs)
        return false;

    if (!m_isInitialized && !rhs.m_isInitialized)
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    if (!m_isInitialized)
        return true;

    if (!rhs.m_isInitialized)
        return false;

    return (m_xOverlap.GetTwoViewXOverlapSpan() < rhs.m_xOverlap.GetTwoViewXOverlapSpan());
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

TrackTwoViewOverlapResult::TrackTwoViewOverlapResult() : m_isInitialized(false), m_matchingScore(0)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackTwoViewOverlapResult::TrackTwoViewOverlapResult(const float matchingScore) : m_isInitialized(true), m_matchingScore(matchingScore)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackTwoViewOverlapResult::TrackTwoViewOverlapResult(const TrackTwoViewOverlapResult &rhs) :
    m_isInitialized(rhs.m_isInitialized), m_matchingScore(rhs.m_matchingScore)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackTwoViewOverlapResult::~TrackTwoViewOverlapResult()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackTwoViewOverlapResult::operator<(const TrackTwoViewOverlapResult &rhs) const
{
    if (this == &rhs)
        return false;

    if (!m_isInitialized && !rhs.m_isInitialized)
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    if (!m_isInitialized)
        return true;

    if (!rhs.m_isInitialized)
        return false;

    return (m_matchingScore < rhs.m_matchingScore);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackTwoViewOverlapResult::operator>(const TrackTwoViewOverlapResult &rhs) const
{
    if (this == &rhs)
        return false;

    return !(*this < rhs);
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackTwoViewOverlapResult &TrackTwoViewOverlapResult::operator=(const TrackTwoViewOverlapResult &rhs)
{
    if (this == &rhs)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    m_isInitialized = rhs.m_isInitialized;
    m_matchingScore = rhs.m_matchingScore;

    return *this;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

TwoViewTransverseOverlapResult::TwoViewTransverseOverlapResult() :
    TrackTwoViewOverlapResult(),
    m_downsamplingFactor(0.f),
    m_nSamplingPoints(0),
    m_nMatchedSamplingPoints(0),
    m_correlationCoefficient(0.f),
    m_twoViewXOverlap(TwoViewXOverlap(0.f, 0.f, 0.f, 0.f))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TwoViewTransverseOverlapResult::TwoViewTransverseOverlapResult(const float matchingScore, const float downsamplingFactor, const unsigned int nSamplingPoints,
    const unsigned int nMatchedSamplingPoints, const float correlationCoefficient, const TwoViewXOverlap &twoViewXOverlap) :
    TrackTwoViewOverlapResult(matchingScore),
    m_downsamplingFactor(downsamplingFactor),
    m_nSamplingPoints(nSamplingPoints),
    m_nMatchedSamplingPoints(nMatchedSamplingPoints),
    m_correlationCoefficient(correlationCoefficient),
    m_twoViewXOverlap(twoViewXOverlap)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TwoViewTransverseOverlapResult::TwoViewTransverseOverlapResult(const TwoViewTransverseOverlapResult &rhs) :
    TrackTwoViewOverlapResult(rhs),
    m_downsamplingFactor(rhs.m_downsamplingFactor),
    m_nSamplingPoints(rhs.m_nSamplingPoints),
    m_nMatchedSamplingPoints(rhs.m_nMatchedSamplingPoints),
    m_correlationCoefficient(rhs.m_correlationCoefficient),
    m_twoViewXOverlap(rhs.m_isInitialized ? rhs.m_twoViewXOverlap : TwoViewXOverlap(0.f, 0.f, 0.f, 0.f))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TwoViewTransverseOverlapResult::~TwoViewTransverseOverlapResult()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewTransverseOverlapResult::operator<(const TwoViewTransverseOverlapResult &rhs) const
{
    if (this == &rhs)
        return false;

    if (!m_isInitialized && !rhs.m_isInitialized)
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    if (!m_isInitialized)
        return true;

    if (!rhs.m_isInitialized)
        return false;

    if (std::fabs(m_matchingScore - rhs.m_matchingScore) > std::numeric_limits<float>::epsilon())
        return (m_matchingScore < rhs.m_matchingScore);

    if (std::fabs(m_correlationCoefficient - rhs.m_correlationCoefficient) > std::numeric_limits<float>::epsilon())
        return (m_correlationCoefficient < rhs.m_correlationCoefficient);

    if (m_nMatchedSamplingPoints != rhs.m_nMatchedSamplingPoints)
        return (m_nMatchedSamplingPoints < rhs.m_nMatchedSamplingPoints);

    if (m_nSamplingPoints != rhs.m_nSamplingPoints)
        return (m_nSamplingPoints < rhs.m_nSamplingPoints);

    if (std::fabs(this->GetLocallyMatchedFraction() - rhs.GetLocallyMatchedFraction()) > std::numeric_limits<float>::epsilon())
        return (this->GetLocallyMatchedFraction() < rhs.GetLocallyMatchedFraction());

    if (std::fabs(m_twoViewXOverlap.GetTwoViewXOverlapSpan() - rhs.m_twoViewXOverlap.GetTwoViewXOverlapSpan()) > std::numeric_limits<float>::epsilon())
        return (m_twoViewXOverlap.GetTwoViewXOverlapSpan() < rhs.m_twoViewXOverlap.GetTwoViewXOverlapSpan());

    if (std::fabs(m_twoViewXOverlap.GetXSpan0() - rhs.m_twoViewXOverlap.GetXSpan0()) > std::numeric_limits<float>::epsilon())
        return (m_twoViewXOverlap.GetXSpan0() < rhs.m_twoViewXOverlap.GetXSpan0());

    return (m_twoViewXOverlap.GetXSpan1() < rhs.m_twoViewXOverlap.GetXSpan1());
}

//------------------------------------------------------------------------------------------------------------------------------------------

TwoViewTransverseOverlapResult &TwoViewTransverseOverlapResult::operator=(const TwoViewTransverseOverlapResult &rhs)
{
    this->TrackTwoViewOverlapResult::operator=(rhs);

    if (rhs.m_isInitialized)
    {
        m_isInitialized = rhs.m_isInitialized;
        m_matchingScore = rhs.m_matchingScore;
        m_downsamplingFactor = rhs.m_downsamplingFactor;
        m_nSamplingPoints = rhs.m_nSamplingPoints;
        m_nMatchedSamplingPoints = rhs.m_nMatchedSamplingPoints;
        m_correlationCoefficient = rhs.m_correlationCoefficient;
        m_twoViewXOverlap = rhs.m_twoViewXOverlap;
    }
    else
    {
        m_isInitialized = false;
        m_matchingScore = 0.f;
        m_downsamplingFactor = 0.f;
        m_nSamplingPoints = 0;
        m_nMatchedSamplingPoints = 0;
        m_correlationCoefficient = 0.f;
        m_twoViewXOverlap = TwoViewXOverlap(0.f, 0.f, 0.f, 0.f);
    }

    return *this;
}

} // namespace lar_content
