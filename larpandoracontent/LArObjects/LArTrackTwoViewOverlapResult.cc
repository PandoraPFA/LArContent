/**
 *  @file   larpandoracontent/LArObjects/LArTrackTwoViewOverlapResult.cc
 *
 *  @brief  Implementation of the lar track overlap result class.
 *
 *  $Log: $
 */

#include "Objects/CaloHit.h"

#include "Pandora/PandoraInternal.h"

#include "larpandoracontent/LArObjects/LArTrackTwoViewOverlapResult.h"

using namespace pandora;

namespace lar_content
{

TrackTwoViewOverlapResult::TrackTwoViewOverlapResult() :
    m_isInitialized(false) 
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackTwoViewOverlapResult::TrackTwoViewOverlapResult(const float matchingScore) :
    m_isInitialized(true),
    m_matchingScore(matchingScore)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackTwoViewOverlapResult::TrackTwoViewOverlapResult(const TrackTwoViewOverlapResult &rhs) :
    m_isInitialized(rhs.m_isInitialized),
    m_matchingScore(rhs.m_matchingScore)
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

    if (!m_isInitialized && !rhs.IsInitialized())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    if (!m_isInitialized)
        return true;

    if (!rhs.IsInitialized())
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
    m_twoViewXOverlap(TwoViewXOverlap(0.f, 0.f, 0.f, 0.f, 0.f))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TwoViewTransverseOverlapResult::TwoViewTransverseOverlapResult(const float matchingScore, const TwoViewXOverlap &twoViewXOverlap) :
    TrackTwoViewOverlapResult(matchingScore),
    m_twoViewXOverlap(twoViewXOverlap)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TwoViewTransverseOverlapResult::TwoViewTransverseOverlapResult(const TwoViewTransverseOverlapResult &rhs) :
    TrackTwoViewOverlapResult(rhs),
    m_twoViewXOverlap(rhs.IsInitialized() ? rhs.GetTwoViewXOverlap() : TwoViewXOverlap(0.f, 0.f, 0.f, 0.f, 0.f))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TwoViewTransverseOverlapResult::~TwoViewTransverseOverlapResult()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TwoViewTransverseOverlapResult &TwoViewTransverseOverlapResult::operator=(const TwoViewTransverseOverlapResult &rhs)
{
    this->TrackTwoViewOverlapResult::operator=(rhs);

    if (rhs.IsInitialized())
    {
        m_isInitialized = rhs.m_isInitialized;
        m_matchingScore = rhs.m_matchingScore;
        m_twoViewXOverlap = rhs.GetTwoViewXOverlap();
    }
    else
    {
        m_matchingScore = 0.f;
        m_twoViewXOverlap = TwoViewXOverlap(0.f, 0.f, 0.f, 0.f, 0.f);
    }

    return *this;
}

} // namespace lar_content
