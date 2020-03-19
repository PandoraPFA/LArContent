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
    m_isInitialized(true) //TODO set to false when alternative constructor is written
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackTwoViewOverlapResult::TrackTwoViewOverlapResult(const TrackTwoViewOverlapResult &rhs) :
    m_isInitialized(rhs.m_isInitialized)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackTwoViewOverlapResult::~TrackTwoViewOverlapResult()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackTwoViewOverlapResult &TrackTwoViewOverlapResult::operator=(const TrackTwoViewOverlapResult &rhs)
{
    if (this == &rhs)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    m_isInitialized = rhs.m_isInitialized;

    return *this;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

TwoViewTransverseOverlapResult::TwoViewTransverseOverlapResult() :
    TrackTwoViewOverlapResult(),
    m_twoViewXOverlap(TwoViewXOverlap(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TwoViewTransverseOverlapResult::TwoViewTransverseOverlapResult(const TwoViewXOverlap &twoViewXOverlap) :
    TrackTwoViewOverlapResult(),
    m_twoViewXOverlap(twoViewXOverlap)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TwoViewTransverseOverlapResult::TwoViewTransverseOverlapResult(const TwoViewTransverseOverlapResult &rhs) :
    TrackTwoViewOverlapResult(rhs),
    m_xOverlap(rhs.IsInitialized() ? rhs.GetTwoViewXOverlap() : TwoViewXOverlap(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f))
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
        m_twoViewXOverlap = rhs.GetTwoViewXOverlap();
    }
    else
    {
        m_twoViewXOverlap = TwoViewXOverlap(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f);
    }

    return *this;
}

} // namespace lar_content
