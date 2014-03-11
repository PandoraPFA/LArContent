/**
 *  @file   LArContent/src/LArObjects/LArTrackOverlapResult.cc
 * 
 *  @brief  Implementation of the lar track overlap result class.
 * 
 *  $Log: $
 */

#include "Pandora/PandoraInternal.h"
#include "Pandora/StatusCodes.h"

#include "LArObjects/LArTrackOverlapResult.h"

using namespace pandora;

namespace lar
{

TrackOverlapResult::TrackOverlapResult() :
    m_isInitialized(false),
    m_nMatchedSamplingPoints(0),
    m_nSamplingPoints(0),
    m_matchedFraction(0.f),
    m_chi2(0.f),
    m_reducedChi2(0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackOverlapResult::TrackOverlapResult(const unsigned int nMatchedSamplingPoints, const unsigned int nSamplingPoints, const float chi2) :
    m_isInitialized(true),
    m_nMatchedSamplingPoints(nMatchedSamplingPoints),
    m_nSamplingPoints(nSamplingPoints),
    m_matchedFraction(0.f),
    m_chi2(chi2),
    m_reducedChi2(0.f)
{
    if ((0 == m_nSamplingPoints) || (m_nMatchedSamplingPoints > m_nSamplingPoints))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    m_matchedFraction = static_cast<float>(m_nMatchedSamplingPoints) / static_cast<float>(m_nSamplingPoints);
    m_reducedChi2 = m_chi2 / static_cast<float>(m_nSamplingPoints);
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackOverlapResult::TrackOverlapResult(const TrackOverlapResult &rhs) :
    m_isInitialized(rhs.m_isInitialized),
    m_nMatchedSamplingPoints(rhs.m_nMatchedSamplingPoints),
    m_nSamplingPoints(rhs.m_nSamplingPoints),
    m_matchedFraction(rhs.m_matchedFraction),
    m_chi2(rhs.m_chi2),
    m_reducedChi2(rhs.m_reducedChi2)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackOverlapResult::operator<(const TrackOverlapResult &rhs) const
{
    if (m_nMatchedSamplingPoints != rhs.m_nMatchedSamplingPoints)
        return (m_nMatchedSamplingPoints < rhs.m_nMatchedSamplingPoints);

    return (m_reducedChi2 > rhs.m_reducedChi2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackOverlapResult::operator>(const TrackOverlapResult &rhs) const
{
    if (this == &rhs)
        return false;

    return !(*this < rhs);
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackOverlapResult &TrackOverlapResult::operator=(const TrackOverlapResult &rhs)
{
    if (this == &rhs)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    m_isInitialized = rhs.m_isInitialized;
    m_nMatchedSamplingPoints = rhs.m_nMatchedSamplingPoints;
    m_nSamplingPoints = rhs.m_nSamplingPoints;
    m_matchedFraction = rhs.m_matchedFraction;
    m_chi2 = rhs.m_chi2;
    m_reducedChi2 = rhs.m_reducedChi2;

    return *this;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

TrackOverlapResult operator+(const TrackOverlapResult &lhs, const TrackOverlapResult &rhs)
{
    if (!lhs.IsInitialized() && !rhs.IsInitialized())
        return TrackOverlapResult();

    return TrackOverlapResult(lhs.GetNMatchedSamplingPoints() + rhs.GetNMatchedSamplingPoints(),
        lhs.GetNSamplingPoints() + rhs.GetNSamplingPoints(), lhs.GetChi2() + rhs.GetChi2());
}

} // namespace lar
