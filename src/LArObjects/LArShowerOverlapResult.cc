/**
 *  @file   LArContent/src/LArObjects/LArShowerOverlapResult.cc
 * 
 *  @brief  Implementation of the lar shower overlap result class.
 * 
 *  $Log: $
 */

#include "Objects/CaloHit.h"

#include "Pandora/PandoraInternal.h"

#include "LArObjects/LArShowerOverlapResult.h"

using namespace pandora;

namespace lar
{

ShowerOverlapResult::ShowerOverlapResult() :
    m_isInitialized(false),
    m_nMatchedSamplingPoints(0),
    m_nSamplingPoints(0),
    m_matchedFraction(0.f),
    m_xOverlap(XOverlap(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

ShowerOverlapResult::ShowerOverlapResult(const unsigned int nMatchedSamplingPoints, const unsigned int nSamplingPoints, const XOverlap &xOverlap) :
    m_isInitialized(true),
    m_nMatchedSamplingPoints(nMatchedSamplingPoints),
    m_nSamplingPoints(nSamplingPoints),
    m_matchedFraction(0.f),
    m_xOverlap(xOverlap)
{
    if ((0 == m_nSamplingPoints) || (m_nMatchedSamplingPoints > m_nSamplingPoints))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    m_matchedFraction = static_cast<float>(m_nMatchedSamplingPoints) / static_cast<float>(m_nSamplingPoints);
}

//------------------------------------------------------------------------------------------------------------------------------------------

ShowerOverlapResult::ShowerOverlapResult(const ShowerOverlapResult &rhs) :
    m_isInitialized(rhs.m_isInitialized),
    m_nMatchedSamplingPoints(rhs.m_nMatchedSamplingPoints),
    m_nSamplingPoints(rhs.m_nSamplingPoints),
    m_matchedFraction(rhs.m_matchedFraction),
    m_xOverlap(rhs.IsInitialized() ? rhs.GetXOverlap() : XOverlap(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

ShowerOverlapResult::~ShowerOverlapResult()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ShowerOverlapResult::operator<(const ShowerOverlapResult &rhs) const
{
    if (this == &rhs)
        return false;

    if (!m_isInitialized && !rhs.IsInitialized())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    if (!m_isInitialized)
        return true;

    if (!rhs.IsInitialized())
        return false;

    return (m_nMatchedSamplingPoints < rhs.m_nMatchedSamplingPoints);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ShowerOverlapResult::operator>(const ShowerOverlapResult &rhs) const
{
    if (this == &rhs)
        return false;

    return !(*this < rhs);
}

//------------------------------------------------------------------------------------------------------------------------------------------

ShowerOverlapResult &ShowerOverlapResult::operator=(const ShowerOverlapResult &rhs)
{
    if (this == &rhs)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    m_isInitialized = rhs.m_isInitialized;
    m_nMatchedSamplingPoints = rhs.m_nMatchedSamplingPoints;
    m_nSamplingPoints = rhs.m_nSamplingPoints;
    m_matchedFraction = rhs.m_matchedFraction;

    if (rhs.IsInitialized())
    {
        m_xOverlap = rhs.GetXOverlap();
    }
    else
    {
        m_xOverlap = XOverlap(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f);
    }

    return *this;
}

} // namespace lar
