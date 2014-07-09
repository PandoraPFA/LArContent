/**
 *  @file   LArContent/src/LArObjects/LArTrackOverlapResult.cc
 * 
 *  @brief  Implementation of the lar track overlap result class.
 * 
 *  $Log: $
 */

#include "Objects/CaloHit.h"

#include "Pandora/PandoraInternal.h"

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

TrackOverlapResult::~TrackOverlapResult()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackOverlapResult::operator<(const TrackOverlapResult &rhs) const
{
    if (this == &rhs)
        return false;

    if (!m_isInitialized && !rhs.IsInitialized())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    if (!m_isInitialized)
        return true;

    if (!rhs.IsInitialized())
        return false;

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

TransverseOverlapResult::TransverseOverlapResult() :
    TrackOverlapResult(),
    m_xOverlap(XOverlap(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TransverseOverlapResult::TransverseOverlapResult(const unsigned int nMatchedSamplingPoints, const unsigned int nSamplingPoints,
        const float chi2, const XOverlap &xOverlap) :
    TrackOverlapResult(nMatchedSamplingPoints, nSamplingPoints, chi2),
    m_xOverlap(xOverlap)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TransverseOverlapResult::TransverseOverlapResult(const TransverseOverlapResult &rhs) :
    TrackOverlapResult(rhs),
    m_xOverlap(rhs.IsInitialized() ? rhs.GetXOverlap() : XOverlap(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TransverseOverlapResult::~TransverseOverlapResult()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TransverseOverlapResult &TransverseOverlapResult::operator=(const TransverseOverlapResult &rhs)
{
    this->TrackOverlapResult::operator=(rhs);
    m_xOverlap = rhs.GetXOverlap();

    return *this;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

TransverseOverlapResult operator+(const TransverseOverlapResult &lhs, const TransverseOverlapResult &rhs)
{
    if (!lhs.IsInitialized() && !rhs.IsInitialized())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    if (!lhs.IsInitialized())
        return rhs;

    if (!rhs.IsInitialized())
        return lhs;

    const TransverseOverlapResult::XOverlap &xOverlapLhs(lhs.GetXOverlap());
    const TransverseOverlapResult::XOverlap &xOverlapRhs(rhs.GetXOverlap());

    const TransverseOverlapResult::XOverlap xOverlapSum(
        std::min(xOverlapLhs.GetUMinX(), xOverlapRhs.GetUMinX()), std::max(xOverlapLhs.GetUMaxX(), xOverlapRhs.GetUMaxX()),
        std::min(xOverlapLhs.GetVMinX(), xOverlapRhs.GetVMinX()), std::max(xOverlapLhs.GetVMaxX(), xOverlapRhs.GetVMaxX()),
        std::min(xOverlapLhs.GetWMinX(), xOverlapRhs.GetWMinX()), std::max(xOverlapLhs.GetWMaxX(), xOverlapRhs.GetWMaxX()),
        xOverlapLhs.GetXOverlapSpan() + xOverlapRhs.GetXOverlapSpan());

    return TransverseOverlapResult(lhs.GetNMatchedSamplingPoints() + rhs.GetNMatchedSamplingPoints(),
        lhs.GetNSamplingPoints() + rhs.GetNSamplingPoints(), lhs.GetChi2() + rhs.GetChi2(), xOverlapSum);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LongitudinalOverlapResult::LongitudinalOverlapResult() :
    TrackOverlapResult(),
    m_innerChi2(0.f),
    m_outerChi2(0.f)
{
}
//------------------------------------------------------------------------------------------------------------------------------------------

LongitudinalOverlapResult::LongitudinalOverlapResult(const TrackOverlapResult trackOverlapResult, const float innerChi2, const float outerChi2) :
    TrackOverlapResult(trackOverlapResult),
    m_innerChi2(innerChi2),
    m_outerChi2(outerChi2)
{

}

//------------------------------------------------------------------------------------------------------------------------------------------

LongitudinalOverlapResult::LongitudinalOverlapResult(const unsigned int nMatchedSamplingPoints, const unsigned int nSamplingPoints,
        const float chi2, const float innerChi2, const float outerChi2) :
    TrackOverlapResult(nMatchedSamplingPoints, nSamplingPoints, chi2),
    m_innerChi2(innerChi2),
    m_outerChi2(outerChi2)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

LongitudinalOverlapResult::LongitudinalOverlapResult(const LongitudinalOverlapResult &rhs) :
    TrackOverlapResult(rhs),
    m_innerChi2(rhs.IsInitialized() ? rhs.GetInnerChi2() : 0.f),
    m_outerChi2(rhs.IsInitialized() ? rhs.GetOuterChi2() : 0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

LongitudinalOverlapResult::~LongitudinalOverlapResult()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

LongitudinalOverlapResult &LongitudinalOverlapResult::operator=(const LongitudinalOverlapResult &rhs)
{
    this->TrackOverlapResult::operator=(rhs);
    m_innerChi2 = rhs.GetInnerChi2();
    m_outerChi2 = rhs.GetOuterChi2();

    return *this;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

FragmentOverlapResult::FragmentOverlapResult() :
    TrackOverlapResult(),
    m_caloHitList(),
    m_clusterList()
{
}
//------------------------------------------------------------------------------------------------------------------------------------------

FragmentOverlapResult::FragmentOverlapResult(const TrackOverlapResult trackOverlapResult, const pandora::CaloHitList &caloHitList,
        const pandora::ClusterList &clusterList) :
    TrackOverlapResult(trackOverlapResult),
    m_caloHitList(caloHitList),
    m_clusterList(clusterList)
{

}

//------------------------------------------------------------------------------------------------------------------------------------------

FragmentOverlapResult::FragmentOverlapResult(const unsigned int nMatchedSamplingPoints, const unsigned int nSamplingPoints,
        const float chi2, const pandora::CaloHitList &caloHitList, const pandora::ClusterList &clusterList) :
    TrackOverlapResult(nMatchedSamplingPoints, nSamplingPoints, chi2),
    m_caloHitList(caloHitList),
    m_clusterList(clusterList)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

FragmentOverlapResult::FragmentOverlapResult(const FragmentOverlapResult &rhs) :
    TrackOverlapResult(rhs),
    m_caloHitList(rhs.IsInitialized() ? rhs.GetFragmentCaloHitList() : CaloHitList()),
    m_clusterList(rhs.IsInitialized() ? rhs.GetFragmentClusterList() : ClusterList())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

FragmentOverlapResult::~FragmentOverlapResult()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

FragmentOverlapResult &FragmentOverlapResult::operator=(const FragmentOverlapResult &rhs)
{
    this->TrackOverlapResult::operator=(rhs);
    m_caloHitList = rhs.GetFragmentCaloHitList();
    m_clusterList = rhs.GetFragmentClusterList();

    return *this;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::HitType FragmentOverlapResult::GetFragmentHitType() const
{
    if (m_caloHitList.empty())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    return (*(m_caloHitList.begin()))->GetHitType();
}

} // namespace lar
