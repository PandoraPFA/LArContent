/**
 *  @file   larpandoracontent/LArObjects/LArTwoViewXOverlap.h
 *
 *  @brief  Header file for the lar x two view overlap class.
 *
 *  $Log: $
 */
#ifndef LAR_TWO_VIEW_X_OVERLAP_H
#define LAR_TWO_VIEW_X_OVERLAP_H 1

#include <algorithm>
#include <cmath>
#include <limits>

namespace lar_content
{

/**
 *  @brief  TwoViewXOverlap class
 */
class TwoViewXOverlap
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  xMin0 min x value in the view 0
     *  @param  xMax0 max x value in the view 0
     *  @param  xMin1 min x value in the view 1
     *  @param  xMax1 max x value in the view 1
     */
    TwoViewXOverlap(const float xMin0, const float xMax0, const float xMin1, const float xMax1);

    /**
     *  @brief  Get the min x value in the view 0
     *
     *  @return the min x value in the view 0
     */
    float GetXMin0() const;

    /**
     *  @brief  Get the max x value in the view 0
     *
     *  @return the max x value in the view 0
     */
    float GetXMax0() const;

    /**
     *  @brief  Get the min x value in the view 1
     *
     *  @return the min x value in the view 1
     */
    float GetXMin1() const;

    /**
     *  @brief  Get the max x value in the view 1
     *
     *  @return the max x value in the view 1
     */
    float GetXMax1() const;

    /**
     *  @brief  Get the x span in the view 0
     *
     *  @return the x span in the view 0
     */
    float GetXSpan0() const;

    /**
     *  @brief  Get the x span in the view 1
     *
     *  @return the x span in the view 1
     */
    float GetXSpan1() const;

    /**
     *  @brief  Get the x overlap span
     *
     *  @return the x overlap span
     */
    float GetTwoViewXOverlapSpan() const;

    /**
     *  @brief  Get the x overlap max X value
     *
     *  @return the x overlap min
     */
    float GetTwoViewXOverlapMin() const;

    /**
     *  @brief  Get the x overlap min X value
     *
     *  @return the x overlap max
     */
    float GetTwoViewXOverlapMax() const;

    /**
     *  @brief  Get the fraction of the view 0 cluster that overlaps in x
     *
     *  @return the view 0 cluster's fractional overlap
     */
    float GetXOverlapFraction0() const;

    /**
     *  @brief  Get the fraction of the view 1 cluster that overlaps in x
     *
     *  @return the view 1 cluster's fractional overlap
     */
    float GetXOverlapFraction1() const;

private:
    float m_xMin0;        ///< The min x value in the view 0
    float m_xMax0;        ///< The max x value in the view 0
    float m_xMin1;        ///< The min x value in the view 1
    float m_xMax1;        ///< The max x value in the view 1
    float m_xOverlapSpan; ///< The x overlap span
};

/**
 *  @brief  x overlap result + operator
 *
 *  @param  lhs the first x overlap result to add
 *  @param  rhs the second x overlap result to add
 */
TwoViewXOverlap operator+(const TwoViewXOverlap &lhs, const TwoViewXOverlap &rhs);

//------------------------------------------------------------------------------------------------------------------------------------------

inline TwoViewXOverlap::TwoViewXOverlap(const float xMin0, const float xMax0, const float xMin1, const float xMax1) :
    m_xMin0(xMin0),
    m_xMax0(xMax0),
    m_xMin1(xMin1),
    m_xMax1(xMax1),
    m_xOverlapSpan(std::min(m_xMax0, m_xMax1) - std::max(m_xMin0, m_xMin1))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TwoViewXOverlap::GetXMin0() const
{
    return m_xMin0;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TwoViewXOverlap::GetXMax0() const
{
    return m_xMax0;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TwoViewXOverlap::GetXMin1() const
{
    return m_xMin1;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TwoViewXOverlap::GetXMax1() const
{
    return m_xMax1;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TwoViewXOverlap::GetXSpan0() const
{
    return std::fabs(m_xMax0 - m_xMin0);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TwoViewXOverlap::GetXSpan1() const
{
    return std::fabs(m_xMax1 - m_xMin1);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TwoViewXOverlap::GetTwoViewXOverlapSpan() const
{
    return m_xOverlapSpan;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TwoViewXOverlap::GetTwoViewXOverlapMin() const
{
    return std::max(m_xMin0, m_xMin1);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TwoViewXOverlap::GetTwoViewXOverlapMax() const
{
    return std::min(m_xMax0, m_xMax1);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TwoViewXOverlap::GetXOverlapFraction0() const
{
    return (std::numeric_limits<float>::epsilon() < this->GetXSpan0()) ? m_xOverlapSpan / this->GetXSpan0() : 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TwoViewXOverlap::GetXOverlapFraction1() const
{
    return (std::numeric_limits<float>::epsilon() < this->GetXSpan1()) ? m_xOverlapSpan / this->GetXSpan1() : 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TwoViewXOverlap operator+(const TwoViewXOverlap &lhs, const TwoViewXOverlap &rhs)
{
    const float xMin0(std::min(lhs.GetXMin0(), rhs.GetXMin0()));
    const float xMax0(std::max(lhs.GetXMax0(), rhs.GetXMax0()));
    const float xMin1(std::min(lhs.GetXMin1(), rhs.GetXMin1()));
    const float xMax1(std::max(lhs.GetXMax1(), rhs.GetXMax1()));

    return TwoViewXOverlap(xMin0, xMax0, xMin1, xMax1);
}

} // namespace lar_content

#endif // #ifndef LAR_TWO_VIEW_X_OVERLAP_H
