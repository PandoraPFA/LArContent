/**
 *  @file   larpandoracontent/LArObjects/LArTwoViewXOverlap.h
 *
 *  @brief  Header file for the lar x overlap class.
 *
 *  $Log: $
 */
#ifndef LAR_TWO_VIEW_X_OVERLAP_H
#define LAR_TWO_VIEW_X_OVERLAP_H 1

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
     *  @param  uMinX min x value in the u view
     *  @param  uMaxX max x value in the u view
     *  @param  vMinX min x value in the v view
     *  @param  vMaxX max x value in the v view
     *  @param  xOverlapSpan the x overlap span
     */
    TwoViewXOverlap(const float uMinX, const float uMaxX, const float vMinX, const float vMaxX, const float xOverlapSpan);

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
     *  @brief  Get the x overlap span
     *
     *  @return the x overlap span
     */
    float GetTwoViewXOverlapSpan() const;

    /**
     *  @brief  Get the fraction of the U cluster that overlaps in x
     *
     *  @return the U cluster's fractional overlap
     */
    float GetXOverlapFractionU() const;

    /**
     *  @brief  Get the fraction of the V cluster that overlaps in x
     *
     *  @return the V cluster's fractional overlap
     */
    float GetXOverlapFractionV() const;

private:
    float       m_uMinX;                        ///< The min x value in the u view
    float       m_uMaxX;                        ///< The max x value in the u view
    float       m_vMinX;                        ///< The min x value in the v view
    float       m_vMaxX;                        ///< The max x value in the v view
    float       m_xOverlapSpan;                 ///< The x overlap span
};

/**
 *  @brief  x overlap result + operator
 *
 *  @param  lhs the first x overlap result to add
 *  @param  rhs the second x overlap result to add
 */
TwoViewXOverlap operator+(const TwoViewXOverlap &lhs, const TwoViewXOverlap &rhs);

//------------------------------------------------------------------------------------------------------------------------------------------

inline TwoViewXOverlap::TwoViewXOverlap(const float uMinX, const float uMaxX, const float vMinX, const float vMaxX,
        const float xOverlapSpan) :
    m_uMinX(uMinX),
    m_uMaxX(uMaxX),
    m_vMinX(vMinX),
    m_vMaxX(vMaxX),
    m_xOverlapSpan(xOverlapSpan)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TwoViewXOverlap::GetUMinX() const
{
    return m_uMinX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TwoViewXOverlap::GetUMaxX() const
{
    return m_uMaxX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TwoViewXOverlap::GetVMinX() const
{
    return m_vMinX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TwoViewXOverlap::GetVMaxX() const
{
    return m_vMaxX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TwoViewXOverlap::GetXSpanU() const
{
    return std::fabs(m_uMaxX - m_uMinX);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TwoViewXOverlap::GetXSpanV() const
{
    return std::fabs(m_vMaxX - m_vMinX);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TwoViewXOverlap::GetTwoViewXOverlapSpan() const
{
    return m_xOverlapSpan;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TwoViewXOverlap::GetXOverlapFractionU() const
{
    return (std::numeric_limits<float>::epsilon() < GetXSpanU()) ? m_xOverlapSpan / GetXSpanU() : 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TwoViewXOverlap::GetXOverlapFractionV() const
{
    return (std::numeric_limits<float>::epsilon() < GetXSpanU()) ? m_xOverlapSpan / GetXSpanV() : 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TwoViewXOverlap operator+(const TwoViewXOverlap &lhs, const TwoViewXOverlap &rhs)
{
    const float uMinX(std::min(lhs.GetUMinX(), rhs.GetUMinX()));
    const float uMaxX(std::max(lhs.GetUMaxX(), rhs.GetUMaxX()));
    const float vMinX(std::min(lhs.GetVMinX(), rhs.GetVMinX()));
    const float vMaxX(std::max(lhs.GetVMaxX(), rhs.GetVMaxX()));
    const float minX(std::max(uMinX, vMinX));
    const float maxX(std::min(uMaxX, vMaxX));
    const float xOverlapSpan(maxX - minX);

    return TwoViewXOverlap(uMinX, uMaxX, vMinX, vMaxX, xOverlapSpan);
}


} // namespace lar_content

#endif // #ifndef LAR_TWO_VIEW_X_OVERLAP_H
