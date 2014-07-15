/**
 *  @file   LArContent/include/LArObjects/LArXOverlap.h
 * 
 *  @brief  Header file for the lar x overlap class.
 * 
 *  $Log: $
 */
#ifndef LAR_X_OVERLAP_H
#define LAR_X_OVERLAP_H 1

namespace lar
{

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
 *  @brief  x overlap result + operator
 * 
 *  @param  lhs the first x overlap result to add
 *  @param  rhs the second x overlap result to add
 */
XOverlap operator+(const XOverlap &lhs, const XOverlap &rhs);

//------------------------------------------------------------------------------------------------------------------------------------------

inline XOverlap::XOverlap(const float uMinX, const float uMaxX, const float vMinX, const float vMaxX,
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

inline float XOverlap::GetUMinX() const
{
    return m_uMinX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float XOverlap::GetUMaxX() const
{
    return m_uMaxX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float XOverlap::GetVMinX() const
{
    return m_vMinX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float XOverlap::GetVMaxX() const
{
    return m_vMaxX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float XOverlap::GetWMinX() const
{
    return m_wMinX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float XOverlap::GetWMaxX() const
{
    return m_wMaxX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float XOverlap::GetXSpanU() const
{
    return std::fabs(m_uMaxX - m_uMinX);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float XOverlap::GetXSpanV() const
{
    return std::fabs(m_vMaxX - m_vMinX);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float XOverlap::GetXSpanW() const
{
    return std::fabs(m_wMaxX - m_wMinX);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float XOverlap::GetXOverlapSpan() const
{
    return m_xOverlapSpan;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline XOverlap operator+(const XOverlap &lhs, const XOverlap &rhs)
{
    const float uMinX(std::min(lhs.GetUMinX(), rhs.GetUMinX()));
    const float uMaxX(std::max(lhs.GetUMaxX(), rhs.GetUMaxX()));
    const float vMinX(std::min(lhs.GetVMinX(), rhs.GetVMinX()));
    const float vMaxX(std::max(lhs.GetVMaxX(), rhs.GetVMaxX()));
    const float wMinX(std::min(lhs.GetWMinX(), rhs.GetWMinX()));
    const float wMaxX(std::max(lhs.GetWMaxX(), rhs.GetWMaxX()));
    const float minX(std::max(uMinX, std::max(vMinX, wMinX)));
    const float maxX(std::min(uMaxX, std::min(vMaxX, wMaxX)));
    const float xOverlapSpan(maxX - minX);

    return XOverlap(uMinX, uMaxX, vMinX, vMaxX, wMinX, wMaxX, xOverlapSpan);
}

} // namespace lar

#endif // #ifndef LAR_X_OVERLAP_H
