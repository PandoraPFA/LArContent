/**
 *  @file   LArContent/include/LArObjects/LArShowerOverlapResult.h
 * 
 *  @brief  Header file for the lar shower overlap result class.
 * 
 *  $Log: $
 */
#ifndef LAR_SHOWER_OVERLAP_RESULT_H
#define LAR_SHOWER_OVERLAP_RESULT_H 1

#include "Pandora/PandoraInputTypes.h"
#include "Pandora/StatusCodes.h"

#include "LArObjects/LArXOverlap.h"

#include <cmath>
#include <vector>

namespace lar_content
{

/**
 *  @brief  ShowerOverlapResult class
 */
class ShowerOverlapResult
{
public:
    /**
     *  @brief  Default constructor
     */
    ShowerOverlapResult();

    /**
     *  @brief  Constructor
     * 
     *  @param  nMatchedSamplingPoints
     *  @param  nSamplingPoints
     */
    ShowerOverlapResult(const unsigned int nMatchedSamplingPoints, const unsigned int nSamplingPoints, const XOverlap &xOverlap);

    /**
     *  @brief  Copy constructor
     * 
     *  @param  rhs
     */
    ShowerOverlapResult(const ShowerOverlapResult &rhs);

    /**
     *  @brief  Destructor
     */
    ~ShowerOverlapResult();

    /**
     *  @brief  Whether the track overlap result has been initialized
     *
     *  @return boolean
     */
    bool IsInitialized() const;

    /**
     *  @brief  Get the number of matched sampling points
     *
     *  @return the number of matched sampling points
     */
    unsigned int GetNMatchedSamplingPoints() const;

    /**
     *  @brief  Get the number of sampling points
     *
     *  @return the number of sampling points
     */
    unsigned int GetNSamplingPoints() const;

    /**
     *  @brief  Get the fraction of sampling points resulting in a match
     *
     *  @return the fraction of sampling points resulting in a match
     */
    float GetMatchedFraction() const;

    /**
     *  @brief  Get the x overlap object
     * 
     *  @return the x overlap object
     */
    const XOverlap &GetXOverlap() const;

    /**
     *  @brief  Track overlap result less than operator
     * 
     *  @param  rhs the track overlap result for comparison
     */
    bool operator<(const ShowerOverlapResult &rhs) const;

    /**
     *  @brief  Track overlap result greater than operator
     * 
     *  @param  rhs the track overlap result for comparison
     */
    bool operator>(const ShowerOverlapResult &rhs) const;

    /**
     *  @brief  Track overlap result assigment operator
     * 
     *  @param  rhs the track overlap result to assign
     */
    ShowerOverlapResult &operator=(const ShowerOverlapResult &rhs);

protected:
    bool            m_isInitialized;                ///< Whether the track overlap result has been initialized
    unsigned int    m_nMatchedSamplingPoints;       ///< The number of matched sampling points
    unsigned int    m_nSamplingPoints;              ///< The number of sampling points
    float           m_matchedFraction;              ///< The fraction of sampling points resulting in a match
    XOverlap        m_xOverlap;                     ///< The x overlap object
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool ShowerOverlapResult::IsInitialized() const
{
    return m_isInitialized;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int ShowerOverlapResult::GetNMatchedSamplingPoints() const
{
    if (m_isInitialized)
        return m_nMatchedSamplingPoints;

    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int ShowerOverlapResult::GetNSamplingPoints() const
{
    if (m_isInitialized)
        return m_nSamplingPoints;

    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ShowerOverlapResult::GetMatchedFraction() const
{
    if (m_isInitialized)
        return m_matchedFraction;

    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const XOverlap &ShowerOverlapResult::GetXOverlap() const
{
    if (m_isInitialized)
        return m_xOverlap;

    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);
}

} // namespace lar_content

#endif // #ifndef LAR_SHOWER_OVERLAP_RESULT_H
