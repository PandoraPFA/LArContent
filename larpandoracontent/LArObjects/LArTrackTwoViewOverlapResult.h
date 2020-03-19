/**
 *  @file   larpandoracontent/LArObjects/LArTrackTwoViewOverlapResult.h
 *
 *  @brief  Header file for the lar track two view overlap result class.
 *
 *  $Log: $
 */
#ifndef LAR_TRACK_TWO_VIEW_OVERLAP_RESULT_H
#define LAR_TRACK_TWO_VIEW_OVERLAP_RESULT_H 1

#include "Pandora/PandoraInputTypes.h"
#include "Pandora/StatusCodes.h"

#include "larpandoracontent/LArObjects/LArTwoViewXOverlap.h"

#include <cmath>
#include <vector>

namespace lar_content
{

/**
 *  @brief  TrackTwoViewOverlapResult class
 */
class TrackTwoViewOverlapResult
{
public:
    /**
     *  @brief  Default constructor
     */
    TrackTwoViewOverlapResult();

    /**
     *  @brief  Copy constructor
     *
     *  @param  rhs
     */
    TrackTwoViewOverlapResult(const TrackTwoViewOverlapResult &rhs);

    /**
     *  @brief  Destructor
     */
    virtual ~TrackTwoViewOverlapResult();

    /**
     *  @brief  Whether the track overlap result has been initialized
     *
     *  @return boolean
     */
    bool IsInitialized() const;

    /**
     *  @brief  Track overlap result assigment operator
     *
     *  @param  rhs the track overlap result to assign
     */
    TrackTwoViewOverlapResult &operator=(const TrackTwoViewOverlapResult &rhs);

protected:
    bool            m_isInitialized;                ///< Whether the track overlap result has been initialized
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  TwoViewTransverseOverlapResult class
 */
class TwoViewTransverseOverlapResult : public TrackTwoViewOverlapResult
{
public:
    /**
     *  @brief  Default constructor
     */
    TwoViewTransverseOverlapResult();

    /**
     *  @brief  Constructor
     *
     *  @param  twoViewXOverlap
     */
    TwoViewTransverseOverlapResult(const TwoViewXOverlap &twoViewXOverlap);

    /**
     *  @brief  Copy constructor
     *
     *  @param  rhs
     */
    TwoViewTransverseOverlapResult(const TwoViewTransverseOverlapResult &rhs);

    /**
     *  @brief  Destructor
     */
    ~TwoViewTransverseOverlapResult();

    /**
     *  @brief  Get the two view x overlap object
     *
     *  @return the two view x overlap object
     */
    const TwoViewXOverlap &GetTwoViewXOverlap() const;

    /**
     *  @brief  Track overlap result assigment operator
     *
     *  @param  rhs the track overlap result to assign
     */
    TwoViewTransverseOverlapResult &operator=(const TwoViewTransverseOverlapResult &rhs);

private:
    TwoViewXOverlap        m_twoViewXOverlap;                     ///< The two view  x overlap object
};

typedef std::vector<TwoViewTransverseOverlapResult> TwoViewTransverseOverlapResultVector;

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline bool TrackTwoViewOverlapResult::IsInitialized() const
{
    return m_isInitialized;
}

} // namespace lar_content

#endif // #ifndef LAR_TRACK_OVERLAP_RESULT_H
