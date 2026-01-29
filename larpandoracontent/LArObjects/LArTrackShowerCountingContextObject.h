/**
 *  @file   larpandoracontent/LArObjects/LArTrackShowerCountingContextObject.h
 *
 *  @brief  Header file for the event context object for track and shower counting
 *
 *  $Log: $
 */
#ifndef LAR_TRACK_SHOWER_COUNTING_CONTEXT_OBJECT_H
#define LAR_TRACK_SHOWER_COUNTING_CONTEXT_OBJECT_H 1

#include <vector>

#include "Objects/EventContext.h"

namespace lar_content
{

/**
 *  @brief  TrackShowerCountingContextObject class
 */
class TrackShowerCountingContextObject : public pandora::EventContextObject
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  nuScores the vector of neutrino class scores
     *  @param  trackScores the vector of track counting scores
     *  @param  showerScores the vector of shower counting scores
     */
    TrackShowerCountingContextObject(const pandora::FloatVector &nuScores, const pandora::FloatVector &trackScores, const pandora::FloatVector &showerScores);

    /**
     *  @brief  Get the neutrino class scores
     *
     *  @return the vector of neutrino class scores
     */
    const pandora::FloatVector &GetNuClassificationScores() const;

    /**
     *  @brief  Get the track counting scores
     *
     *  @return the vector of track counting scores
     */
    const pandora::FloatVector &GetTrackCountingScores() const;

    /**
     *  @brief  Get the shower counting scores
     *
     *  @return the vector of shower counting scores
     */
    const pandora::FloatVector &GetShowerCountingScores() const;

private:
    pandora::FloatVector m_nuScores;     ///< Vector of neutrino class scores
    pandora::FloatVector m_trackScores;  ///< Vector of track counting scores
    pandora::FloatVector m_showerScores; ///< Vector of shower counting scores
};

//-----------------------------------------------------------------------------------------------------------------------------------------

inline TrackShowerCountingContextObject::TrackShowerCountingContextObject(
    const pandora::FloatVector &nuScores, const pandora::FloatVector &trackScores, const pandora::FloatVector &showerScores) :
    m_nuScores(nuScores),
    m_trackScores(trackScores),
    m_showerScores(showerScores)
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::FloatVector &TrackShowerCountingContextObject::GetNuClassificationScores() const
{
    return m_nuScores;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::FloatVector &TrackShowerCountingContextObject::GetTrackCountingScores() const
{
    return m_trackScores;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::FloatVector &TrackShowerCountingContextObject::GetShowerCountingScores() const
{
    return m_showerScores;
}

} // namespace lar_content

#endif
