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
    TrackShowerCountingContextObject(const std::vector<float> &nuScores, const std::vector<float> &trackScores, const std::vector<float> &showerScores);

    /**
     *  @brief  Get the neutrino class scores
     *
     *  @return the vector of neutrino class scores
     */
    const std::vector<float> &GetNuClassificationScores() const;

    /**
     *  @brief  Get the track counting scores
     *
     *  @return the vector of track counting scores
     */
    const std::vector<float> &GetTrackCountingScores() const;

    /**
     *  @brief  Get the shower counting scores
     *
     *  @return the vector of shower counting scores
     */
    const std::vector<float> &GetShowerCountingScores() const;

private:

    std::vector<float> m_nuScores;      ///< Vector of neutrino class scores 
    std::vector<float> m_trackScores;   ///< Vector of track counting scores
    std::vector<float> m_showerScores;  ///< Vector of shower counting scores
};

//-----------------------------------------------------------------------------------------------------------------------------------------

inline TrackShowerCountingContextObject::TrackShowerCountingContextObject(const std::vector<float> &nuScores, const std::vector<float> &trackScores, const std::vector<float> &showerScores) :
    m_nuScores(nuScores),
    m_trackScores(trackScores),
    m_showerScores(showerScores)
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

inline const std::vector<float> &TrackShowerCountingContextObject::GetNuClassificationScores() const
{
    return m_nuScores;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

inline const std::vector<float> &TrackShowerCountingContextObject::GetTrackCountingScores() const
{
    return m_trackScores;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

inline const std::vector<float> &TrackShowerCountingContextObject::GetShowerCountingScores() const
{
    return m_showerScores;
}

}

#endif 
