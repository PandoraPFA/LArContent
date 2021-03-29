/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/TrackHitsBaseTool.h
 *
 *  @brief  Header file for the track hits base tool.
 *
 *  $Log: $
 */
#ifndef TRACK_HITS_BASE_TOOL_H
#define TRACK_HITS_BASE_TOOL_H 1

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArThreeDReco/LArHitCreation/HitCreationBaseTool.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  TrackHitsBaseTool class
 */
class TrackHitsBaseTool : public HitCreationBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    TrackHitsBaseTool();

    virtual void Run(ThreeDHitCreationAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pPfo,
        const pandora::CaloHitVector &inputTwoDHits, ProtoHitVector &protoHitVector);

protected:
    typedef std::map<pandora::HitType, TwoDSlidingFitResult> MatchedSlidingFitMap;

    /**
     *  @brief  Calculate 3D hits from an input list of 2D hits
     *
     *  @param  pAlgorithm the hit creation algorithm
     *  @param  inputTwoDHits the input vector of 2D hits
     *  @param  matchedSlidingFitMap the group of sliding fit results
     *  @param  protoHitVector to receive the new three dimensional proto hits
     */
    virtual void GetTrackHits3D(const pandora::CaloHitVector &inputTwoDHits, const MatchedSlidingFitMap &matchedSlidingFitMap,
        ProtoHitVector &protoHitVector) const = 0;

    /**
     *  @brief  Calculate sliding fit results for clusters from each view
     *
     *  @param  pPfo  the input particle flow object
     *  @param  matchedSlidingFitMap  the group of sliding fit results
     */
    virtual void BuildSlidingFitMap(const pandora::ParticleFlowObject *const pPfo, MatchedSlidingFitMap &matchedSlidingFitMap) const;

    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int m_minViews;         ///< The minimum number of views required for building hits
    unsigned int m_slidingFitWindow; ///< The layer window for the sliding linear fits
};

} // namespace lar_content

#endif // #ifndef TRACK_HITS_BASE_TOOL_H
