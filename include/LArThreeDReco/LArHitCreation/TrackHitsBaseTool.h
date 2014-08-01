/**
 *  @file   LArContent/include/LArThreeDReco/LArHitCreation/TrackHitsBaseTool.h
 *
 *  @brief  Header file for the track hits base tool.
 *
 *  $Log: $
 */
#ifndef TRACK_HITS_BASE_TOOL_H
#define TRACK_HITS_BASE_TOOL_H 1

#include "LArObjects/LArTwoDSlidingFitResult.h"

#include "LArThreeDReco/LArHitCreation/HitCreationBaseTool.h"

namespace lar
{

class ThreeDHitCreationAlgorithm;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  TrackHitsBaseTool class
 */
class TrackHitsBaseTool : public HitCreationBaseTool
{
public:
    virtual void Run(ThreeDHitCreationAlgorithm *pAlgorithm, const pandora::ParticleFlowObject *const pPfo, const pandora::CaloHitList &inputTwoDHits,
        pandora::CaloHitList &newThreeDHits);

protected:
    typedef std::map<pandora::HitType, TwoDSlidingFitResult> MatchedSlidingFitMap;

    /**
     *  @brief  Calculate sliding fit results for clusters from each view
     *
     *  @param  pPfo  the input particle flow object
     *  @param  matchedSlidingFitMap  the group of sliding fit results
     */
    virtual void BuildSlidingFitMap(const pandora::ParticleFlowObject *const pPfo, MatchedSlidingFitMap &matchedSlidingFitMap) const;

    /**
     *  @brief  Calculate 3D hits from an input list of 2D hits
     *
     *  @param  pAlgorithm  the hit creation algorithm
     *  @param  inputTwoDHits  the input list of 2D hits
     *  @param  matchedSlidingFitMap  the group of sliding fit results
     *  @param  newThreeDHits  the output list of 3D hits
     *  @param  omittedTwoDHits  the output list of omitted 2D hits
     */
    virtual void CreateThreeDHits(ThreeDHitCreationAlgorithm *pAlgorithm, const pandora::CaloHitList &inputTwoDHits,
        const MatchedSlidingFitMap &matchedSlidingFitMap, pandora::CaloHitList &newThreeDHits) const = 0;

    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int    m_minViews;                 ///< The minimum number of views required for building hits
    unsigned int    m_slidingFitWindow;         ///< The layer window for the sliding linear fits
    float           m_chiSquaredCut;            ///< The chi squared cut (accept only values below the cut value)
};

} // namespace lar

#endif // #ifndef TRACK_HITS_BASE_TOOL_H
