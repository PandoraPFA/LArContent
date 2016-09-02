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

class ThreeDHitCreationAlgorithm;

//------------------------------------------------------------------------------------------------------------------------------------------

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

    virtual void Run(ThreeDHitCreationAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pPfo, const pandora::CaloHitList &inputTwoDHits,
        pandora::CaloHitList &newThreeDHits);

protected:
    typedef std::unordered_map<pandora::HitType, TwoDSlidingFitResult, std::hash<unsigned int> > MatchedSlidingFitMap;

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
    virtual void CreateThreeDHits(ThreeDHitCreationAlgorithm *const pAlgorithm, const pandora::CaloHitList &inputTwoDHits,
        const MatchedSlidingFitMap &matchedSlidingFitMap, pandora::CaloHitList &newThreeDHits) const = 0;

    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int    m_minViews;                 ///< The minimum number of views required for building hits
    unsigned int    m_slidingFitWindow;         ///< The layer window for the sliding linear fits
    float           m_chiSquaredCut;            ///< The chi squared cut (accept only values below the cut value)
};

} // namespace lar_content

#endif // #ifndef TRACK_HITS_BASE_TOOL_H
