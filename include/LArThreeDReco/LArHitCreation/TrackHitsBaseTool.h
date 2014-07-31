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

#include "LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.h"

namespace lar
{

/**
 *  @brief  TrackHitsBaseTool class
 */
class TrackHitsBaseTool : public HitCreationTool
{
public:
    void Run(ThreeDHitCreationAlgorithm *pAlgorithm, const pandora::ParticleFlowObject *const pPfo, const pandora::CaloHitList &inputTwoDHits,
        pandora::CaloHitList &newThreeDHits);

protected:
    typedef std::map<pandora::HitType, TwoDSlidingFitResult> MatchedSlidingFitMap;

    /**
     *  @brief  Calculate sliding fit results for clusters from each view
     *
     *  @param  pPfo  the input particle flow object
     *  @param  matchedSlidingFitMap  the group of sliding fit results
     */
    void BuildSlidingFitMap(const pandora::ParticleFlowObject *const pPfo, MatchedSlidingFitMap &matchedSlidingFitMap) const;

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

    /**
     *  @brief  Get the three dimensional position using a provided two dimensional calo hit and candidate fit positions from the other two views
     *
     *  @param  pCaloHit2D address of the two dimensional calo hit
     *  @param  hitType1 the hit type identifying the first view
     *  @param  hitType2 the hit type identifying the second view
     *  @param  fitPosition1 the candidate sliding fit position in the first view
     *  @param  fitPosition2 the candidate sliding fit position in the second view
     *  @param  position3D to receive the three dimensional position
     *  @param  chiSquared to receive the chi squared value
     */
    void GetPosition3D(const pandora::CaloHit *const pCaloHit2D, const pandora::HitType hitType1, const pandora::HitType hitType2,
        const pandora::CartesianVector &fitPosition1, const pandora::CartesianVector &fitPosition2, pandora::CartesianVector &position3D, float &chiSquared) const;

    /**
     *  @brief  Get the three dimensional position using a provided two dimensional calo hit and candidate fit positions from the other two views
     *
     *  @param  pCaloHit2D address of the two dimensional calo hit
     *  @param  hitType the hit type identifying the other view
     *  @param  fitPosition the candidate sliding fit position in the other view
     *  @param  position3D to receive the three dimensional position
     *  @param  chiSquared to receive the chi squared value
     */
    void GetPosition3D(const pandora::CaloHit *const pCaloHit2D, const pandora::HitType hitType, const pandora::CartesianVector &fitPosition,
        pandora::CartesianVector &position3D, float &chiSquared) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int    m_slidingFitWindow;         ///< The layer window for the sliding linear fits
    bool            m_useChiSquaredApproach;    ///< Whether to obtain y, z positions via chi2 approach, or projected position approach
    float           m_chiSquaredCut;            ///< The chi squared cut (accept only values below the cut value)
};

} // namespace lar

#endif // #ifndef TRACK_HITS_BASE_TOOL_H
