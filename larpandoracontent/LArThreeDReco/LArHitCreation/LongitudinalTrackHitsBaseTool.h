/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/LongitudinalTrackHitsBaseTool.h
 *
 *  @brief  Header file for the longitudinal track hits base tool.
 *
 *  $Log: $
 */
#ifndef LONGITUDINAL_TRACK_HITS_BASE_TOOL_H
#define LONGITUDINAL_TRACK_HITS_BASE_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArHitCreation/TrackHitsBaseTool.h"

namespace lar_content
{

/**
 *  @brief  LongitudinalTrackHitsBaseTool class
 */
class LongitudinalTrackHitsBaseTool : public TrackHitsBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    LongitudinalTrackHitsBaseTool();

protected:
    /**
     *  @brief  Get the three dimensional position using a provided two dimensional calo hit and sliding linear fits in the other two views
     *
     *  @param  matchedSlidingFitMap map of sliding fit results from each view
     *  @param  vtx3D the 3D vertex position
     *  @param  end3D the 3D end position
     *  @param  protoHit to receive the populated proto hit
     */
    virtual void GetLongitudinalTrackHit3D(const MatchedSlidingFitMap &matchedSlidingFitMap, const pandora::CartesianVector &vtx3D,
        const pandora::CartesianVector &end3D, ProtoHit &protoHit) const = 0;

    virtual void GetTrackHits3D(const pandora::CaloHitVector &inputTwoDHits, const MatchedSlidingFitMap &matchedSlidingFitMap,
        ProtoHitVector &protoHitVector) const;

    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

private:
    /**
     *  @brief Get reconstructed vertex and end positions for this 3D track
     *
     *  @param inputSlidingFitMap input map of sliding fit results from each view
     *  @param outputSlidingFitMap output map of clean sliding fit results from each view
     *  @param outputVtx3D reconstructed start position of 3D track
     *  @param outputEnd3D reconstructed end position of 3D track
     */
    void GetVertexAndEndPositions(const MatchedSlidingFitMap &inputSlidingFitMap, MatchedSlidingFitMap &outputSlidingFitMap,
        pandora::CartesianVector &outputVtx3D, pandora::CartesianVector &outputEnd3D) const;

    /**
     *  @brief Combine two 2D coordinates to give a 3D coordinate
     *
     *  @param hitType1 the view corresponding to the first position
     *  @param hitType2 the view corresponding to the second position
     *  @param vtx1 the first position
     *  @param vtx2 the second position
     *  @param bestVtx the combined vertex position
     *  @param bestChi2 the chi-squared from the combination
     */
    void UpdateBestPosition(const pandora::HitType hitType1, const pandora::HitType hitType2, const pandora::CartesianVector &vtx1,
        const pandora::CartesianVector &vtx2, pandora::CartesianVector &bestVtx, float &bestChi2) const;

    float  m_vtxDisplacementCutSquared;   ///<
    float  m_minTrackLengthSquared;       ///<
};

} // namespace lar_content

#endif // #ifndef LONGITUDINAL_TRACK_HITS_BASE_TOOL_H
