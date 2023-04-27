/**
 *  @file   larpandoracontent/LArShowerRefinement/ShowerSpineFinderTool.h
 *
 *  @brief  Header file for the peak direction finder tool class.
  *
 *  $Log: $
 */
#ifndef LAR_SHOWER_SPINE_FINDER_TOOL_H
#define LAR_SHOWER_SPINE_FINDER_TOOL_H 1

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

namespace lar_content
{

class ShowerSpineFinderTool : public pandora::AlgorithmTool
{
public:
    ShowerSpineFinderTool();

    pandora::StatusCode Run(const pandora::CartesianVector &nuVertex3D, const pandora::CaloHitList *const pViewHitList, 
        const pandora::HitType hitType, const pandora::CartesianVector &peakDirection, 
        pandora::CaloHitList &unavailableHitList, pandora::CaloHitList &showerSpineHitList);

private:
    void FindShowerSpine(const pandora::CaloHitList *const pViewShowerHitList, const pandora::CartesianVector &nuVertex2D, 
        const pandora::CartesianVector &initialDirection, pandora::CaloHitList &unavailableHitList, pandora::CaloHitList &showerSpineHitList) const;

    bool CollectSubsectionHits(const TwoDSlidingFitResult &extrapolatedFit, const pandora::CartesianVector &extrapolatedStartPosition, 
        const pandora::CartesianVector &extrapolatedEndPosition, const pandora::CartesianVector &extrapolatedDirection, const bool isEndDownstream, 
        const pandora::CaloHitList *const pViewHitList, pandora::CartesianPointVector &runningFitPositionVector, pandora::CaloHitList &unavailableHitList, 
        pandora::CaloHitList &showerSpineHitList) const;

    bool IsCloseToLine(const pandora::CartesianVector &hitPosition, const pandora::CartesianVector &lineStart, 
        const pandora::CartesianVector &lineDirection, const float distanceToLine) const;

    void CollectConnectedHits(const pandora::CaloHitList &collectedHits, const pandora::CartesianVector &extrapolatedStartPosition, 
        const pandora::CartesianVector &extrapolatedDirection, pandora::CartesianPointVector &runningFitPositionVector, 
        pandora::CaloHitList &showerSpineHitList) const;

    float GetClosestDistance(const pandora::CartesianVector &position, const pandora::CartesianPointVector &testPositions) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_growingFitInitialLength;
    float m_initialFitDistanceToLine;
    unsigned int m_minInitialHitsFound;
    unsigned int m_maxFittingHits;
    unsigned int m_localSlidingFitWindow;
    float m_growingFitSegmentLength;
    unsigned int m_highResolutionSlidingFitWindow;
    float m_distanceToLine;
    float m_hitConnectionDistance;
};

//------------------------------------------------------------------------------------------------------------------------------------------
} // namespace lar_content

#endif // #ifndef LAR_SHOWER_SPINE_FINDER_TOOL_H
