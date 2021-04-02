/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/DeltaRayRemovalTool.h
 *
 *  @brief  Header file for the delta ray removal tool class.
 *
 *  $Log: $
 */
#ifndef DELTA_RAY_REMOVAL_TOOL_H
#define DELTA_RAY_REMOVAL_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/RemovalBaseTool.h"
#include "larpandoracontent/LArThreeDReco/LArCosmicRay/ThreeViewDeltaRayMatchingAlgorithm.h"

namespace lar_content
{
/**
 *  @brief  DeltaRayRemovalTool class
 */
class DeltaRayRemovalTool : public RemovalBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    DeltaRayRemovalTool();

private:
    bool Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor);
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);    

    /**
     *  @brief  Remove hits from cosmic ray clusters that belong to a child delta ray 
     *
     *  @param  pAlgorithm address of the calling algorithm 
     *  @param  elementList the tensor element list
     *  @param  changesMade receive boolean indicating whether clusters in element have been modified
     */
    void ExamineConnectedElements(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::ElementList &elementList,
        bool &changesMade) const;

    /**
     *  @brief  Determine whether element satifies simple checks
     *
     *  @param  pAlgorithm address of the calling algorithm 
     *  @param  element the tensor element
     *  @param  hitType the hit type of the cluster under investigation
     *
     *  @return  whether the checks pass
     */
    virtual bool PassElementChecks(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::Element &element, const pandora::HitType &hitType) const;

    /**
     *  @brief  Determine whether the cosmic ray cluster under investigation has delta ray contamination
     *
     *  @param  pAlgorithm address of the calling algorithm 
     *  @param  element the tensor element
     *  @param  hitType the hit type of the view under investigation
     *
     *  @return  whether the cosmic ray cluster is likely to contain delta ray hits that can be removed
     */
    bool IsContaminated(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::Element &element, const pandora::HitType &hitType) const;

    /**
     *  @brief  Remove collected delta ray hits from the cosmic ray pfo
     *
     *  @param  pAlgorithm address of the calling algorithm 
     *  @param  element the tensor element
     *  @param  hitType the hit type of the cluster under investigation
     *  @param  deltaRayHits the list of delta ray hits to remove
     */
    void SplitMuonCluster(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::Element &element,
        const pandora::HitType &hitType, const pandora::CaloHitList &deltaRayHits) const;

    float m_minSeparation;                    ///< The minimum delta ray - parent cosmic ray separation required to investigate cosmic ray cluster
    unsigned int m_slidingFitWindow;          ///< The sliding fit window used in cosmic ray parameterisations
    float m_minDeviationFromTransverse;       ///< The minimum deviation from transverse required to avoid mistakes
    float m_contaminationWindow;              ///< The distance in which to search for delta ray contamination in the cosmic ray track
    unsigned int m_significantHitThreshold;   ///< The threshold number of hits which define significant contimination
    float m_minDistanceFromMuon; ///< The minimum distance of a hit from the cosmic ray track required for removal
    float m_maxDistanceToCollected; ///< The maximim distance of a hit from the projected delta ray hits required for removal
};

} // namespace lar_content

#endif // #ifndef DELTA_RAY_REMOVAL_TOOL_H