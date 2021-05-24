/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/RemovalBaseTool.h
 *
 *  @brief  Header file for the removal base tool class.
 *
 *  $Log: $
 */
#ifndef REMOVAL_BASE_TOOL_H
#define REMOVAL_BASE_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/ThreeViewDeltaRayMatchingAlgorithm.h"

namespace lar_content
{
/**
 *  @brief  RemovalBaseTool class
 */
class RemovalBaseTool : public DeltaRayTensorTool
{
public:
    typedef std::vector<pandora::HitType> HitTypeVector;

    /**
     *  @brief  Default constructor
     */
    RemovalBaseTool();

protected:
    virtual bool Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor) = 0;
    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle) = 0;

    /**
     *  @brief  Determine whether the matched clusters suggest that the delta ray is at the endpoint of the cosmic ray (and is likely to be a michel)
     *
     *  @param  element the tensor element
     *  @param  ignoreHitType whether to ignore the cluster under investigation
     *  @param  hitTypeToIgnore the hit type to ignore
     *
     *  @return  whether the delta ray is at the endpoint of the cosmic ray
     */
    bool IsMuonEndpoint(const TensorType::Element &element, const bool ignoreHitType, const pandora::HitType hitTypeToIgnore = pandora::TPC_VIEW_U) const;

    /**
     *  @brief  Determine whether the input element is the best to use to modify the contaminated cluster (best is defined by the total hit count)
     *
     *  @param  element the tensor element
     *  @param  hitType the hit type of the cluster under investigation
     *  @param  elementList the tensor element list
     *
     *  @return  whether the input element is the best element
     */
    bool IsBestElement(const TensorType::Element &element, const pandora::HitType hitType, const TensorType::ElementList &elementList) const;

    /**
     *  @brief  Determine whether element satifies simple checks
     *
     *  @param  element the tensor element
     *  @param  hitType the hit type of the cluster under investigation
     *
     *  @return  whether the checks pass
     */
    virtual bool PassElementChecks(const TensorType::Element &element, const pandora::HitType hitType) const = 0;

    /**
     *  @brief  Whether a given position is close to a defined line 
     *
     *  @param  hitPosition the input position
     *  @param  lineStart the start position of the line
     *  @param  lineEnd the end position of the line
     *  @param  distanceToLine the definition of close
     *
     *  @return  whether the position is close to the definied line
     */
    bool IsCloseToLine(const pandora::CartesianVector &hitPosition, const pandora::CartesianVector &lineStart,
        const pandora::CartesianVector &lineEnd, const float distanceToLine) const;

    /**
     *  @brief  Whether the projection of a given position lies on a defined line
     *
     *  @param  lowerBoundary the start position of the line
     *  @param  upperBoundary the end position of the line
     *  @param  point the input position
     *
     *  @return  whether the position lies between the two points
     */
    bool IsInLineSegment(const pandora::CartesianVector &lowerBoundary, const pandora::CartesianVector &upperBoundary,
        const pandora::CartesianVector &point) const;

    /**
     *  @brief  Collect the hits that are closest to and can be projected onto a defined line 
     *
     *  @param  pCluster the address of the input cluster containing the hits to be investigated
     *  @param  lowerBoundary the start position of the line
     *  @param  upperBoundary the end position of the line
     *  @param  collectedHits the collected hits
     */
    void FindExtrapolatedHits(const pandora::Cluster *const pCluster, const pandora::CartesianVector &lowerBoundary,
        const pandora::CartesianVector &upperBoundary, pandora::CaloHitList &collectedHits) const;

    /**
     *  @brief  Use two views of a delta ray pfo to calculate projected positions in a given third view
     *
     *  @param  element the tensor element
     *  @param  hitType the view to be projected into
     *  @param  projectedPositions the output list of projected positions
     *
     *  @return  a status code reflecting whether the procedure ran smoothly and if the outcome is good
     */
    pandora::StatusCode ProjectDeltaRayPositions(
        const TensorType::Element &element, const pandora::HitType hitType, pandora::CartesianPointVector &projectedPositions) const;

    float m_minSeparation;  ///< The minimum delta ray - parent muon cluster separation required to investigate a delta/cosmic ray cluster
    float m_distanceToLine; ///< The maximum perpendicular distance of a position to a line for it to be considered close
};

} // namespace lar_content

#endif // #ifndef REMOVAL_BASE_TOOL_H
