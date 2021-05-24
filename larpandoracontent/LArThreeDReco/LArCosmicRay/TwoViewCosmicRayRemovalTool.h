/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/TwoViewCosmicRayRemovalTool.h
 *
 *  @brief  Header file for the cosmic ray removal tool class
 *
 *  $Log: $
 */
#ifndef TWO_VIEW_COSMIC_RAY_REMOVAL_TOOL_H
#define TWO_VIEW_COSMIC_RAY_REMOVAL_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/TwoViewDeltaRayMatchingAlgorithm.h"

namespace lar_content
{
/**
 *  @brief  TwoViewCosmicRayRemovalTool class
 */
class TwoViewCosmicRayRemovalTool : public DeltaRayMatrixTool
{
public:
    typedef std::vector<pandora::HitType> HitTypeVector;
    /**
     *  @brief  Default constructor
     */
    TwoViewCosmicRayRemovalTool();

private:
    bool Run(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, MatrixType &overlapMatrix);
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Remove hits from delta ray clusters that belong to the parent muon
     *
     *  @param  elementList the matrix element list
     *
     *  @return  whether any clusters have been modified
     */
    bool RemoveCosmicRayHits(const MatrixType::ElementList &elementList) const;

    /**
     *  @brief  Determine whether element satifies simple checks
     *
     *  @param  element the matrix element
     *  @param  hitType the hit type of the cluster under investigation
     *
     *  @return  whether the checks pass
     */
    bool PassElementChecks(const MatrixType::Element &element, const pandora::HitType hitType) const;

    /**
     *  @brief  Determine whether the matched clusters suggest that the delta ray is at the endpoint of the cosmic ray (and is likely to be a michel)
     *
     *  @param  element the matrix element
     *  @param  ignoreHitType whether to ignore the cluster under investigation
     *  @param  hitTypeToIgnore the hit type to ignore
     *
     *  @return  whether the delta ray is at the endpoint of the cosmic ray
     */
    bool IsMuonEndpoint(const MatrixType::Element &element, const pandora::HitType hitType) const;

    /**
     *  @brief  Determine whether the cluster under investigation has muon contamination
     *
     *  @param  element the matrix element
     *  @param  hitType the hit type of the cluster under investigation
     *
     *  @return  whether the cluster contains muon hits to remove
     */
    bool IsContaminated(const MatrixType::Element &element, const pandora::HitType hitType) const;

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
     *  @brief  Determine whether the input element is the best to use to modify the contaminated cluster (best is defined by the total hit count)
     *
     *  @param  element the matrix element
     *  @param  hitType the hit type of the cluster under investigation
     *  @param  elementList the matrix element list
     *
     *  @return  whether the input element is the best element
     */
    bool IsBestElement(const MatrixType::Element &element, const pandora::HitType hitType, const MatrixType::ElementList &elementList) const;

    /**
     *  @brief  Collect hits in the delta ray cluster that lie close to calculated projected positions forming a seed to later grow
     *
     *  @param  element the matrix element
     *  @param  hitType the hit type of the cluster under investigation
     *  @param  collectedHits the output list of identified delta ray hits
     */
    void CreateSeed(const MatrixType::Element &element, const pandora::HitType hitType, pandora::CaloHitList &collectedHits) const;

    /**
     *  @brief  Examine remaining hits in the delta ray cluster adding them to the delta ray seed or parent muon if appropriate and a separate list if not
     *
     *  @param  element the matrix element
     *  @param  hitType the hit type of the cluster under investigation
     *  @param  collectedHits the list of identified delta ray hits
     *  @param  deltaRayRemantHits the list of remainder hits 
     *
     *  @return  whether the muon projection mechanics were successful - abort process if not
     */
    pandora::StatusCode GrowSeed(const MatrixType::Element &element, const pandora::HitType hitType, pandora::CaloHitList &collectedHits,
        pandora::CaloHitList &deltaRayRemantHits) const;

    /**
     *  @brief  Collect hits from the delta ray cluster to either keep (demandCloserToCollected == true) or separate into a new shower (demandCloserToCollected == false)
     *
     *  @param  positionOnMuon the parameterised muon position
     *  @param  muonDirection the parameterised muon direction
     *  @param  pDeltaRayCluster the delta ray cluster under investigation
     *  @param  minDistanceFromMuon the minimum distance of a hit from the muon track for it to not belong to the muon
     *  @param  demandCloserToCollected whether to demand a hit be closer to the collected hits than to the muon hits for it to be collected
     *  @param  protectedHits the hits that are protected from being collected
     *  @param  collectedHits the output list of collected hits
     */
    void CollectHitsFromDeltaRay(const pandora::CartesianVector &positionOnMuon, const pandora::CartesianVector &muonDirection,
        const pandora::Cluster *const pDeltaRayCluster, const float &minDistanceFromMuon, const bool demandCloserToCollected,
        const pandora::CaloHitList &protectedHits, pandora::CaloHitList &collectedHits) const;

    /**
     *  @brief  Fragment the delta ray cluster, refining the matched delta ray cluster perhaps creating significant remnants in the process
     *
     *  @param  element the matrix element
     *  @param  hitType the hit type of the cluster under investigation
     *  @param  collectedHits the list of identified delta ray hits
     *  @param  deltaRayRemantHits the list of remainder hits 
     */
    void SplitDeltaRayCluster(const MatrixType::Element &element, const pandora::HitType hitType, pandora::CaloHitList &collectedHits,
        pandora::CaloHitList &deltaRayRemnantHits) const;

    /**
     *  @brief  Reculster the remnant cluster, merging insignificant created clusters into the parent muon (if proximity checks pass)
     *
     *  @param  hitType the hit type of the cluster under investigation
     *  @param  pMuonCluster the address of the parent muon cluster
     *  @param  pDeltaRayRemnant the address of the delta ray remnant
     *  @param  clusterVector a vector containing the address of created/modified clusters for bookeeping purposes
     *  @param  pfoVector a vector containing the address of the pfo to which a modified muon cluster belongs for bookeeping purposes
     */
    void ReclusterRemnant(const pandora::HitType hitType, const pandora::Cluster *const pMuonCluster,
        const pandora::Cluster *const pDeltaRayRemnant, pandora::ClusterVector &clusterVector, pandora::PfoVector &pfoVector) const;

    /**
     *  @brief  Use two views of a delta ray pfo to calculate projected positions in a given third view
     *
     *  @param  element the matrix element
     *  @param  hitType the view to be projected into
     *  @param  projectedPositions the output list of projected positions
     *
     *  @return  a status code reflecting whether the procedure ran smoothly and if the outcome is good
     */
    pandora::StatusCode ProjectDeltaRayPositions(
        const MatrixType::Element &element, const pandora::HitType hitType, pandora::CartesianPointVector &projectedPositions) const;

    float m_minSeparation;           ///< The minimum delta ray - parent muon cluster separation required to investigate delta ray cluster
    unsigned int m_slidingFitWindow; ///< The sliding fit window used in cosmic ray parameterisations
    float m_distanceToLine;          ///< The maximum perpendicular distance of a position to a line for it to be considered close
    float m_minContaminationLength; ///< The minimum projected length of a delta ray cluster onto the muon track for it to be considered contaminated
    float m_maxDistanceToHit;             ///< The maximum distance of a hit from the cosmic ray track for it to be added into the CR
    unsigned int m_minRemnantClusterSize; ///< The minimum hit number of a remnant cluster for it to be considered significant
    float m_maxDistanceToTrack; ///< The minimum distance of an insignificant remnant cluster from the cosmic ray track for it to be merged into the CR
};

} // namespace lar_content

#endif // #ifndef TWO_VIEW_COSMIC_RAY_REMOVAL_TOOL_H
