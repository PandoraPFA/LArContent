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
     *  @param  pAlgorithm address of the calling algorithm 
     *  @param  elementList the tensor element list
     *  @param  changesMade receive boolean indicating whether clusters in element have been modified
     */
    void RemoveMuonHits(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const MatrixType::ElementList &elementList, bool &changesMade) const;

    /**
     *  @brief  Determine whether element satifies simple checks
     *
     *  @param  element the tensor element
     *  @param  hitType the hit type of the cluster under investigation
     *
     *  @return  whether the checks pass
     */    
    bool PassElementChecks(const MatrixType::Element &element, const pandora::HitType &hitType) const;    

    /**
     *  @brief  Determine whether the cluster under investigation has muon contamination
     *
     *  @param  element the tensor element
     *  @param  hitType the hit type of the cluster under investigation
     *
     *  @return  whether the cluster contains muon hits to remove
     */        
    bool IsContaminated(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const MatrixType::Element &element, const pandora::HitType &hitType) const;

    /**
     *  @brief  Collect hits in the delta ray cluster that lie close to calculated projected positions forming a seed to later grow
     *
     *  @param  element the tensor element
     *  @param  hitType the hit type of the cluster under investigation
     *  @param  collectedHits the output list of identified delta ray hits
     */
    void CreateSeed(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const MatrixType::Element &element, const pandora::HitType &hitType,
        pandora::CaloHitList &collectedHits) const;

    /**
     *  @brief  Examine remaining hits in the delta ray cluster adding them to the delta ray seed or parent muon if appropriate and a separate list if not
     *
     *  @param  element the tensor element
     *  @param  hitType the hit type of the cluster under investigation
     *  @param  collectedHits the list of identified delta ray hits
     *  @param  deltaRayRemantHits the list of remainder hits 
     *
     *  @return  whether the muon projection mechanics were successful - abort process if not
     */   
    pandora::StatusCode GrowSeed(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const MatrixType::Element &element, const pandora::HitType &hitType,
        pandora::CaloHitList &collectedHits, pandora::CaloHitList &deltaRayRemantHits) const;

    /**
     *  @brief  Fragment the delta ray cluster, refining the matched delta ray cluster perhaps creating significant remnants in the process
     *
     *  @param  pAlgorithm address of the calling algorithm 
     *  @param  element the tensor element
     *  @param  hitType the hit type of the cluster under investigation
     *  @param  collectedHits the list of identified delta ray hits
     *  @param  deltaRayRemantHits the list of remainder hits 
     */       
    void SplitCluster(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const MatrixType::Element &element, const pandora::HitType &hitType,
        pandora::CaloHitList &collectedHits, pandora::CaloHitList &deltaRayRemnantHits) const;

    /**
     *  @brief  Reculster the remnant cluster, merging insignificant created clusters into the parent muon (if proximity checks pass)
     *
     *  @param  pAlgorithm address of the calling algorithm 
     *  @param  hitType the hit type of the cluster under investigation
     *  @param  pMuonCluster the address of the parent muon cluster
     *  @param  pDeltaRayRemnant the address of the delta ray remnant
     *  @param  clusterVector a vector containing the address of created/modified clusters for bookeeping purposes
     *  @param  pfoVector a vector containing the address of the pfo to which a modified muon cluster belongs for bookeeping purposes
     */        
    void FragmentRemnant(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const pandora::HitType &hitType, const pandora::Cluster *const pMuonCluster,
        const pandora::Cluster *const pDeltaRayRemnant, pandora::ClusterVector &clusterVector, pandora::PfoVector &pfoVector) const;

        /**
     *  @brief  Return
     *
     *  @param  pAlgorithm address of the calling algorithm 
     *  @param  elementList the tensor element list
     *  @param  changesMade receive boolean indicating whether clusters in element have been modified
     */    
    pandora::StatusCode GetMuonCluster(const MatrixType::Element &element, const pandora::HitType &hitType, const pandora::Cluster *&pMuonCluster) const;    

    /**
     *  @brief  Determine whether the input element is the best to use to modify the contaminated cluster (best is defined by the total hit count)
     *
     *  @param  element the tensor element
     *  @param  hitType the hit type of the cluster under investigation
     *  @param  elementList the tensor element list
     *
     *  @return  whether the input element is the best element
     */
    bool IsBestElement(const MatrixType::Element &element, const pandora::HitType &hitType, const MatrixType::ElementList &elementList) const;

   
    bool IsCloseToLine(const pandora::CartesianVector &hitPosition, const pandora::CartesianVector &lineStart, const pandora::CartesianVector &lineEnd,
        const float distanceToLine) const;     

    void FindExtrapolatedHits(const pandora::Cluster *const pCluster, const pandora::CartesianVector &lowerBoundary, const pandora::CartesianVector &upperBoundary,
	    pandora::CaloHitList &collectedHits) const;

    bool IsInLineSegment(const pandora::CartesianVector &lowerBoundary, const pandora::CartesianVector &upperBoundary, const pandora::CartesianVector &point) const;

    pandora::StatusCode ProjectDeltaRayPositions(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const MatrixType::Element &element,
        const pandora::HitType &hitType, pandora::CartesianPointVector &projectedPositions) const;

    const pandora::Cluster *GetCluster(const MatrixType::Element &element, const pandora::HitType &hitType) const;

    float m_minSeparation;    ///< The minimum delta ray - parent muon cluster separation
};

} // namespace lar_content

#endif // #ifndef TWO_VIEW_COSMIC_RAY_REMOVAL_TOOL_H
