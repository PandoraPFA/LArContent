/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/CosmicRayRemovalTool.h
 *
 *  @brief  Header file for the cosmic ray removal tool class.
 *
 *  $Log: $
 */

#ifndef COSMIC_RAY_REMOVAL_TOOL_H
#define COSMIC_RAY_REMOVAL_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/RemovalBaseTool.h"
#include "larpandoracontent/LArThreeDReco/LArCosmicRay/ThreeViewDeltaRayMatchingAlgorithm.h"

namespace lar_content
{
/**
 *  @brief  CosmicRayRemovalTool class
 */
class CosmicRayRemovalTool : public RemovalBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    CosmicRayRemovalTool();

private:
    bool Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor);    
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Remove hits from delta ray clusters that belong to the parent muon
     *
     *  @param  elementList the tensor element list
     *  @param  changesMade receive boolean indicating whether clusters in element have been modified
     */
    void ExamineConnectedElements(const TensorType::ElementList &elementList, bool &changesMade) const;

    /**
     *  @brief  Determine whether element satifies simple checks
     *
     *  @param  element the tensor element
     *  @param  hitType the hit type of the cluster under investigation
     *
     *  @return  whether the checks pass
     */    
    virtual bool PassElementChecks(const TensorType::Element &element, const pandora::HitType &hitType) const;

    /**
     *  @brief  Determine whether the cluster under investigation has muon contamination
     *
     *  @param  element the tensor element
     *  @param  hitType the hit type of the cluster under investigation
     *
     *  @return  whether the cluster contains muon hits to remove
     */        
    bool IsContaminated(const TensorType::Element &element, const pandora::HitType &hitType) const;

    /**
     *  @brief  Collect hits in the delta ray cluster that lie close to calculated projected positions forming a seed to later grow
     *
     *  @param  element the tensor element
     *  @param  hitType the hit type of the cluster under investigation
     *  @param  collectedHits the output list of identified delta ray hits
     */
    void CreateSeed(const TensorType::Element &element, const pandora::HitType &hitType, pandora::CaloHitList &collectedHits) const;

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
    pandora::StatusCode GrowSeed(const TensorType::Element &element, const pandora::HitType &hitType, pandora::CaloHitList &collectedHits,
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
    void CollectHitsFromDeltaRay(const pandora::CartesianVector &positionOnMuon, const pandora::CartesianVector &muonDirection, const pandora::Cluster *const pDeltaRayCluster,
        const float &minDistanceFromMuon, const bool demandCloserToCollected, const pandora::CaloHitList &protectedHits, pandora::CaloHitList &collectedHits) const;

    /**
     *  @brief  Fragment the delta ray cluster adding hits to the muon track and perhaps creating significant remnants in the process
     *
     *  @param  element the tensor element
     *  @param  hitType the hit type of the cluster under investigation
     *  @param  collectedHits the list of identified delta ray hits
     *  @param  deltaRayRemantHits the list of remainder hits 
     */       
    void SplitDeltaRayCluster(const TensorType::Element &element, const pandora::HitType &hitType, pandora::CaloHitList &collectedHits,
        pandora::CaloHitList &deltaRayRemnantHits) const;

    /**
     *  @brief  Reculster a given remnant cluster, merging insignificant created clusters into the parent muon (if proximity checks pass)
     *
     *  @param  hitType the hit type of the cluster under investigation
     *  @param  pMuonCluster the address of the parent muon cluster
     *  @param  pDeltaRayRemnant the address of the delta ray remnant
     *  @param  clusterVector a vector containing the address of created/modified clusters for bookeeping purposes
     *  @param  pfoVector a vector containing the address of the pfo to which a modified muon cluster belongs for bookeeping purposes
     */        
    void FragmentRemnant(const pandora::HitType &hitType, const pandora::Cluster *const pMuonCluster,
        const pandora::Cluster *const pDeltaRayRemnant, pandora::ClusterVector &clusterVector, pandora::PfoVector &pfoVector) const;
    
    float m_minSeparation;                    ///< The minimum delta ray - parent muon cluster separation required to investigate delta ray cluster
    unsigned int m_slidingFitWindow;          ///< The sliding fit window used in cosmic ray parameterisations
    float m_minContaminationLength;           ///< The minimum projected length of a delta ray cluster onto the muon track for it to be considered contaminated
    float m_maxDistanceToHit;                 ///< The maximum distance of a hit from the cosmic ray track for it to be added into the CR
    unsigned int m_minRemnantClusterSize;     ///< The minimum hit number of a remnant cluster for it to be considered significant
    float m_maxDistanceToTrack;               ///< The minimum distance of an insignificant remnant cluster from the cosmic ray track for it to be merged into the CR
};

} // namespace lar_content

#endif // #ifndef COSMIC_RAY_REMOVAL_TOOL_H
