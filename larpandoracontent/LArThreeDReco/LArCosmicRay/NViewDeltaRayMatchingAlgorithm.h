/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/NViewDeltaRayMatching.h
 *
 *  @brief  Header file for the n view delta ray matching class.
 *
 *  $Log: $
 */
#ifndef LAR_N_VIEW_DELTA_RAY_MATCHING_ALGORITHM_H
#define LAR_N_VIEW_DELTA_RAY_MATCHING_ALGORITHM_H 1

#include "larpandoracontent/LArThreeDReco/LArThreeDBase/NViewMatchingAlgorithm.h"

#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

namespace lar_content
{
/**
 *  @brief  NViewDeltaRayMatchingAlgorithm class
 */
template<typename T>    
class NViewDeltaRayMatchingAlgorithm : public NViewMatchingAlgorithm<T>
{
public:
    typedef std::map<const pandora::CaloHit*, const pandora::Cluster*> HitToClusterMap;
    typedef std::map<const pandora::Cluster*, const pandora::ParticleFlowObject*> ClusterToPfoMap;
    typedef std::map<const pandora::Cluster*, pandora::ClusterList> ClusterProximityMap;

    typedef KDTreeLinkerAlgo<const pandora::CaloHit*, 2> HitKDTree2D;
    typedef KDTreeNodeInfoT<const pandora::CaloHit*, 2> HitKDNode2D;
    typedef std::vector<HitKDNode2D> HitKDNode2DList;

    /**
     *  @brief  Default constructor
     */
    NViewDeltaRayMatchingAlgorithm();
    
    /**
     *  @brief  To check whether a given cluster meets the requirements to be added into the matching container (tensor/matrix)
     *
     *  @param  pCluster the address of the input cluster
     *
     *  @return  whether the checks were met
     */
    virtual bool DoesClusterPassTesorThreshold(const pandora::Cluster *const pCluster) const = 0;

    /**
     *  @brief  Return the cluster of the common cosmic ray pfo in a given view (function demands there to be only one common CR pfo) 
     *
     *  @param  commonMuonPfoList the element's list of common muon pfos
     *  @param  hitType the specified view
     *  @param  pMuonCluster the output address of the cluster
     *
     *  @return  a status code reflecting if one and only one cosmic ray cluster was found
     */
    pandora::StatusCode GetMuonCluster(const pandora::PfoList &commonMuonPfoList, const pandora::HitType &hitType, const pandora::Cluster *&pMuonCluster) const;

    /**
     *  @brief  To determine how well three clusters (one in each view) map onto one another expressing this in terms of a chi-squared like parameter
     *
     *  @param  pCluster1 the first cluster
     *  @param  pCluster2 the second cluster
     *  @param  pCluster3 the third cluster
     *  @param  reducedChiSquared the reduced chi squared
     *
     *  @return  a status code reflecting whether the matching procedure ran smoothly and if the outcome is good 
     */
    pandora::StatusCode PerformThreeViewMatching(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, const pandora::Cluster *const pCluster3,
        float &reducedChiSquared) const;

    /**
     *  @brief  To determine how well three clusters (one in each view) map onto one another expressing this in terms of a chi-squared like parameter
     *
     *  @param  pClusterU the U cluster (if the xOverlap object is to be retained this must be the u cluster - for labels to make sense)
     *  @param  pClusterV the V cluster
     *  @param  pClusterW the W cluster
     *  @param  chiSquaredSum the sum of the chi-squared values of the sampled points
     *  @param  nSamplingPoints the number of sampled points
     *  @param  nMatchedSamplingPoints the number of matched sampled points
     *  @param  xOverlap the xOverlap object  
     *
     *  @return  a status code reflecting whether the matching procedure ran smoothly and if the outcome is good 
     */
    pandora::StatusCode PerformThreeViewMatching(const pandora::Cluster *const pClusterU, const pandora::Cluster *const pClusterV, const pandora::Cluster *const pClusterW,
        float &chiSquaredSum, unsigned int &nSamplingPoints, unsigned int &nMatchedSamplingPoints, XOverlap &XOverlap) const;

    /**
     *  @brief  To determine how well three CaloHitLists (one for each view) map onto one another expressing this in terms of a chi-squared like parameter
     *
     *  @param  clusterU the U CaloHitList (if the xOverlap object is to be retained this must be the u cluster - for labels to make sense)
     *  @param  clusterV the V CaloHitList
     *  @param  clusterW the W CaloHitList
     *  @param  chiSquaredSum the sum of the chi-squared values of the sampled points
     *  @param  nSamplingPoints the number of sampled points
     *  @param  nMatchedSamplingPoints the number of matched sampled points
     *  @param  xOverlap the xOverlap object  
     *
     *  @return  a status code reflecting whether the matching procedure ran smoothly and if the outcome is good 
     */
    pandora::StatusCode PerformThreeViewMatching(const pandora::CaloHitList &clusterU, const pandora::CaloHitList &clusterV, const pandora::CaloHitList &clusterW,
        float &chiSquaredSum, unsigned int &nSamplingPoints, unsigned int &nMatchedSamplingPoints, XOverlap &XOverlap) const;

    /**
     *  @brief  Use two views of a cosmic ray pfo to calculate projected positions in a given the third view
     *
     *  @param  thirdViewHitType the view to be projected into
     *  @param  pParentMuon the input cosmic ray pfo
     *  @param  projectedPositions the output projected positions
     *
     *  @return  a status code reflecting whether the procedure ran smoothly and if the outcome is good 
     */
    pandora::StatusCode ProjectMuonPositions(const pandora::HitType &thirdViewHitType, const pandora::ParticleFlowObject *const pParentMuon,
        pandora::CartesianPointVector &projectedPositions) const;

    /**
     *  @brief  Use two clusters from different views to calculate projected positions in the remaining third view
     *
     *  @param  pCluster1 the address of the first cluster
     *  @param  pCluster2 the address of the second cluster
     *  @param  projectedPositions the output projected positions
     *
     *  @return  a status code reflecting whether the procedure ran smoothly and if the outcome is good 
     */
    pandora::StatusCode GetProjectedPositions(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2,
        pandora::CartesianPointVector &projectedPositions) const;    

    /**
     *  @brief  In one view, pull out any hits from a cosmic ray cluster that belong to the child delta ray cluster 
     *
     *  @param  pCluster1 the address of a delta ray cluster in a view that is to go unmodified
     *  @param  pCluster2 the address of a delta ray cluster in the other view that is to unmodified 
     *  @param  pThirdViewCluster the address of the delta ray cluster in the view in which the hit removal process will run
     *  @param  pParentMuon the address of the parent cosmic ray pfo
     *  @param  minDistanceFromMuon the minimum distance of a hit from the cosmic ray track required for removal
     *  @param  maxDistanceToCollected the maximim distance of a hit from the projected delta ray hits required for removal
     *  @param  collectedHits the list of hits to be removed from the cosmic ray 
     *
     *  @return  a status code reflecting whether the procedure ran smoothly and if the outcome is good 
     */
    pandora::StatusCode CollectHitsFromMuon(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2,
        const pandora::Cluster *const pThirdViewCluster, const pandora::ParticleFlowObject *const pParentMuon, const float minDistanceFromMuon, 
        const float maxDistanceToCollected, pandora::CaloHitList &collectedHits) const;

    /**
     *  @brief  Parameterise the projection of a cosmic ray track in order to avoid poor/sparse projections
     *
     *  @param  pParentMuon the address of the cosmic ray pfo
     *  @param  pDeltaRayCluster the address of the delta ray cluster in the projected view
     *  @param  positionOnMuon the output position to localise the parameterisation in space
     *  @param  muonDirection the output cosmic ray direction
     *
     *  @return  a status code reflecting whether the procedure ran smoothly and if the outcome is good 
     */
    pandora::StatusCode ParameteriseMuon(const pandora::ParticleFlowObject *const pParentMuon, const pandora::Cluster *const pDeltaRayCluster,
        pandora::CartesianVector &positionOnMuon, pandora::CartesianVector &muonDirection) const;

    /**
     *  @brief  Parameterise the projection of a cosmic ray track in order to avoid poor/sparse projections
     *
     *  @param  pParentMuon the address of the cosmic ray pfo
     *  @param  deltaRayProjectedPositions the projected positions of the delta ray
     *  @param  hitType the view in which the projection is made
     *  @param  positionOnMuon the output position to localise the parameterisation in space
     *  @param  muonDirection the output cosmic ray direction
     *
     *  @return  a status code reflecting whether the procedure ran smoothly and if the outcome is good 
     */
    pandora::StatusCode ParameteriseMuon(const pandora::ParticleFlowObject *const pParentMuon, const pandora::CartesianPointVector &deltaRayProjectedPositions,
        const pandora::HitType &hitType, pandora::CartesianVector &positionOnMuon, pandora::CartesianVector &muonDirection) const;

    /**
     *  @brief  In one view, pull out any hits from a cosmic ray cluster that belong to the child delta ray cluster 
     *
     *  @param  positionOnMuon the parameterised cosmic ray position 
     *  @param  muonDirection the parameterised cosmic ray direction 
     *  @param  pMuon the address of the parent cosmic ray pfo
     *  @param  deltaRayProjectedPositions the projected positions of the delta ray
     *  @param  minDistanceFromMuon the minimum distance of a hit from the cosmic ray track required for removal
     *  @param  maxDistanceToCollected the maximim distance of a hit from the projected delta ray hits required for removal
     *  @param  collectedHits the list of hits to be removed from the cosmic ray 
     */
    void CollectHitsFromMuon(const pandora::CartesianVector &positionOnMuon, const pandora::CartesianVector &muonDirection,
        const pandora::Cluster *const pMuonCluster, const pandora::CartesianPointVector &deltaRayProjectedPositions, const float &minDistanceFromMuon,
        const float maxDistanceToCollected, pandora::CaloHitList &collectedHits) const;

    /**
     *  @brief  Move a list of hits from a cosmic ray cluster into the given child delta ray cluster
     *
     *  @param  clusterListName the pandora list to which the cosmic ray and delta ray clusters belong
     *  @param  pMuonCluster the address of the cosmic ray cluster
     *  @param  collectedHits the list of hits to reassign 
     *  @param  pDeltaRayCluster the address of the delta ray cluster (may be a nullptr if cluster is yet to be made)
     */
    void SplitMuonCluster(const std::string &clusterListName, const pandora::Cluster *const pMuonCluster, const pandora::CaloHitList &collectedHits,
        const pandora::Cluster *&pDeltaRayCluster) const;    

    /**
     *  @brief  Create delta ray pfos maxmising completeness by searching for and merging in any stray clusters
     *
     *  @param  protoParticleVector the proto particle vector
     *
     *  @return  whether any pfos were created
     */ 
    bool CreatePfos(ProtoParticleVector &protoParticleVector);

    /**
     *  @brief  Add a new cluster to algorithm ownership maps and, if it a delta ray cluster, to the underlying matches container (tensor/matrix)
     *
     *  @param  newClusterVector the vector of clusters to add - the order must match the pfoVector
     *  @param  pfoVector the vector of cosmic ray pfos to which the new clusters belong (nullptr for delta ray cluster)
     */ 
    void UpdateForNewClusters(const pandora::ClusterVector &newClusterVector, const pandora::PfoVector &pfoVector);

    /**
     *  @brief  Add a new cluster to algorithm ownership maps
     *
     *  @param  newClusterVector the vector of clusters to add - the order must match the pfoVector
     *  @param  pfoVector the vector of cosmic ray pfos to which the new clusters belong (nullptr for delta ray cluster)
     */     
    void UpdateContainers(const pandora::ClusterVector &newClusterVector, const pandora::PfoVector &pfoVector);
    
    void UpdateUponDeletion(const pandora::Cluster *const pDeletedCluster);
    void SelectInputClusters(const pandora::ClusterList *const pInputClusterList, pandora::ClusterList &selectedClusterList) const;
    void PrepareInputClusters(pandora::ClusterList &preparedClusterList);

protected:
    /**
     *  @brief  Populate the hit to cluster map  
     *
     *  @param  hitType the hit type of the map to populate
     */
    void FillHitToClusterMap(const pandora::HitType &hitType);

    /**
     *  @brief  Add the hits of a given cluster to the hit to cluster map
     *
     *  @param  hitType the hit type of the input cluster and of the map to add to
     */
    void AddToClusterMap(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Populate the cluster proximity map which maps clusters to neighbouring clusters
     *
     *  @param  hitType the hit type of the map to populate
     */
    void FillClusterProximityMap(const pandora::HitType &hitType);

    /**
     *  @brief  Build the KD tree
     *
     *  @param  hitType the hit type of the tree to build
     */
    void BuildKDTree(const pandora::HitType &hitType);

    /**
     *  @brief  Add a cluster to the cluster proximity map
     *
     *  @param  hitType the hit type of the input cluster and of the map to add to
     */
    void AddToClusterProximityMap(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Populate the cluster to pfo map which maps clusters to their cosmic ray pfo owner (if appropriate)
     *
     *  @param  hitType the hit type of the map to populate
     */
    void FillClusterToPfoMap(const pandora::HitType &hitType);

    /**
     *  @brief  Fill the stray cluster list with clusters that do not pass the tensor threshold requirement 
     *
     *  @param  hitType the hit type of the list to fill
     */
    void FillStrayClusterList(const pandora::HitType &hitType);
    
    /**
     *  @brief  Use the cluster proximity map to travel along paths of nearby clusters finding the cosmic ray clusters on which they terminate 
     *
     *  @param  pCluster the address of the input cluster
     *  @param  consideredClusters the list of investigated clusters
     *  @param  nearbyMuonPfos the output list of the cosmic ray pfos to which the nearby cosmic ray clusters belong
     */
    void GetNearbyMuonPfos(const pandora::Cluster *const pCluster, pandora::ClusterList &consideredClusters, pandora::PfoList &nearbyMuonPfos) const;

    /**
     *  @brief  Calculate the xSpan of a list of CaloHits
     *
     *  @param  caloHitList the input list of CaloHits
     *  @param  xMin the output minimum x coordinate
     *  @param  xMax the output maximum x coordinate
     */
    void GetClusterSpanX(const pandora::CaloHitList &caloHitList, float &xMin, float &xMax) const;

    /**
     *  @brief  Calculate the zSpan of a list of CaloHits in a specified x range
     *
     *  @param  caloHitList the input list of CaloHits
     *  @param  xMin the minimum x coordinate of the region of interest
     *  @param  xMax the maximum x coordinate of the region of interest
     *  @param  zMin the output minimum z coordinate
     *  @param  zMax the output maximum z coordinate
     */
    pandora::StatusCode GetClusterSpanZ(const pandora::CaloHitList &caloHitList, const float xMin, const float xMax, float &zMin, float &zMax) const;

    /**
     *  @brief  Collect the stray clusters that are close to a specified cluster and that lie within a given x range
     *
     *  @param  pClusterToEnlarge the specified cluster 
     *  @param  rangeMinX the minimum x coordinate of the region of interest
     *  @param  rangeMaxX the maximum x coordinate of the region of interest
     *  @param  collectedClusters the list of collected stray clusters
     */
    void CollectStrayClusters(const pandora::Cluster *const pClusterToEnlarge, const float rangeMinX, const float rangeMaxX, pandora::ClusterList &collectedClusters);

    /**
     *  @brief  Merge in the collected stray clusters of a given delta ray cluster
     *
     *  @param  pClusterToEnlarge the delta ray cluster to enlarge 
     *  @param  collectedClusters the list of collected stray clusters
     */
    void AddInStrayClusters(const pandora::Cluster *const pClusterToEnlarge, const pandora::ClusterList &collectedClusters);

    void TidyUp();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_muonPfoListName; ///< The list of reconstructed cosmic ray pfos
    HitToClusterMap m_hitToClusterMapU; ///< The mapping of hits to the clusters to which they belong (in the U view)
    HitToClusterMap m_hitToClusterMapV; ///< The mapping of hits to the clusters to which they belong (in the V view)
    HitToClusterMap m_hitToClusterMapW; ///< The mapping of hits to the clusters to which they belong (in the W view)
    HitKDTree2D m_kdTreeU; ///< The KD tree (in the U view)    
    HitKDTree2D m_kdTreeV; ///< The KD tree (in the V view)
    HitKDTree2D m_kdTreeW; ///< The KD tree (in the W view)
    ClusterProximityMap m_clusterProximityMapU; ///< The mapping of clusters to their neighbouring clusters (in the U view)
    ClusterProximityMap m_clusterProximityMapV; ///< The mapping of clusters to their neighbouting clusters (in the V view)
    ClusterProximityMap m_clusterProximityMapW; ///< The mapping of clusters to their neighbouring clusters (in the W view)
    ClusterToPfoMap m_clusterToPfoMapU; ///< The mapping of cosmic ray U clusters to the cosmic ray pfos to which they belong 
    ClusterToPfoMap m_clusterToPfoMapV; ///< The mapping of cosmic ray V clusters to the cosmic ray pfos to which they belong 
    ClusterToPfoMap m_clusterToPfoMapW; ///< The mapping of cosmic ray W clusters to the cosmic ray pfos to which they belong 
    pandora::ClusterList m_strayClusterListU; ///< The list of U clusters that do not pass the tensor threshold requirement 
    pandora::ClusterList m_strayClusterListV; ///< The list of V clusters that do not pass the tensor threshold requirement 
    pandora::ClusterList m_strayClusterListW; ///< The list of W clusters that do not pass the tensor threshold requirement     
    float m_searchRegion1D; ///< Search region, applied to each dimension, for look-up from kd-tree    
    float m_pseudoChi2Cut; ///< Pseudo chi2 cut for three view matching 
    float m_xOverlapWindow; ///< The maximum allowed displacement in x position
    float m_minMatchedFraction; ///< The threshold matched fraction of sampling points for a good match  
    unsigned int m_minMatchedPoints; ///< The threshold number of matched sampling points for a good match
    unsigned int m_minProjectedPositions; ///< The threshold number of projected points for a good projection
    float m_maxCosmicRayHitFraction; ///< The maximum allowed fraction of hits to be removed from the cosmic ray track
    float m_strayClusterSeparation; ///< The maximum allowed separation of a stray cluster and a delta ray cluster for merge
};
    
} // namespace lar_content

#endif // #ifndef LAR_N_VIEW_DELTA_RAY_MATCHING_ALGORITHM_H
