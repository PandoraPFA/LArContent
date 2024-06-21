/**
 *  @file   larpandoracontent/LArThreeDReco/LArEventBuilding/EventSlicingTool.h
 *
 *  @brief  Header file for the event slicing tool class.
 *
 *  $Log: $
 */
#ifndef LAR_EVENT_SLICING_TOOL_H
#define LAR_EVENT_SLICING_TOOL_H 1

#include "larpandoracontent/LArControlFlow/EventSlicingBaseTool.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingConeFitResult.h"

#include <unordered_map>

namespace lar_content
{

typedef std::map<pandora::HitType, std::string> HitTypeToNameMap;

template <typename, unsigned int>
class KDTreeLinkerAlgo;
template <typename, unsigned int>
class KDTreeNodeInfoT;

class SimpleCone;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  EventSlicingTool class
 */
class EventSlicingTool : public EventSlicingBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    EventSlicingTool();

    void RunSlicing(const pandora::Algorithm *const pAlgorithm, const HitTypeToNameMap &caloHitListNames,
        const HitTypeToNameMap &clusterListNames, SliceList &sliceList);

private:
    /**
     *  @brief  Copy all the input hits in an event into a single slice
     *
     *  @param  pAlgorithm the address of the calling algorithm
     *  @param  caloHitListNames the hit type to calo hit list name map
     *  @param  sliceList the slice list to receive the single new slice
     */
    void CopyAllHitsToSingleSlice(const pandora::Algorithm *const pAlgorithm, const HitTypeToNameMap &caloHitListNames, SliceList &sliceList) const;

    typedef std::unordered_map<const pandora::Cluster *, const pandora::ParticleFlowObject *> ClusterToPfoMap;

    /**
     *  @brief  Get the 3D clusters from a specified list of pfos, storing the 3D clusters in the provided list and populating
     *          a map from 3D cluster to parent pfo
     *
     *  @param  pAlgorithm the address of the calling algorithm
     *  @param  pfoListName the pfo list name
     *  @param  clusters3D to receive the list of 3D clusters
     *  @param  clusterToPfoMap to receive the mapping from 3D clusters to parent pfos
     */
    void GetThreeDClusters(const pandora::Algorithm *const pAlgorithm, const std::string &pfoListName, pandora::ClusterList &clusters3D,
        ClusterToPfoMap &clusterToPfoMap) const;

    typedef std::vector<pandora::ClusterVector> ClusterSliceList;

    /**
     *  @brief  Divide the provided lists of 3D track and shower clusters into slices
     *
     *  @param  trackClusters3D the list of 3D track clusters
     *  @param  showerClusters3D the list of 3D shower clusters
     *  @param  clusterSliceList to receive the list of 3D clusters, divided into slices (one 3D cluster list per slice)
     */
    void GetClusterSliceList(
        const pandora::ClusterList &trackClusters3D, const pandora::ClusterList &showerClusters3D, ClusterSliceList &clusterSliceList) const;

    /**
     *  @brief  Collect all clusters associated with a provided cluster
     *
     *  @param  pClusterInSlice the address of the cluster already in a slice
     *  @param  candidateClusters the list of candidate clusters
     *  @param  trackFitResults the map of sliding fit results for track candidate clusters
     *  @param  showerConeFitResults the map of sliding const fit results for shower candidate clusters
     *  @param  clusterSlice the cluster slice
     *  @param  usedClusters the list of clusters already added to slices
     */
    void CollectAssociatedClusters(const pandora::Cluster *const pClusterInSlice, const pandora::ClusterVector &candidateClusters,
        const ThreeDSlidingFitResultMap &trackFitResults, const ThreeDSlidingConeFitResultMap &showerConeFitResults,
        pandora::ClusterVector &clusterSlice, pandora::ClusterSet &usedClusters) const;

    /**
     *  @brief  Compare the provided clusters to assess whether they are associated via pointing (checks association "both ways")
     *
     *  @param  pClusterInSlice address of a cluster already in the slice
     *  @param  pCandidateCluster address of the candidate cluster
     *  @param  trackFitResults the map of sliding fit results for track candidate clusters
     *
     *  @return whether an addition to the cluster slice should be made
     */
    bool PassPointing(const pandora::Cluster *const pClusterInSlice, const pandora::Cluster *const pCandidateCluster,
        const ThreeDSlidingFitResultMap &trackFitResults) const;

    /**
     *  @brief  Compare the provided clusters to assess whether they are associated via pointing
     *
     *  @param  pClusterInSlice address of a cluster already in the slice
     *  @param  pCandidateCluster address of the candidate cluster
     *
     *  @return whether an addition to the cluster slice should be made
     */
    bool PassProximity(const pandora::Cluster *const pClusterInSlice, const pandora::Cluster *const pCandidateCluster) const;

    /**
     *  @brief  Compare the provided clusters to assess whether they are associated via cone fits to the shower cluster (single "direction" check)
     *
     *  @param  pClusterInSlice address of a cluster already in the slice
     *  @param  pCandidateCluster address of the candidate cluster
     *  @param  showerConeFitResults the map of sliding cone fit results for shower candidate clusters
     *
     *  @return whether an addition to the cluster slice should be made
     */
    bool PassShowerCone(const pandora::Cluster *const pConeCluster, const pandora::Cluster *const pNearbyCluster,
        const ThreeDSlidingConeFitResultMap &showerConeFitResults) const;

    /**
     *  @brief  Check closest approach metrics for a pair of pointing clusters
     *
     *  @param  cluster1 the first pointing cluster
     *  @param  cluster2 the second pointing cluster
     *
     *  @return whether the pointing clusters are declared to be in the same slice
     */
    bool CheckClosestApproach(const LArPointingCluster &cluster1, const LArPointingCluster &cluster2) const;

    /**
     *  @brief  Check closest approach metrics for a pair of pointing cluster vertices
     *
     *  @param  vertex1 the first pointing cluster vertex
     *  @param  vertex2 the second pointing cluster vertex
     *
     *  @return whether the pointing clusters are declared to be in the same slice
     */
    bool CheckClosestApproach(const LArPointingCluster::Vertex &vertex1, const LArPointingCluster::Vertex &vertex2) const;

    /**
     *  @brief  Check whether a pair of pointing clusters are nodally associated
     *
     *  @param  cluster1 the first pointing cluster
     *  @param  cluster2 the second pointing cluster
     *
     *  @return whether the pointing clusters are declared to be in the same slice
     */
    bool IsNode(const LArPointingCluster &cluster1, const LArPointingCluster &cluster2) const;

    /**
     *  @brief  Check whether a pair of pointing clusters are consistent with an emission
     *
     *  @param  cluster1 the first pointing cluster
     *  @param  cluster2 the second pointing cluster
     *
     *  @return whether the pointing clusters are declared to be in the same slice
     */
    bool IsEmission(const LArPointingCluster &cluster1, const LArPointingCluster &cluster2) const;

    typedef std::unordered_map<const pandora::Cluster *, unsigned int> ClusterToSliceIndexMap;

    /**
     *  @brief  Create new slices for each of the groupings of 3D clusters in the provided cluster slice list
     *
     *  @param  clusterSliceList the list of 3D clusters, divided into slices (one 3D cluster list per slice)
     *  @param  sliceList the slice list to receive the new slices
     *  @param  clusterToSliceIndexMap to receive the mapping from 3D clusters to index in the slice list
     */
    void CreateSlices(const ClusterSliceList &clusterSliceList, SliceList &sliceList, ClusterToSliceIndexMap &clusterToSliceIndexMap) const;

    /**
     *  @brief  Use 3D clusters in the cluster slice list, find their parent pfos and assign all hits in all 2D clusters in the pfos
     *          to the relevant slice in the output slice list
     *
     *  @param  clusterToSliceIndexMap the mapping from 3D clusters to index in the slice list
     *  @param  clusterToPfoMap the mapping from 3D clusters to parent pfos
     *  @param  sliceList the list containing slices to be populated with 2D hits
     *  @param  assignedClusters to receive the list of 2D clusters with hits assigned to slices
     */
    void CopyPfoHitsToSlices(const ClusterToSliceIndexMap &clusterToSliceIndexMap, const ClusterToPfoMap &clusterToPfoMap,
        SliceList &sliceList, pandora::ClusterSet &assignedClusters) const;

    /**
     *  @brief  Get the list of 2D clusters with hits yets to be assigned to slices
     *
     *  @param  pAlgorithm the address of the calling algorithm
     *  @param  clusterListNames the hit type to cluster list name map
     *  @param  assignedClusters the list of 2D clusters with hits assigned to slices
     *  @param  remainingClusters to receive the list of 2D clusters with hits yet to be assigned to slices
     */
    void GetRemainingClusters(const pandora::Algorithm *const pAlgorithm, const HitTypeToNameMap &clusterListNames,
        const pandora::ClusterSet &assignedClusters, pandora::ClusterList &remainingClusters) const;

    /**
     *  @brief  Get the list of 2D clusters (from a named 2D cluster list) with hits yets to be assigned to slices
     *
     *  @param  pAlgorithm the address of the calling algorithm
     *  @param  clusterListName the cluster list name
     *  @param  assignedClusters the list of 2D clusters with hits assigned to slices
     *  @param  remainingClusters to receive the list of 2D clusters with hits yet to be assigned to slices
     */
    void GetRemainingClusters(const pandora::Algorithm *const pAlgorithm, const std::string &clusterListName,
        const pandora::ClusterSet &assignedClusters, pandora::ClusterList &remainingClusters) const;

    /**
     *  @brief  Use the list of remaining 2D clusters to assign all remaining 2D hits to existing slices in the slice list
     *
     *  @param  remainingClusters the list of 2D clusters with hits yet to be assigned to slices
     *  @param  clusterToSliceIndexMap the mapping from 3D clusters to index in the slice list
     *  @param  sliceList the list containing slices to be populated with 2D hits
     */
    void AssignRemainingHitsToSlices(
        const pandora::ClusterList &remainingClusters, const ClusterToSliceIndexMap &clusterToSliceIndexMap, SliceList &sliceList) const;

    typedef KDTreeLinkerAlgo<const pandora::CartesianVector *, 2> PointKDTree2D;
    typedef KDTreeNodeInfoT<const pandora::CartesianVector *, 2> PointKDNode2D;
    typedef std::vector<PointKDNode2D> PointKDNode2DList;

    typedef std::list<const pandora::CartesianVector *> PointList;
    typedef std::unordered_map<const pandora::CartesianVector *, unsigned int> PointToSliceIndexMap;

    /**
     *  @brief  Use projections of 3D hits already assigned to slices to populate kd trees to aid assignment of remaining clusters
     *
     *  @param  sliceList the slice list
     *  @param  pointsU to receive the points in the u view
     *  @param  pointsV to receive the points in the v view
     *  @param  pointsW to receive the points in the w view
     *  @param  pointToSliceIndexMap to receive the mapping from points to slice index
     */
    void GetKDTreeEntries2D(const SliceList &sliceList, PointList &pointsU, PointList &pointsV, PointList &pointsW,
        PointToSliceIndexMap &pointToSliceIndexMap) const;

    /**
     *  @brief  Use 2D hits already assigned to slices to populate kd trees to aid assignment of remaining clusters
     *
     *  @param  clusterToSliceIndexMap the 3D cluster to slice index map
     *  @param  pointsU to receive the points in the u view
     *  @param  pointsV to receive the points in the v view
     *  @param  pointsW to receive the points in the w view
     *  @param  pointToSliceIndexMap to receive the mapping from points to slice index
     */
    void GetKDTreeEntries3D(const ClusterToSliceIndexMap &clusterToSliceIndexMap, PointList &pointsU, PointList &pointsV,
        PointList &pointsW, PointToSliceIndexMap &pointToSliceIndexMap) const;

    /**
     *  @brief  Use the provided kd tree to efficiently identify the most appropriate slice for the provided 2D cluster
     *
     *  @param  pCluster2D the address of the 2D cluster
     *  @param  kdTree the kd tree
     *
     *  @return the nearest-neighbour point identified by the kd tree
     */
    const PointKDNode2D *MatchClusterToSlice(const pandora::Cluster *const pCluster2D, PointKDTree2D &kdTree) const;

    /**
     *  @brief  Sort points (use Z, followed by X, followed by Y)
     *
     *  @param  pLhs address of first point
     *  @param  pRhs address of second point
     */
    static bool SortPoints(const pandora::CartesianVector *const pLhs, const pandora::CartesianVector *const pRhs);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_trackPfoListName;  ///< The name of the input track pfo list
    std::string m_showerPfoListName; ///< The name of the input shower pfo list

    unsigned int m_minHitsPer3DCluster;     ///< The minimum number of hits in a 3D cluster to warrant consideration in slicing
    unsigned int m_min3DHitsToSeedNewSlice; ///< The minimum number of hits in a 3D cluster to seed a new slice
    unsigned int m_halfWindowLayers;        ///< The number of layers to use for half-window of sliding fit

    bool m_usePointingAssociation;         ///< Whether to use pointing association
    float m_minVertexLongitudinalDistance; ///< Pointing association check: min longitudinal distance cut
    float m_maxVertexLongitudinalDistance; ///< Pointing association check: max longitudinal distance cut
    float m_maxVertexTransverseDistance;   ///< Pointing association check: max transverse distance cut
    float m_vertexAngularAllowance;        ///< Pointing association check: pointing angular allowance in degrees
    float m_maxClosestApproach;            ///< Pointing association: max distance of closest approach between straight line fits
    float m_maxInterceptDistance;          ///< Pointing association: max distance from cluster vertex to point of closest approach

    bool m_useProximityAssociation;  ///< Whether to use proximity association
    float m_maxHitSeparationSquared; ///< Proximity association: max distance allowed between the closest pair of hits

    bool m_useShowerConeAssociation; ///< Whether to use shower cone association
    unsigned int m_nConeFitLayers;   ///< The number of layers over which to sum fitted direction to obtain cone fit
    unsigned int m_nConeFits;        ///< The number of cone fits to perform, spread roughly uniformly along the shower length
    float m_coneLengthMultiplier;    ///< The cone length multiplier to use when calculating bounded cluster fractions
    float m_maxConeLength;           ///< The maximum allowed cone length to use when calculating bounded cluster fractions
    float m_coneTanHalfAngle1;       ///< The cone tan half angle to use when calculating bounded cluster fractions 1
    float m_coneBoundedFraction1;    ///< The minimum cluster bounded fraction for association 1
    float m_coneTanHalfAngle2;       ///< The cone tan half angle to use when calculating bounded cluster fractions 2
    float m_coneBoundedFraction2;    ///< The minimum cluster bounded fraction for association 2

    bool m_use3DProjectionsInHitPickUp; ///< Whether to include 3D cluster projections when assigning remaining clusters to slices
};

} // namespace lar_content

#endif // #ifndef LAR_EVENT_SLICING_TOOL_H
