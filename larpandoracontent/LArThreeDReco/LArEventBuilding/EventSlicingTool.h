/**
 *  @file   LArContent/LArThreeDReco/LArEventBuilding/EventSlicingTool.h
 * 
 *  @brief  Header file for the event slicing tool class.
 * 
 *  $Log: $
 */
#ifndef LAR_EVENT_SLICING_TOOL_H
#define LAR_EVENT_SLICING_TOOL_H 1

#include "larpandoracontent/LArContent/LArUtility/NeutrinoParentAlgorithm.h"

#include <unordered_map>

namespace lar_content
{

template<typename, unsigned int> class KDTreeLinkerAlgo;
template<typename, unsigned int> class KDTreeNodeInfoT;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  EventSlicingTool class
 */
class EventSlicingTool : public SlicingTool
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm tool
     */
    class Factory : public pandora::AlgorithmToolFactory
    {
    public:
        pandora::AlgorithmTool *CreateAlgorithmTool() const;
    };

    /**
     *  @brief  Default constructor
     */
    EventSlicingTool();

    typedef NeutrinoParentAlgorithm::SliceList SliceList;
    typedef NeutrinoParentAlgorithm::HitTypeToNameMap HitTypeToNameMap;

    void Slice(const NeutrinoParentAlgorithm *const pAlgorithm, const HitTypeToNameMap &caloHitListNames, const HitTypeToNameMap &clusterListNames,
        SliceList &sliceList);

private:
    typedef std::unordered_map<const pandora::Cluster*, const pandora::ParticleFlowObject*> ClusterToPfoMap;

    /**
     *  @brief  Get the 3D clusters from a specified list of pfos, storing the 3D clusters in the provided list and populating
     *          a map from 3D cluster to parent pfo
     *
     *  @param  pAlgorithm the address of the parent algorithm
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
    void GetClusterSliceList(const pandora::ClusterList &trackClusters3D, const pandora::ClusterList &showerClusters3D,
        ClusterSliceList &clusterSliceList) const;

    /**
     *  @brief  Compare the provided cluster slice and list of 3D track clusters. Add next appropriate (pointing) cluster to the slice.
     *
     *  @param  trackCandidates the sorted list of 3D candidate track clusters
     *  @param  slidingFitResultMap the map from 3D track clusters to 3D sliding fit results
     *  @param  clusterSlice the cluster slice
     *  @param  usedClusters the list of clusters already added to slices
     * 
     *  @return whether an addition to the cluster slice has been made
     */
    bool AddNextPointing(const pandora::ClusterVector &trackCandidates, const ThreeDSlidingFitResultMap &slidingFitResultMap,
        pandora::ClusterVector &clusterSlice, pandora::ClusterList &usedClusters) const;

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

    /**
     *  @brief  Compare the provided cluster slice and list of 3D clusters. Add next appropriate (nearby) cluster to the slice.
     *
     *  @param  clusterCandidates the sorted list of 3D candidate clusters
     *  @param  clusterSlice the cluster slice
     *  @param  usedClusters the list of clusters already added to slices
     * 
     *  @return whether an addition to the cluster slice has been made
     */
    bool AddNextProximity(const pandora::ClusterVector &clusterCandidates, pandora::ClusterVector &clusterSlice,
        pandora::ClusterList &usedClusters) const;

    /**
     *  @brief  Check separation of hits in two 3D clusters
     *
     *  @param  pCluster1 the address of the first cluster
     *  @param  pCluster2 the address of the second cluster
     * 
     *  @return whether the clusters are declared to be in the same slice
     */
    bool CheckHitSeparation(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2) const;

    typedef std::unordered_map<const pandora::Cluster*, unsigned int> ClusterToSliceIndexMap;

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
    void CopyPfoHitsToSlices(const ClusterToSliceIndexMap &clusterToSliceIndexMap, const ClusterToPfoMap &clusterToPfoMap, SliceList &sliceList,
        pandora::ClusterList &assignedClusters) const;

    /**
     *  @brief  Get the list of 2D clusters with hits yets to be assigned to slices
     *
     *  @param  pAlgorithm the address of the parent algorithm
     *  @param  clusterListNames the hit type to cluster list name map
     *  @param  assignedClusters the list of 2D clusters with hits assigned to slices
     *  @param  remainingClusters to receive the list of 2D clusters with hits yet to be assigned to slices
     */
    void GetRemainingClusters(const pandora::Algorithm *const pAlgorithm, const HitTypeToNameMap &clusterListNames,
        const pandora::ClusterList &assignedClusters, pandora::ClusterList &remainingClusters) const;

    /**
     *  @brief  Get the list of 2D clusters (from a named 2D cluster list) with hits yets to be assigned to slices
     *
     *  @param  pAlgorithm the address of the parent algorithm
     *  @param  clusterListName the cluster list name
     *  @param  assignedClusters the list of 2D clusters with hits assigned to slices
     *  @param  remainingClusters to receive the list of 2D clusters with hits yet to be assigned to slices
     */
    void GetRemainingClusters(const pandora::Algorithm *const pAlgorithm, const std::string &clusterListName,
        const pandora::ClusterList &assignedClusters, pandora::ClusterList &remainingClusters) const;

    /**
     *  @brief  Use the list of remaining 2D clusters to assign all remaining 2D hits to existing slices in the slice list
     *
     *  @param  remainingClusters the list of 2D clusters with hits yet to be assigned to slices
     *  @param  clusterToSliceIndexMap the mapping from 3D clusters to index in the slice list
     *  @param  sliceList the list containing slices to be populated with 2D hits
     */
    void AssignRemainingHitsToSlices(const pandora::ClusterList &remainingClusters, const ClusterToSliceIndexMap &clusterToSliceIndexMap,
        SliceList &sliceList) const;

    typedef KDTreeLinkerAlgo<const pandora::CartesianVector*, 2> PointKDTree2D;
    typedef KDTreeNodeInfoT<const pandora::CartesianVector*, 2> PointKDNode2D;
    typedef std::vector<PointKDNode2D> PointKDNode2DList;

    typedef std::unordered_set<const pandora::CartesianVector*> PointList;
    typedef std::unordered_map<const pandora::CartesianVector*, unsigned int> PointToSliceIndexMap;

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

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string     m_trackPfoListName;                 ///< The name of the input track pfo list
    std::string     m_showerPfoListName;                ///< The name of the input shower pfo list

    unsigned int    m_halfWindowLayers;                 ///< The number of layers to use for half-window of sliding fit

    float           m_minVertexLongitudinalDistance;    ///< Pointing association check: min longitudinal distance cut
    float           m_maxVertexLongitudinalDistance;    ///< Pointing association check: max longitudinal distance cut
    float           m_maxVertexTransverseDistance;      ///< Pointing association check: max transverse distance cut
    float           m_vertexAngularAllowance;           ///< Pointing association check: pointing angular allowance in degrees

    float           m_maxClosestApproach;               ///< Pointing association: max distance of closest approach between straight line fits
    float           m_maxInterceptDistance;             ///< Pointing association: max distance from cluster vertex to point of closest approach

    float           m_maxHitSeparationSquared;          ///< Proximity association: max distance allowed between the closest pair of hits

    bool            m_use3DProjectionsInHitPickUp;      ///< Whether to include 3D cluster projections when assigning remaining clusters to slices
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *EventSlicingTool::Factory::CreateAlgorithmTool() const
{
    return new EventSlicingTool();
}

} // namespace lar_content

#endif // #ifndef LAR_EVENT_SLICING_TOOL_H
