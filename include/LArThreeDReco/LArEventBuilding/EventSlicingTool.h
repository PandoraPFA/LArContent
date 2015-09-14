/**
 *  @file   LArContent/include/LArThreeDReco/LArEventBuilding/EventSlicingTool.h
 * 
 *  @brief  Header file for the event slicing tool class.
 * 
 *  $Log: $
 */
#ifndef LAR_EVENT_SLICING_TOOL_H
#define LAR_EVENT_SLICING_TOOL_H 1

#include "LArUtility/NeutrinoParentAlgorithm.h"

#include <unordered_map>

namespace lar_content
{

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

    typedef NeutrinoParentAlgorithm::SliceList SliceList;
    typedef NeutrinoParentAlgorithm::HitTypeToNameMap HitTypeToNameMap;

    void Slice(NeutrinoParentAlgorithm *const pAlgorithm, const HitTypeToNameMap &caloHitListNames, const HitTypeToNameMap &clusterListNames,
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

    typedef std::vector<pandora::ClusterList> ClusterSliceList;

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
     *  @brief  Copy all the input hits in an event into a single slice
     *
     *  @param  pAlgorithm the address of the parent algorithm
     *  @param  caloHitListNames the hit type to calo hit list name map
     *  @param  sliceList the slice list to receive the single new slice
     */
    void CopyAllHitsToSingleSlice(const pandora::Algorithm *const pAlgorithm, const HitTypeToNameMap &caloHitListNames,
        SliceList &sliceList) const;

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

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string     m_trackPfoListName;         ///< The name of the input track pfo list
    std::string     m_showerPfoListName;        ///< The name of the input shower pfo list
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *EventSlicingTool::Factory::CreateAlgorithmTool() const
{
    return new EventSlicingTool();
}

} // namespace lar_content

#endif // #ifndef LAR_EVENT_SLICING_TOOL_H
