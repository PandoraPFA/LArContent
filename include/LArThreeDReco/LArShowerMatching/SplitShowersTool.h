/**
 *  @file   LArContent/include/LArThreeDReco/LArShowerMatching/SplitShowersTool.h
 * 
 *  @brief  Header file for the split showers tool class.
 * 
 *  $Log: $
 */
#ifndef SPLIT_SHOWERS_TOOL_H
#define SPLIT_SHOWERS_TOOL_H 1

#include "LArThreeDReco/LArShowerMatching/ThreeDShowersAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  SplitShowersTool class
 */
class SplitShowersTool : public ShowerTensorTool
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
    SplitShowersTool();

    bool Run(ThreeDShowersAlgorithm *pAlgorithm, TensorType &overlapTensor);

private:
    /**
     *  @brief  Find split showers, using information from the overlap tensor
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  overlapTensor the overlap tensor
     *  @param  clusterMergeMap to receive the list of cluster merges
     */
    void FindSplitShowers(ThreeDShowersAlgorithm *pAlgorithm, const TensorType &overlapTensor, ClusterMergeMap &clusterMergeMap) const;

    /**
     *  @brief  Whether a provided (iterator to a) tensor element passes the selection cuts for undershoots identification
     * 
     *  @param  eIter the iterator to the tensor element
     *  @param  usedClusters the list of used clusters
     */
    bool PassesElementCuts(TensorType::ElementList::const_iterator eIter, const pandora::ClusterList &usedClusters) const;

    /**
     *  @brief  Select elements representing possible components of interest due to undershoots in clustering
     * 
     *  @param  eIter iterator to a candidate element
     *  @param  elementList the provided element list
     *  @param  usedClusters the list of used clusters
     *  @param  iteratorList to receive a list of iterators to relevant elements
     */
    void SelectTensorElements(TensorType::ElementList::const_iterator eIter, const TensorType::ElementList &elementList,
        const pandora::ClusterList &usedClusters, IteratorList &iteratorList) const;

    /**
     *  @brief  Get cluster merges specific elements of the tensor
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  iteratorList list of iterators to relevant tensor elements
     *  @param  usedClusters the list of used clusters
     *  @param  clusterMergeMap to be populated with cluster merges
     */
    void FindShowerMerges(ThreeDShowersAlgorithm *pAlgorithm, const IteratorList &iteratorList, pandora::ClusterList &usedClusters,
        ClusterMergeMap &clusterMergeMap) const;

    /**
     *  @brief  Check the clusters in a provided cluster list are in suitable proximity for merging
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  clusterList the cluster list
     */
    bool CheckClusterProximities(ThreeDShowersAlgorithm *pAlgorithm, const pandora::ClusterList &clusterList) const;

    /**
     *  @brief  Check the consistency of the clusters in a provided cluster list with the event vertex, if available
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  clusterList the cluster list
     */
    bool CheckClusterVertexRelations(ThreeDShowersAlgorithm *pAlgorithm, const pandora::ClusterList &clusterList) const;

    /**
     *  @brief  Check the consistency of the split positions in the provided u, v and w cluster lists
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  clusterListU the u cluster list
     *  @param  clusterListV the v cluster list
     *  @param  clusterListW the w cluster list
     */
    bool CheckClusterSplitPositions(ThreeDShowersAlgorithm *pAlgorithm, const pandora::ClusterList &clusterListU,
        const pandora::ClusterList &clusterListV, const pandora::ClusterList &clusterListW) const;

    /**
     *  @brief  Get the x coordinate representing the midpoint between two clusters (hypothesis: clusters represent a split shower)
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  pClusterA the address of cluster A
     *  @param  pClusterB the address of cluster B
     *  @param  splitXPosition to receive the split position estimate
     *  @param  overlapX to receive the overlap estimate
     */
    void GetSplitXDetails(ThreeDShowersAlgorithm *pAlgorithm, pandora::Cluster *const pClusterA, pandora::Cluster *const pClusterB,
        float &splitXPosition, float &overlapX) const;

    /**
     *  @brief  Populate the cluster merge map, based on the information contained in the provided cluster list
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  clusterList the cluster list
     *  @param  clusterMergeMap to receive the populated cluster merge map
     */
    void SpecifyClusterMerges(ThreeDShowersAlgorithm *pAlgorithm, const pandora::ClusterList &clusterList, ClusterMergeMap &clusterMergeMap) const;

    /**
     *  @brief  Apply the changes cached in a cluster merge map and update the tensor accordingly
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  clusterMergeMap the cluster merge map
     * 
     *  @return whether changes to the tensor have been made
     */
    bool ApplyChanges(ThreeDShowersAlgorithm *pAlgorithm, const ClusterMergeMap &clusterMergeMap) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int    m_nCommonClusters;                  ///< The number of common clusters
    float           m_minMatchedFraction;               ///< The min matched sampling point fraction for use as a key tensor element
    unsigned int    m_minMatchedSamplingPoints;         ///< The min number of matched sampling points for use as a key tensor element

    bool            m_checkClusterProximities;          ///< Whether to check the proximities of the candidate split shower clusters
    float           m_maxClusterSeparation;             ///< The maximum separation for clusters to be merged

    bool            m_checkClusterVertexRelations;      ///< Whether to check the consistency of the clusters with the event vertex
    float           m_minVertexLongitudinalDistance;    ///< Vertex association check: min longitudinal distance cut
    float           m_maxVertexLongitudinalDistance;    ///< Vertex association check: max longitudinal distance cut
    float           m_maxVertexTransverseDistance;      ///< Vertex association check: max transverse distance cut
    float           m_vertexAngularAllowance;           ///< Vertex association check: pointing angular allowance in degrees
    unsigned int    m_maxVertexAssociations;            ///< The maximum number of vertex associations for clusters to be merged

    bool            m_checkClusterSplitPositions;       ///< Whether to check the cluster split positions, if there are splits in multiple views
    float           m_vetoMergeXDifference;             ///< The x distance between split positions in two views below which may refuse a merge
    float           m_vetoMergeXOverlap;                ///< The x overlap between candidate cluster sliding fits below which may refuse a merge
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *SplitShowersTool::Factory::CreateAlgorithmTool() const
{
    return new SplitShowersTool();
}

} // namespace lar_content

#endif // #ifndef SPLIT_SHOWERS_TOOL_H
