/**
 *  @file   LArContent/include/LArThreeDReco/LArTrackMatching/ThreeDKinkBaseTool.h
 * 
 *  @brief  Header file for the three d kink base tool
 * 
 *  $Log: $
 */
#ifndef THREE_D_KINK_BASE_TOOL_H
#define THREE_D_KINK_BASE_TOOL_H 1

#include "LArThreeDReco/LArTrackMatching/ThreeDTransverseTracksAlgorithm.h"

namespace lar
{

/**
 *  @brief  ThreeDKinkBaseTool class
 */
class ThreeDKinkBaseTool : public TensorManipulationTool
{
public:
    /**
     *  @brief  Constructor
     * 
     *  @param  nCommonClusters the number of common clusters to select
     */
    ThreeDKinkBaseTool(const unsigned int nCommonClusters);

protected:
    typedef std::map<pandora::Cluster*, pandora::CartesianPointList> SplitPositionMap;
    typedef std::map<pandora::Cluster*, pandora::ClusterList> ClusterMergeMap;
    
    /**
     *  @brief  Modification class
     */
    class Modification
    {
    public:
        SplitPositionMap            m_splitPositionMap;     ///< The split position map
        ClusterMergeMap             m_clusterMergeMap;      ///< The cluster merge map
        pandora::ClusterList        m_affectedClusters;     ///< The list of affected clusters
    };

    typedef std::vector<Modification> ModificationList;

    /**
     *  @brief  Whether a provided (iterator to a) tensor element passes the selection cuts for overshoot identification
     * 
     *  @param  eIter the iterator to the tensor element
     *  @param  usedClusters the list of used clusters
     */
    virtual bool PassesElementCuts(TensorType::ElementList::const_iterator eIter, const pandora::ClusterList &usedClusters) const;

    typedef std::vector<TensorType::ElementList::const_iterator> IteratorList;

    /**
     *  @brief  Get modification objects for a specific elements of the tensor, identifying required splits and merges for clusters
     * 
     *  @param  iteratorList list of iterators to relevant tensor elements
     *  @param  modificationList to be populated with modifications
     */
    virtual void GetIteratorListModifications(const IteratorList &iteratorList, ModificationList &modificationList) const = 0;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int    m_nCommonClusters;                  ///< The number of common clusters
    float           m_minMatchedFraction;               ///< The min matched sampling point fraction for use as a key tensor element
    unsigned int    m_minMatchedSamplingPoints;         ///< The min number of matched sampling points for use as a key tensor element
    float           m_minLongitudinalImpactParameter;   ///< The min longitudinal impact parameter for connecting accompanying clusters

private:
    bool Run(ThreeDTransverseTracksAlgorithm *pAlgorithm, TensorType &overlapTensor);

    /**
     *  @brief  Get modification objects, identifying required splits and merges for clusters
     * 
     *  @param  overlapTensor the overlap tensor
     *  @param  modificationList to be populated with modifications
     */
    void GetModifications(const TensorType &overlapTensor, ModificationList &modificationList) const;

    /**
     *  @brief  Select elements representing possible components of interest due to overshoots or undershoots in clustering
     * 
     *  @param  eIter iterator to a candidate element
     *  @param  elementList the provided element list
     *  @param  usedClusters the list of used clusters
     *  @param  iteratorList to receive a list of iterators to relevant elements
     */
    void SelectTensorElements(TensorType::ElementList::const_iterator eIter, const TensorType::ElementList &elementList,
        const pandora::ClusterList &usedClusters, IteratorList &iteratorList) const;

    /**
     *  @brief  Apply the changes cached in a modification list and update the tensor accordingly
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  modificationList the modification list
     * 
     *  @return whether changes to the tensor have been made
     */
    bool ApplyChanges(ThreeDTransverseTracksAlgorithm *pAlgorithm, const ModificationList &modificationList) const;

    /**
     *  @brief  Merge clusters together
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  clusterMergeMap the cluster merge map
     * 
     *  @return whether changes to the tensor have been made
     */
    bool MakeClusterMerges(ThreeDTransverseTracksAlgorithm *pAlgorithm, const ClusterMergeMap &clusterMergeMap) const;

    /**
     *  @brief  Make cluster splits
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  splitPositionMap the split position map
     * 
     *  @return whether changes to the tensor have been made
     */
    bool MakeClusterSplits(ThreeDTransverseTracksAlgorithm *pAlgorithm, const SplitPositionMap &splitPositionMap) const;

    /**
     *  @brief  Make a cluster split
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  splitPosition the split position
     *  @param  pCurrentCluster the cluster to split
     *  @param  pLowXCluster to receive the low x cluster
     *  @param  pHighXCluster to receive the high x cluster
     */
    void MakeClusterSplit(ThreeDTransverseTracksAlgorithm *pAlgorithm, const pandora::CartesianVector &splitPosition,
        pandora::Cluster *&pCurrentCluster, pandora::Cluster *&pLowXCluster, pandora::Cluster *&pHighXCluster) const;

    /**
     *  @brief  Sort split position cartesian vectors by increasing x coordinate
     * 
     *  @param  lhs the first cartesian vector
     *  @param  rhs the second cartesian vector
     */
    static bool SortSplitPositions(const pandora::CartesianVector &lhs, const pandora::CartesianVector &rhs);
};

} // namespace lar

#endif // #ifndef THREE_D_KINK_BASE_TOOL_H
