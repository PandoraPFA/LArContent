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

namespace lar
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
     *  @param  clusterMergeMap to be populated with cluster merges
     */
    void FindShowerMerges(ThreeDShowersAlgorithm *pAlgorithm, const IteratorList &iteratorList, ClusterMergeMap &clusterMergeMap) const;

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
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *SplitShowersTool::Factory::CreateAlgorithmTool() const
{
    return new SplitShowersTool();
}

} // namespace lar

#endif // #ifndef SPLIT_SHOWERS_TOOL_H
