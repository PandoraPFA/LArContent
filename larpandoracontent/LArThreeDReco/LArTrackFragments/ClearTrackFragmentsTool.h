/**
 *  @file   larpandoracontent/LArThreeDReco/LArTrackFragments/ClearTrackFragmentsTool.h
 *
 *  @brief  Header file for the clear track fragments tool class.
 *
 *  $Log: $
 */
#ifndef CLEAR_TRACK_FRAGMENTS_TOOL_H
#define CLEAR_TRACK_FRAGMENTS_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArTrackFragments/ThreeViewTrackFragmentsAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  ClearTrackFragmentsTool class
 */
class ClearTrackFragmentsTool : public FragmentTensorTool
{
public:
    /**
     *  @brief  Default constructor
     */
    ClearTrackFragmentsTool();

    bool Run(ThreeViewTrackFragmentsAlgorithm *const pAlgorithm, TensorType &overlapTensor);

private:
    /**
     *  @brief  Find suitable matching track fragments in the overlap tensor to use for 3D particle creation,
     *          return value indicates whether particles are made
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  overlapTensor the overlap tensor
     *
     *  @return boolean
     */
    bool FindTrackFragments(ThreeViewTrackFragmentsAlgorithm *const pAlgorithm, const TensorType &overlapTensor) const;

    /**
     *  @brief  Get the list of elements connected to a given cluster and check its suitability (no ambiguities)
     *
     *  @param  overlapTensor the overlap tensor
     *  @param  pCluster address of the key cluster
     *  @param  elementList to receive the element list
     *
     *  @return boolean
     */
    bool GetAndCheckElementList(const TensorType &overlapTensor, const pandora::Cluster *const pCluster, TensorType::ElementList &elementList) const;

    /**
     *  @brief  Check whether the overlap result passes matched sampling point and number of matched hit checks
     *
     *  @param  overlapResult the overlap result
     *
     *  @return boolean
     */
    bool CheckOverlapResult(const TensorType::OverlapResult &overlapResult) const;

    /**
     *  @brief  Select a list of clear track-like elements from a set of connected tensor elements
     *
     *  @param  elementList the full list of connected tensor elements
     *  @param  iteratorList to receive a list of iterators to long track-like elements
     */
    void SelectClearElements(const TensorType::ElementList &elementList, IteratorList &iteratorList) const;

    /**
     *  @brief  Process a tensor element, reclustering the fragments as required
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  overlapTensor the overlap tensor
     *  @param  overlapResult the overlap result
     *  @param  pFragmentCluster to receive the address of the new fragment cluster
     */
    void ProcessTensorElement(ThreeViewTrackFragmentsAlgorithm *const pAlgorithm, const TensorType &overlapTensor,
        const TensorType::OverlapResult &overlapResult, const pandora::Cluster *&pFragmentCluster) const;

    /**
     *  @brief  Rearrange the hits in a cluster from the fragment list, using the Pandora fragmentation mechanism
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  pCluster address of the input cluster
     *  @param  daughterHits the full list of hits to place in the new fragment cluster
     *  @param  separateHits the full list of hits that are not to be placed in the new fragment cluster
     *  @param  deletedClusters to receive the set of deleted clusters
     *  @param  badClusters the set of clusters that should not be dereferenced
     *  @param  pFragmentCluster to receive the address of the new fragment cluster
     */
    void Recluster(ThreeViewTrackFragmentsAlgorithm *const pAlgorithm, const pandora::Cluster *const pCluster,
        const pandora::CaloHitList &daughterHits, const pandora::CaloHitList &separateHits, pandora::ClusterSet &deletedClusters,
        pandora::ClusterSet &badClusters, const pandora::Cluster *&pFragmentCluster) const;

    /**
     *  @brief  Rebuild clusters after fragmentation
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  modifiedClusters the list of clusters to rebuild
     *  @param  newClusters the list of new clusters
     */
    void RebuildClusters(ThreeViewTrackFragmentsAlgorithm *const pAlgorithm, const pandora::ClusterList &modifiedClusters,
        pandora::ClusterList &newClusters) const;

    /**
     *  @brief  Get a list of the tensor key clusters for which tensor elements have been impacted by fragmentation operations
     *
     *  @param  overlapTensor the overlap tensor
     *  @param  clustersToRemoveFromTensor the list of clusters removed from the tensor by fragmentation operations
     *  @param  affectedKeyClusters to receive the list of tensor key clusters that have been affected by fragmentation operations
     */
    void GetAffectedKeyClusters(const TensorType &overlapTensor, const pandora::ClusterList &clustersToRemoveFromTensor,
        pandora::ClusterList &affectedKeyClusters) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_minMatchedSamplingPointFraction; ///< The minimum fraction of matched sampling points
    unsigned int m_minMatchedHits;           ///< The minimum number of matched calo hits
};

} // namespace lar_content

#endif // #ifndef CLEAR_TRACK_FRAGMENTS_TOOL_H
