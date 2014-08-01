/**
 *  @file   LArContent/include/LArThreeDReco/LArTrackFragments/ClearTrackFragmentsTool.h
 *
 *  @brief  Header file for the clear track fragments tool class.
 *
 *  $Log: $
 */
#ifndef CLEAR_TRACK_FRAGMENTS_TOOL_H
#define CLEAR_TRACK_FRAGMENTS_TOOL_H 1

#include "LArThreeDReco/LArTrackFragments/ThreeDTrackFragmentsAlgorithm.h"

namespace lar
{

/**
 *  @brief  ClearTrackFragmentsTool class
 */
class ClearTrackFragmentsTool : public FragmentTensorTool
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

    bool Run(ThreeDTrackFragmentsAlgorithm *pAlgorithm, TensorType &overlapTensor);

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
    bool FindTrackFragments(ThreeDTrackFragmentsAlgorithm *pAlgorithm, const TensorType &overlapTensor) const;

    /**
     *  @brief  Get the list of elements connected to a given cluster and check its suitability (no ambiguities)
     *
     *  @param  overlapTensor the overlap tensor
     *  @param  pCluster address of the key cluster
     *  @param  elementList to receive the element list
     *
     *  @return boolean
     */
    bool GetAndCheckElementList(const TensorType &overlapTensor, pandora::Cluster *pCluster, TensorType::ElementList &elementList) const;

    /**
     *  @brief  Check the list of hits, stored in tensor elements, for ambiguities
     *
     *  @param  elementList the element list
     *
     *  @return boolean
     */
    bool CheckForHitAmbiguities(const TensorType::ElementList &elementList) const;

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
     *  @param  overlapResult the overlap result
     *  @param  modifiedClusters to receive the list of modified clusters
     *  @param  deletedClusters to receive the list of deleted clusters
     *  @param  pFragmentCluster to receive the address of the new fragment cluster
     */
    void ProcessTensorElement(ThreeDTrackFragmentsAlgorithm *pAlgorithm, const TensorType::OverlapResult &overlapResult,
        pandora::ClusterList &modifiedClusters, pandora::ClusterList &deletedClusters, pandora::Cluster *&pFragmentCluster) const;

    /**
     *  @brief  Rearrange the hits in a cluster from the fragment list, using the Pandora fragmentation mechanism
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  pCluster address of the input cluster
     *  @param  daughterHits the full list of hits to place in the new fragment cluster
     *  @param  separateHits the full list of hits that are not to be placed in the new fragment cluster
     *  @param  deletedClusters to receive the list of deleted clusters
     *  @param  pFragmentCluster to receive the address of the new fragment cluster
     */
    void Recluster(ThreeDTrackFragmentsAlgorithm *pAlgorithm, pandora::Cluster *pCluster, const pandora::CaloHitList &daughterHits,
        const pandora::CaloHitList &separateHits, pandora::ClusterList &deletedClusters, pandora::Cluster *&pFragmentCluster) const;

    /**
     *  @brief  Rebuild clusters after fragmentation
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  modifiedClusters the list of clusters to rebuild
     *  @param  deletedClusters the list of vetoed clusters
     *  @param  newClusters the list of new clusters
     */
    void RebuildClusters(ThreeDTrackFragmentsAlgorithm *pAlgorithm, const pandora::ClusterList &modifiedClusters,
        const pandora::ClusterList &deletedClusters, pandora::ClusterList &newClusters) const;

    /**
     *  @brief  Update the tensor following the fragmentation operations performed by this tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  overlapTensor the overlap tensor
     *  @param  unavailableClusters the list of clusters now unavailable for future particle reconstruction
     *  @param  newAvailableClusters the list of clusters newly made available for future particle reconstruction
     */
    void UpdateTensor(ThreeDTrackFragmentsAlgorithm *pAlgorithm, const TensorType &overlapTensor,
        const pandora::ClusterList &unavailableClusters, const pandora::ClusterList &newAvailableClusters) const;

    /**
     *  @brief  Get a list of the tensor key clusters for which tensor elements have been impacted by fragmentation operations
     *
     *  @param  overlapTensor the overlap tensor
     *  @param  unavailableClusters the list of clusters now unavailable for future particle reconstruction
     *  @param  affectedKeyClusters to receive the list of tensor key clusters that have been affected by fragmentation operations
     */
    void GetAffectedKeyClusters(const TensorType &overlapTensor, const pandora::ClusterList &unavailableClusters,
        pandora::ClusterList &affectedKeyClusters) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float               m_minMatchedSamplingPointFraction;  ///< The minimum fraction of matched sampling points
    unsigned int        m_minMatchedHits;                   ///< The minimum number of matched calo hits
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *ClearTrackFragmentsTool::Factory::CreateAlgorithmTool() const
{
    return new ClearTrackFragmentsTool();
}

} // namespace lar

#endif // #ifndef CLEAR_TRACK_FRAGMENTS_TOOL_H
