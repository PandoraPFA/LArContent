/**
 *  @file   LArContent/include/LArThreeDReco/LArGenericTrackMatching/ClearTrackFragmentsTool.h
 *
 *  @brief  Header file for the clear track fragments tool class.
 *
 *  $Log: $
 */
#ifndef CLEAR_TRACK_FRAGMENTS_TOOL_H
#define CLEAR_TRACK_FRAGMENTS_TOOL_H 1

#include "LArThreeDReco/LArThreeDBase/ThreeDFragmentsBaseAlgorithm.h"

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

    bool Run(ThreeDFragmentsBaseAlgorithm *pAlgorithm, TensorType &overlapTensor);

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
    bool FindTrackFragments(ThreeDFragmentsBaseAlgorithm *pAlgorithm, const TensorType &overlapTensor) const;

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
     *  @brief  Process a tensor element, reclustering the fragments as required
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  overlapResult the overlap result
     *  @param  unavailableClusters to receive the list of clusters now unavailable for future particle reconstruction
     *  @param  newlyAvailableClusters to receive the list of clusters newly made available for future particle reconstruction
     *  @param  pFragmentCluster to receive the address of the new fragment cluster
     */
    void ProcessTensorElement(ThreeDFragmentsBaseAlgorithm *pAlgorithm, const TensorType::OverlapResult &overlapResult,
        pandora::ClusterList &unavailableClusters, pandora::ClusterList &newlyAvailableClusters, pandora::Cluster *&pFragmentCluster) const;

    /**
     *  @brief  Rearrange the hits in a cluster from the fragment list, using the Pandora fragmentation mechanism
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  pCluster address of the input cluster
     *  @param  daughterHits the full list of hits to place in the new fragment cluster
     *  @param  separateHits the full list of hits that are not to be placed in the new fragment cluster
     *  @param  newlyAvailableClusters to receive the list of clusters newly made available for future particle reconstruction
     *  @param  pFragmentCluster to receive the address of the new fragment cluster
     */
    void Recluster(ThreeDFragmentsBaseAlgorithm *pAlgorithm, pandora::Cluster *pCluster, const pandora::CaloHitList &daughterHits,
        const pandora::CaloHitList &separateHits, pandora::ClusterList &newlyAvailableClusters, pandora::Cluster *&pFragmentCluster) const;

    /**
     *  @brief  Update the tensor following the fragmentation operations performed by this tool
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  overlapTensor the overlap tensor
     *  @param  unavailableClusters the list of clusters now unavailable for future particle reconstruction
     *  @param  newlyAvailableClusters the list of clusters newly made available for future particle reconstruction
     */
    void UpdateTensor(ThreeDFragmentsBaseAlgorithm *pAlgorithm, const TensorType &overlapTensor, const pandora::ClusterList &unavailableClusters,
        const pandora::ClusterList &newlyAvailableClusters) const;

    /**
     *  @brief  Get a list of the tensor key clusters for which tensor elements have been impacted by fragmentation operations
     * 
     *  @param  overlapTensor the overlap tensor
     *  @param  unavailableClusters the list of clusters now unavailable for future particle reconstruction
     *  @param  affectedKeyClusters to receive the list of tensor key clusters that have been affected by fragmentation operations
     */
    void GetAffectedKeyClusters(const TensorType &overlapTensor, const pandora::ClusterList &unavailableClusters, pandora::ClusterList &affectedKeyClusters) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int        m_minMatchedSamplingPoints;         ///< The minimum number of matched sampling points
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
