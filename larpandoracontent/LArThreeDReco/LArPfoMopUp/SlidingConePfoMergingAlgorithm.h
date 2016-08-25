/**
 *  @file   larpandoracontent/LArThreeDReco/LArPfoMopUp/SlidingConePfoMergingAlgorithm.h
 * 
 *  @brief  Header file for the sliding cone pfo merging algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_SLIDING_CONE_PFO_MERGING_ALGORITHM_H
#define LAR_SLIDING_CONE_PFO_MERGING_ALGORITHM_H 1

#include "larpandoracontent/LArThreeDReco/LArPfoMopUp/PfoMergingBaseAlgorithm.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  SlidingConePfoMergingAlgorithm class
 */
class SlidingConePfoMergingAlgorithm : public PfoMergingBaseAlgorithm
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

    /**
     *  @brief  Default constructor
     */
    SlidingConePfoMergingAlgorithm();

protected:
    pandora::StatusCode Run();

    typedef std::unordered_map<const pandora::Cluster*, const pandora::ParticleFlowObject*> ClusterToPfoMap;

    /**
     *  @brief  Get all 3d clusters contained in the input pfo lists, extract 3d track and shower clusters and a mapping from cluster to pfo
     * 
     *  @param  trackClusters3D to receive the list of 3d track clusters
     *  @param  showerClusters3D to receive the list of 3d shower clusters
     *  @param  clusterToPfoMap to receive the mapping from 3d cluster to pfo
     */
    void GetThreeDClusters(pandora::ClusterVector &trackClusters3D, pandora::ClusterVector &showerClusters3D, ClusterToPfoMap &clusterToPfoMap) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector       m_inputPfoListNames;            ///< The input pfo list names
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *SlidingConePfoMergingAlgorithm::Factory::CreateAlgorithm() const
{
    return new SlidingConePfoMergingAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_SLIDING_CONE_PFO_MERGING_ALGORITHM_H
