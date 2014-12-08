/**
 *  @file   LArContent/include/LArTwoDReco/LArClusterMopUp/ProximityBasedMergingAlgorithm.h
 * 
 *  @brief  Header file for the proximity based cluster merging algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_PROXIMITY_BASED_MERGING_ALGORITHM_H
#define LAR_PROXIMITY_BASED_MERGING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArTwoDReco/LArClusterMopUp/ClusterMopUpAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  ProximityBasedMergingAlgorithm class
 */

class ProximityBasedMergingAlgorithm : public ClusterMopUpAlgorithm
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
    ProximityBasedMergingAlgorithm();

private:
    void ClusterMopUp(const pandora::ClusterList &pfoClusters, const pandora::ClusterList &remnantClusters, const ClusterToListNameMap &clusterToListNameMap) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ProximityBasedMergingAlgorithm::Factory::CreateAlgorithm() const
{
    return new ProximityBasedMergingAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_PROXIMITY_BASED_MERGING_ALGORITHM_H
