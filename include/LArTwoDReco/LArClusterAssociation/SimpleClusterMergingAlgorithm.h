/**
 *  @file   LArContent/include/LArTwoDReco/LArClusterAssociation/SimpleClusterMergingAlgorithm.h
 *
 *  @brief  Header file for the simple cluster merging algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_SIMPLE_CLUSTER_MERGING_ALGORITHM_H
#define LAR_SIMPLE_CLUSTER_MERGING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArTwoDReco/LArClusterAssociation/ClusterMergingAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  SimpleClusterMergingAlgorithm class
 */
class SimpleClusterMergingAlgorithm : public ClusterMergingAlgorithm
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

private:
    void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;
    void PopulateClusterMergeMap(const pandora::ClusterVector &clusterVector, ClusterMergeMap &clusterMergeMap) const;

    /**
     *  @brief Decide whether two clusters are associated
     *
     *  @param pClusterI the address of the first cluster
     *  @param pClusterJ the address of the second cluster
     *
     *  @return boolean
     */
    bool IsAssociated(const pandora::Cluster* const pClusterI, const pandora::Cluster* const pClusterJ) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int m_minCaloHitsPerCluster;    ///< The min number of calo hits per candidate cluster
    float        m_maxClusterSeparation;     ///< Maximum distance at which clusters can be joined
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *SimpleClusterMergingAlgorithm::Factory::CreateAlgorithm() const
{
    return new SimpleClusterMergingAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_SIMPLE_CLUSTER_MERGING_ALGORITHM_H
