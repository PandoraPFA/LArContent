/**
 *  @file   LArContent/include/LArTwoDReco/LArCosmicRay/CosmicRayShowerMergingAlgorithm.h
 *
 *  @brief  Header file for the cosmic ray shower merging algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_COSMIC_RAY_SHOWER_MERGING_ALGORITHM_H
#define LAR_COSMIC_RAY_SHOWER_MERGING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArTwoDReco/LArClusterAssociation/ClusterMergingAlgorithm.h"

namespace lar
{

/**
 *  @brief  CosmicRayShowerMergingAlgorithm class
 */
class CosmicRayShowerMergingAlgorithm : public ClusterMergingAlgorithm
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
     *  @brief Separate clusters into cosmic rays, shower seeds and shower fragments
     *
     *  @param inputVector the input cluster vector
     *  @param seedVector the output cluster vector of shower seeds
     *  @param nonSeedVector the output cluster vector of shower fragments
     */
    void SeparateClusters(const pandora::ClusterVector &inputVector, pandora::ClusterVector &seedVector,
        pandora::ClusterVector &nonSeedVector) const;

    /**
     *  @brief Make associations between the clusters in a single cluster vector
     *
     *  @param clusterVector the input cluster vector
     *  @param outputMergemap the map of cluster associations
     */
    void FillClusterMergeMap(const pandora::ClusterVector &clusterVector, ClusterMergeMap &outputMergeMap) const;

    /**
     *  @brief Make associations between the clusters in two cluster vectors
     *
     *  @param firstVector the first cluster vector
     *  @param secondVector the second cluster vector
     *  @param outputMergeMap the map of cluster associations
     */
    void FillClusterMergeMap(const pandora::ClusterVector &firstVector, const  pandora::ClusterVector &secondVector,
        ClusterMergeMap &outputMergeMap) const;

    /**
     *  @brief Collect up groups of associated clusters
     *
     *  @param seedVector the vector of seed clusters, to begin forming groups
     *  @param vetoVector the vector of veto clusters, cannot be in any of the groups
     *  @param inputMergeMap the input map of cluster associations
     *  @param outputMergeMap the output map of grouped clusters
     */
    void FillClusterMergeMap(const pandora::ClusterVector &seedVector, const pandora::ClusterVector &vetoVector,
                             const ClusterMergeMap &inputMergeMap, ClusterMergeMap &outputMergeMap) const;

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


    float m_minClusterLength;         ///<
    float m_maxClusterLength;         ///<

    float m_maxSeedDisplacement;      ///<
    float m_maxNonSeedDisplacement;   ///<
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CosmicRayShowerMergingAlgorithm::Factory::CreateAlgorithm() const
{
    return new CosmicRayShowerMergingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_COSMIC_RAY_SHOWER_MERGING_ALGORITHM_H
