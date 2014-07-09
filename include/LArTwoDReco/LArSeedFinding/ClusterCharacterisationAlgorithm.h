/**
 *  @file   LArContent/include/LArTwoDReco/LArSeedFinding/ClusterCharacterisationAlgorithm.h
 * 
 *  @brief  Header file for the cluster characterisation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_CLUSTER_CHARACTERISATION_ALGORITHM_H
#define LAR_CLUSTER_CHARACTERISATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArTwoDReco/LArSeedFinding/SeedBranchGrowingAlgorithm.h"

namespace lar
{

/**
 *  @brief  ClusterCharacterisationAlgorithm class
 */
class ClusterCharacterisationAlgorithm : public SeedBranchGrowingAlgorithm
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
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Get the next seed candidate, using a list of available candidates and a list of those already used
     * 
     *  @param  pClusterList the list of available seed candidates
     *  @param  usedClusters the list of candidates already considered
     *  @param  pSeedCluster to receive the address of the next seed candidate
     * 
     *  @return whether a seed candidate has been found
     */
    bool GetNextSeedCandidate(const pandora::ClusterList *const pClusterList, const pandora::ClusterList &usedClusters,
        pandora::Cluster *&pSeedCluster) const;

    /**
     *  @brief  Get the seed association list for a given vector of particle seed candidates
     * 
     *  @param  particleSeedVector the particle seed vector
     *  @param  pClusterList the address of the input cluster list
     *  @param  seedAssociationList to receive the populated seed association list
     */
    void GetSeedAssociationList(const pandora::ClusterVector &particleSeedVector, const pandora::ClusterList *const pClusterList,
        SeedAssociationList &seedAssociationList) const;

    /**
     *  @brief  Check a provided seed association list for consistency, making changes as required
     * 
     *  @param  seedIter iterator to an element in the input seed association list
     *  @param  finalSeedAssociationList to receive the output seed association list
     */
    void CheckSeedAssociationList(SeedAssociationList::const_iterator seedIter, SeedAssociationList &finalSeedAssociationList) const;

    /**
     *  @brief  Get a figure of merit representing the consistency of the provided seed associated list
     * 
     *  @param  seedAssociationList the seed association list
     * 
     *  @return the figure of merit
     */
    float GetFigureOfMerit(const SeedAssociationList &seedAssociationList) const;

    /**
     *  @brief  Custom sorting for clusters to determine order in which seeds are considered
     *
     *  @param  pLhs address of first cluster
     *  @param  pRhs address of second cluster
     */
    static bool SortClusters(const pandora::Cluster *const pLhs, const pandora::Cluster *const pRhs);

    std::string     m_inputClusterListName;     ///< The name of the input cluster list
    unsigned int    m_minCaloHitsPerCluster;    ///< The minimum number of calo hits per (seed or branch) cluster
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ClusterCharacterisationAlgorithm::Factory::CreateAlgorithm() const
{
    return new ClusterCharacterisationAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_CLUSTER_CHARACTERISATION_ALGORITHM_H
