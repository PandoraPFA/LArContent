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

    typedef std::map<const pandora::Cluster*, unsigned int> ClusterInfoMap;

    /**
     *  @brief  Get the list of all input pfos (possibly from a number of separate named lists)
     *
     *  @param  pfoList to receive the input pfo list
     */
    void GetInputPfoList(pandora::PfoList &pfoList) const;

    /**
     *  @brief  Store the number of calo hits per cluster in a cluster info map
     *
     *  @param  pCluster address of the relevant cluster
     *  @param  clusterInfoMap the cluster info map
     */
    void StoreNCaloHitsPerCluster(const pandora::Cluster *const pCluster, ClusterInfoMap &clusterInfoMap) const;

    /**
     *  @brief  Store the number of branches per cluster in a cluster info map
     *
     *  @param  pCluster address of the relevant cluster
     *  @param  branchList the branch list
     *  @param  clusterInfoMap the cluster info map
     */
    void StoreNBranchesPerCluster(const pandora::Cluster *const pCluster, const pandora::ClusterVector &branchList, ClusterInfoMap &clusterInfoMap) const;

    /**
     *  @brief  Remove pfos containing clusters to which many shower branches have been added
     *
     *  @param  pClusterList address of the cluster list
     *  @param  pfoList the pfo list
     *  @param  nCaloHitsPerCluster the cluster info map
     *  @param  nBranchesPerCluster the cluster info map
     */
    void RemoveShowerPfos(const pandora::ClusterList *const pClusterList, pandora::PfoList &pfoList, const ClusterInfoMap &nCaloHitsPerCluster,
        const ClusterInfoMap &nBranchesPerCluster) const;

    /**
     *  @brief  Get the address of the pfo containing a specified cluster
     *
     *  @param  pCluster address of the relevant cluster
     *  @param  pfoList the list of all input pfos
     *  @param  pTargetPfo to receive the address of the target pfo
     */
    void FindTargetPfo(pandora::Cluster *const pCluster, const pandora::PfoList &pfoList, pandora::Pfo *&pTargetPfo) const;

    std::string             m_inputClusterListName;     ///< The name of the input cluster list
    pandora::StringVector   m_inputPfoListNames;        ///< The names of the input pfo lists

    unsigned int            m_minCaloHitsPerCluster;    ///< The minimum number of calo hits per (seed or branch) cluster
    bool                    m_shouldRemoveShowerPfos;   ///< Whether to delete any existing pfos to which many shower branches have been added

    unsigned int            m_showerLikeNBranches;      ///< The minimum number of branches before cluster is declared shower like
    float                   m_showerLikeCaloHitRatio;   ///< The minimum ratio of final to original calo hits before cluster is declared shower like
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ClusterCharacterisationAlgorithm::Factory::CreateAlgorithm() const
{
    return new ClusterCharacterisationAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_CLUSTER_CHARACTERISATION_ALGORITHM_H
