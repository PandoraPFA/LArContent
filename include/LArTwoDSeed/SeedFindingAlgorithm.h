/**
 *  @file   LArContent/include/LArTwoDSeed/SeedFindingAlgorithm.h
 * 
 *  @brief  Header file for the seed finding algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_SEED_FINDING_ALGORITHM_H
#define LAR_SEED_FINDING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  SeedFindingAlgorithm class
 */
class SeedFindingAlgorithm : public pandora::Algorithm
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
     *  @brief  ParticleSeed class
     */
    class ParticleSeed
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  pCluster address of the cluster
         */
        ParticleSeed(pandora::Cluster *pCluster);

        /**
         *  @brief  Add a list of clusters to the particle seed
         * 
         *  @param  clusterList the cluster list
         */
        void AddClusterList(const pandora::ClusterList &clusterList);

        /**
         *  @brief  Get the particle seed cluster list
         * 
         *  @return the particle seed cluster list
         */
        const pandora::ClusterList &GetClusterList() const;

    private:
        pandora::ClusterList    m_clusterList;      ///< The cluster list
    };

    typedef std::vector<ParticleSeed*> ParticleSeedVector;

    /**
     *  @brief  Populate particle seed vector with candidate particle seeds
     * 
     *  @param  vertexSeedClusters the vertex seed cluster list
     *  @param  nonSeedClusters the non seed cluster list
     *  @param  particleSeedVector to receive the populated particle seed vector
     */
    void GetParticleSeeds(const pandora::ClusterList &vertexSeedClusters, const pandora::ClusterList &nonSeedClusters,
        ParticleSeedVector &particleSeedVector) const;

    /**
     *  @brief  Find associated particle seeds
     * 
     *  @param  pParticleSeed address of the particle seed
     *  @param  candidateSeeds the list of candidate associated seeds
     *  @param  associatedSeeds to receive the list of associated seeds
     */
    void FindAssociatedSeeds(ParticleSeed *const pParticleSeed, ParticleSeedVector &candidateSeeds, ParticleSeedVector &associatedSeeds) const;

    /**
     *  @brief  Determine whether two particle seeds are associated
     * 
     *  @param  pParticleSeedI address of particle seed I
     *  @param  pParticleSeedJ address of particle seed J
     * 
     *  @return boolean
     */
    bool AreParticleSeedsAssociated(const ParticleSeed *const pParticleSeedI, const ParticleSeed *const pParticleSeedJ) const;

     /**
     *  @brief  Determine whether two particle seed clusters are associated
     * 
     *  @param  pClusterI address of particle seed cluster I
     *  @param  pClusterJ address of particle seed cluster J
     * 
     *  @return boolean
     */
    bool AreSeedClustersAssociated(const pandora::Cluster *const pClusterI, const pandora::Cluster *const pClusterJ) const;

    /**
     *  @brief  Make cluster merges
     * 
     *  @param  particleSeedVector the particle seed list
     *  @param  finalClusterList to receive the final cluster list
     */
    void MakeClusterMerges(const ParticleSeedVector &particleSeedVector, pandora::ClusterList &finalClusterList) const;

    /**
     *  @brief  Sort particle seed by layer span of constituent clusters
     * 
     *  @param  pLhs address of first particle seed
     *  @param  pRhs address of second particle seed
     */
    static bool SortByLayerSpan(const ParticleSeed *const pLhs, const ParticleSeed *const pRhs);

    std::string         m_seedClusterListName;      ///< The seed cluster list name
    std::string         m_nonSeedClusterListName;   ///< The non seed cluster list name

    int                 m_initialLengthCut;         ///< 
    int                 m_finalLengthCut;           ///< 
    int                 m_initialChangeIter;        ///< 
    int                 m_finalChangeIter;          ///< 
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *SeedFindingAlgorithm::Factory::CreateAlgorithm() const
{
    return new SeedFindingAlgorithm();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline SeedFindingAlgorithm::ParticleSeed::ParticleSeed(pandora::Cluster *pCluster)
{
    if (!m_clusterList.insert(pCluster).second)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_ALREADY_PRESENT);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SeedFindingAlgorithm::ParticleSeed::AddClusterList(const pandora::ClusterList &clusterList)
{
    for (pandora::ClusterList::const_iterator iter = clusterList.begin(), iterEnd = clusterList.end(); iter != iterEnd; ++iter)
    {
        if (!m_clusterList.insert(*iter).second)
            throw pandora::StatusCodeException(pandora::STATUS_CODE_ALREADY_PRESENT);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::ClusterList &SeedFindingAlgorithm::ParticleSeed::GetClusterList() const
{
    return m_clusterList;
}

} // namespace lar

#endif // #ifndef LAR_SEED_FINDING_ALGORITHM_H
