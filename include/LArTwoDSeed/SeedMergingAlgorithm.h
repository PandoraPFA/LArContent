/**
 *  @file   LArContent/include/LArTwoDSeed/SeedMergingAlgorithm.h
 * 
 *  @brief  Header file for the seed merging algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_SEED_MERGING_ALGORITHM_H
#define LAR_SEED_MERGING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  SeedMergingAlgorithm class
 */
class SeedMergingAlgorithm : public pandora::Algorithm
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
     *  @param  pClusterList the address of the candidate cluster list
     *  @param  particleSeedVector to receive the populated particle seed vector
     */
    void GetParticleSeeds(const pandora::ClusterList *const pClusterList, ParticleSeedVector &particleSeedVector) const;

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
    bool AreSeedClustersAssociated(pandora::Cluster *const pClusterI, pandora::Cluster *const pClusterJ) const;

    /**
     *  @brief  Make cluster merges
     * 
     *  @param  particleSeedVector the particle seed list
     */
    void MakeClusterMerges(const ParticleSeedVector &particleSeedVector) const;

    /**
     *  @brief  Sort particle seed by layer span of constituent clusters
     * 
     *  @param  pLhs address of first particle seed
     *  @param  pRhs address of second particle seed
     */
    static bool SortByLayerSpan(const ParticleSeed *const pLhs, const ParticleSeed *const pRhs);
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *SeedMergingAlgorithm::Factory::CreateAlgorithm() const
{
    return new SeedMergingAlgorithm();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline SeedMergingAlgorithm::ParticleSeed::ParticleSeed(pandora::Cluster *pCluster)
{
    if (!m_clusterList.insert(pCluster).second)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_ALREADY_PRESENT);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void SeedMergingAlgorithm::ParticleSeed::AddClusterList(const pandora::ClusterList &clusterList)
{
    for (pandora::ClusterList::const_iterator iter = clusterList.begin(), iterEnd = clusterList.end(); iter != iterEnd; ++iter)
    {
        if (!m_clusterList.insert(*iter).second)
            throw pandora::StatusCodeException(pandora::STATUS_CODE_ALREADY_PRESENT);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::ClusterList &SeedMergingAlgorithm::ParticleSeed::GetClusterList() const
{
    return m_clusterList;
}

} // namespace lar

#endif // #ifndef LAR_SEED_MERGING_ALGORITHM_H
