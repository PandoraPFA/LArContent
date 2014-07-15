/**
 *  @file   LArContent/include/LArTwoDReco/LArSeedFinding/SeedGrowingAlgorithm.h
 * 
 *  @brief  Header file for the seed growing algorithm base class.
 * 
 *  $Log: $
 */
#ifndef LAR_SEED_GROWING_ALGORITHM_H
#define LAR_SEED_GROWING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  SeedGrowingAlgorithm class
 */
class SeedGrowingAlgorithm : public pandora::Algorithm
{
protected:
    /**
     *  @brief  AssociationType enum
     */
    enum AssociationType
    {
        NONE = 0,
        SINGLE_ORDER = 1,
        STANDARD = 2,
        STRONG = 3
    };

    /**
     *  @brief  Association class
     */
    class Association
    {
    public:
        /**
         *  @brief  Default constructor
         */
        Association();

        /**
         *  @brief  Constructor
         * 
         *  @param  order the association order
         *  @param  type the association type
         */
        Association(const unsigned int order, const AssociationType type);

        /**
         *  @brief  Set association order
         * 
         *  @param  order the association order
         */
        void SetOrder(const unsigned int order);

        /**
         *  @brief  Set association type
         * 
         *  @param  associationType the association type
         */
        void SetType(const AssociationType associationType);

        /**
         *  @brief  Get association order
         * 
         *  @return the association order
         */
        unsigned int GetOrder() const;

        /**
         *  @brief  Get association type
         * 
         *  @return the association type
         */
        AssociationType GetType() const;

    private:
        unsigned int    m_order;
        AssociationType m_type;
    };

    typedef std::map<pandora::Cluster*, Association> ClusterAssociationMap;
    typedef std::map<pandora::Cluster*, ClusterAssociationMap> ClusterUsageMap;

    /**
     *  @brief  Determine whether two clusters are associated
     * 
     *  @param  pClusterSeed address of cluster seed (may be daughter of primary seed)
     *  @param  pCluster address of cluster
     * 
     *  @return the association type
     */
    virtual AssociationType AreClustersAssociated(const pandora::Cluster *const pClusterSeed, const pandora::Cluster *const pCluster) const = 0;

    /**
     *  @brief  Find clusters associated with a particle seed
     * 
     *  @param  pParticleSeed address of the particle seed
     *  @param  candidateClusters list of clusters which may be associated with seed
     *  @param  forwardUsageMap the particle seed usage map
     *  @param  backwardUsageMap the cluster usage map
     */
    void FindAssociatedClusters(pandora::Cluster *const pParticleSeed, pandora::ClusterVector &candidateClusters,
        ClusterUsageMap &forwardUsageMap, ClusterUsageMap &backwardUsageMap) const;

    typedef std::map<pandora::Cluster*, pandora::ClusterVector> SeedAssociationList;

    /**
     *  @brief  Identify cluster merges
     * 
     *  @param  particleSeedVector the list of all particle seeds
     *  @param  backwardUsageMap the map from cluster to particle seed associations
     *  @param  seedAssociationList to receive the populated seed association list
     */
    void IdentifyClusterMerges(const pandora::ClusterVector &particleSeedVector, const ClusterUsageMap &backwardUsageMap,
        SeedAssociationList &seedAssociationList) const;

    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline SeedGrowingAlgorithm::Association::Association() :
    m_order(std::numeric_limits<unsigned int>::max()),
    m_type(NONE)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline SeedGrowingAlgorithm::Association::Association(const unsigned int order, const AssociationType type) :
    m_order(order),
    m_type(type)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void SeedGrowingAlgorithm::Association::SetOrder(const unsigned int order)
{
    m_order = order;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void SeedGrowingAlgorithm::Association::SetType(const AssociationType associationType)
{
    m_type = associationType;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int SeedGrowingAlgorithm::Association::GetOrder() const
{
    return m_order;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline SeedGrowingAlgorithm::AssociationType SeedGrowingAlgorithm::Association::GetType() const
{
    return m_type;
}

} // namespace lar

#endif // #ifndef LAR_SEED_GROWING_ALGORITHM_H
