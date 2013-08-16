/**
 *  @file   LArContent/include/LArReclustering/ShowerRebuildingAlgorithm.h
 * 
 *  @brief  Header file for the shower rebuilding algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_SHOWER_REBUILDING_ALGORITHM_H
#define LAR_SHOWER_REBUILDING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  ShowerRebuildingAlgorithm class
 */
class ShowerRebuildingAlgorithm : public pandora::Algorithm
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

    /*
     *  @brief  RebuildThreeDShowers
     */
    void RebuildThreeDShowers() const;

    /*
     *  @brief  RebuildThreeDShower
     * 
     *  @param  pPfo
     *  @param  pSeedCluster
     *  @param  pfoVector
     *  @param  seedMerges
     *  @param  nonSeedMerges
     */
    void RebuildThreeDShower(const pandora::ParticleFlowObject *const pPfo, pandora::Cluster *pSeedCluster, pandora::PfoVector &pfoVector,
        pandora::ClusterList &seedMerges, pandora::ClusterList &nonSeedMerges) const;

    /*
     *  @brief  GetSeedClusters
     * 
     *  @param  pPfo
     *  @param  pSeedCluster
     *  @param  pfoVector
     *  @param  seedMerges
     *  @param  nonSeedMerges
     */
    void GetSeedClusters(const pandora::ParticleFlowObject *const pPfo, pandora::Cluster *&pSeedClusterU, pandora::Cluster *&pSeedClusterV,
        pandora::Cluster *&pSeedClusterW) const;

    /*
     *  @brief  GetAssociatedClusters
     * 
     *  @param  pCluster
     *  @param  seedMerges
     *  @param  nonSeedMerges
     */
    void GetAssociatedClusters(const pandora::Cluster *const pCluster, pandora::ClusterList &seedMerges, pandora::ClusterList &nonSeedMerges) const;

    /*
     *  @brief  PerformClusterMerges
     * 
     *  @param  pCluster
     *  @param  clusterList
     *  @param  clusterListName1
     *  @param  clusterListName2
     */
    void PerformClusterMerges(pandora::Cluster *pCluster, const pandora::ClusterList &clusterList, const std::string &clusterListName1,
        const std::string &clusterListName2) const;

    /*
     *  @brief  RebuildTwoDShowers
     * 
     *  @param  seedClusterListName
     *  @param  nonSeedClusterListName
     */
    void RebuildTwoDShowers(const std::string &seedClusterListName, const std::string &nonSeedClusterListName) const;

    /*
     *  @brief  RestoreLoneShowers
     */
    void RestoreLoneShowers() const;

    /*
     *  @brief  Sort pfos by number of hits
     * 
     *  @param  pLhs address of first pfo
     *  @param  pRhs address of second pfo
     */
    static bool SortPfosByNHits(const pandora::ParticleFlowObject *const pLhs, const pandora::ParticleFlowObject *const pRhs);

    /*
     *  @brief  Sort pfos by number of hits
     * 
     *  @param  pLhs address of first cluster
     *  @param  pRhs address of second cluster
     */
    static bool SortClustersByNHits(const pandora::Cluster *const pLhs, const pandora::Cluster *const pRhs);

    std::string     m_particleListName;             ///< 

    std::string     m_seedClusterListNameU;         ///< 
    std::string     m_seedClusterListNameV;         ///< 
    std::string     m_seedClusterListNameW;         ///< 

    std::string     m_nonSeedClusterListNameU;      ///< 
    std::string     m_nonSeedClusterListNameV;      ///< 
    std::string     m_nonSeedClusterListNameW;      ///< 
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ShowerRebuildingAlgorithm::Factory::CreateAlgorithm() const
{
    return new ShowerRebuildingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_SHOWER_REBUILDING_ALGORITHM_H
