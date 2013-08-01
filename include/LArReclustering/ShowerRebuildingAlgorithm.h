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
     *  @param  pfoVector
     *  @param  pSeedCluster
     *  @param  pSeedCluster2
     *  @param  pSeedCluster3
     *  @param  seedClusterListName
     *  @param  nonSeedClusterListName
     */
    void RebuildThreeDShower(pandora::PfoVector &pfoVector, pandora::Cluster *pSeedCluster, pandora::Cluster *pSeedCluster2, pandora::Cluster *pSeedCluster3,
        const std::string &seedClusterListName, const std::string &nonSeedClusterListName) const;

    /*
     *  @brief  GetSeedClusters
     * 
     *  @param  pPfo
     *  @param  pClusterU
     *  @param  pClusterV
     *  @param  pClusterW
     */
    void GetSeedClusters(const pandora::ParticleFlowObject *const pPfo, pandora::Cluster *&pClusterU, pandora::Cluster *&pClusterV, pandora::Cluster *&pClusterW) const;

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
     *  @brief  Sort pfos by descending energy 
     * 
     *  @param  pLhs address of first pfo
     *  @param  pRhs address of second pfo
     */
    static bool SortByNHits(const pandora::ParticleFlowObject *const pLhs, const pandora::ParticleFlowObject *const pRhs);

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
