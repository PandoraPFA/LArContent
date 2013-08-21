/**
 *  @file   LArContent/include/LArTwoDSeed/SeedLengthGrowingAlgorithm.h
 * 
 *  @brief  Header file for the seed length growing algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_SEED_LENGTH_GROWING_ALGORITHM_H
#define LAR_SEED_LENGTH_GROWING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArObjects/LArPointingCluster.h"

#include "LArTwoDSeed/SeedGrowingAlgorithm.h"

namespace lar
{

/**
 *  @brief  SeedLengthGrowingAlgorithm class
 */
class SeedLengthGrowingAlgorithm : public SeedGrowingAlgorithm
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

    void GetCandidateClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;
    AssociationType AreClustersAssociated(const pandora::Cluster *const pClusterSeed, const pandora::Cluster *const pCluster) const;

    typedef std::map<const pandora::Cluster*, LArPointingCluster> PointingClusterMap;
    PointingClusterMap  m_pointingClusterMap;       ///< The pointing cluster map
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *SeedLengthGrowingAlgorithm::Factory::CreateAlgorithm() const
{
    return new SeedLengthGrowingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_SEED_LENGTH_GROWING_ALGORITHM_H
