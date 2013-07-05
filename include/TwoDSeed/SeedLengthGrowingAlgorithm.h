/**
 *  @file   SeedLengthGrowingAlgorithm.h
 * 
 *  @brief  Header file for the seed length growing algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_SEED_LENGTH_GROWING_ALGORITHM_H
#define LAR_SEED_LENGTH_GROWING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArPointingCluster.h"
#include "SeedGrowingAlgorithm.h"

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
