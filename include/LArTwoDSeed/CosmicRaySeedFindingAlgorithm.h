/**
 *  @file   LArContent/include/LArTwoDSeed/CosmicRaySeedFindingAlgorithm.h
 *
 *  @brief  Header file for the cosmic-ray seed-finding algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_COSMIC_RAY_SEED_FINDING_ALGORITHM_H
#define LAR_COSMIC_RAY_SEED_FINDING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArTwoDSeed/SeedFindingBaseAlgorithm.h"

namespace lar
{

/**
 *  @brief  CosmicRaySeedFindingAlgorithm class
 */
class CosmicRaySeedFindingAlgorithm : public SeedFindingBaseAlgorithm
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
    void GetSeedClusterList(const pandora::ClusterVector &candidateClusters, pandora::ClusterList &seedClusterList) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CosmicRaySeedFindingAlgorithm::Factory::CreateAlgorithm() const
{
    return new CosmicRaySeedFindingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_COSMIC_RAY_SEED_FINDING_ALGORITHM_H
