/**
 *  @file   LArContent/include/LArTwoDReco/LArCosmicRay/CosmicRayShowerMatchingAlgorithm.h
 * 
 *  @brief  Header file for the cosmic ray shower matching algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_COSMIC_RAY_SHOWER_MATCHING_ALGORITHM_H
#define LAR_COSMIC_RAY_SHOWER_MATCHING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  CosmicRayShowerMatchingAlgorithm class
 */
class CosmicRayShowerMatchingAlgorithm : public pandora::Algorithm
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
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CosmicRayShowerMatchingAlgorithm::Factory::CreateAlgorithm() const
{
    return new CosmicRayShowerMatchingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_COSMIC_RAY_SHOWER_MATCHING_ALGORITHM_H
