/**
 *  @file   LArContent/include/LArThreeDReco/LArEventBuilding/CosmicRayBuildingAlgorithm.h
 * 
 *  @brief  Header file for the cosmic-ray building algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_COSMIC_RAY_BUILDING_ALGORITHM_H
#define LAR_COSMIC_RAY_BUILDING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  CosmicRayBuildingAlgorithm class
 */
class CosmicRayBuildingAlgorithm : public pandora::Algorithm
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

inline pandora::Algorithm *CosmicRayBuildingAlgorithm::Factory::CreateAlgorithm() const
{
    return new CosmicRayBuildingAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_COSMIC_RAY_BUILDING_ALGORITHM_H
