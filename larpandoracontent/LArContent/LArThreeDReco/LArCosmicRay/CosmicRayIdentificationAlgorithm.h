/**
 *  @file   LArContent/include/LArThreeDReco/LArCosmicRay/CosmicRayIdentificationAlgorithm.h
 * 
 *  @brief  Header file for the cosmic ray identification algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_COSMIC_RAY_IDENTIFICATION_ALGORITHM_H
#define LAR_COSMIC_RAY_IDENTIFICATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  CosmicRayIdentificationAlgorithm class
 */
class CosmicRayIdentificationAlgorithm : public pandora::Algorithm
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

    std::string     m_inputPfoListName;             ///< The input pfo list name
    std::string     m_outputPfoListName;            ///< The output pfo list name
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CosmicRayIdentificationAlgorithm::Factory::CreateAlgorithm() const
{
    return new CosmicRayIdentificationAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_COSMIC_RAY_IDENTIFICATION_ALGORITHM_H
