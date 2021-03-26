/**
 *  @file   larpandoracontent/LArCheating/CheatingCosmicRayRemovalAlgorithm.h
 *
 *  @brief  Header file for the cheating cosmic ray removal algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_COSMIC_RAY_REMOVAL_ALGORITHM_H
#define LAR_CHEATING_COSMIC_RAY_REMOVAL_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  CheatingCosmicRayRemovalAlgorithm::Algorithm class
 */
class CheatingCosmicRayRemovalAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CheatingCosmicRayRemovalAlgorithm() = default;

private:
    pandora::StatusCode Run();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_inputCaloHitListName;  ///< Input calo hit list name
    std::string m_outputCaloHitListName; ///< Output calo hit list name
    std::string m_mcParticleListName;    ///< MC Particle list name
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_COSMIC_RAY_REMOVAL_ALGORITHM_H
