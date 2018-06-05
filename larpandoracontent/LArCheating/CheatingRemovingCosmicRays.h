/**
 *  @file   larpandoracontent/LArCheating/CheatingRemovingCosmicRays.h
 * 
 *  @brief  Header file for the cheating removing cosmic rays algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_CHEATING_REMOVING_COSMIC_RAYS_H
#define LAR_CHEATING_REMOVING_COSMIC_RAYS_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  CheatingRemovingCosmicRays::Algorithm class
 */
class CheatingRemovingCosmicRays : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CheatingRemovingCosmicRays() = default;

private:
    pandora::StatusCode Run();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_inputCaloHitListName;     ///< Input calo hit list name
    std::string m_outputCaloHitListName;    ///< Output calo hit list name
    std::string m_mcParticleListName;       ///< MC Particle list name
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_REMOVING_COSMIC_RAYS_H
