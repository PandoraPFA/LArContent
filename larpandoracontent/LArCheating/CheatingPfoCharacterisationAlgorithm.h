/**
 *  @file   larpandoracontent/LArCheating/CheatingPfoCharacterisationAlgorithm.h
 * 
 *  @brief  Header file for the cheating pfo characterisation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_CHEATING_PFO_CHARACTERISATION_ALGORITHM_H
#define LAR_CHEATING_PFO_CHARACTERISATION_ALGORITHM_H 1

#include "larpandoracontent/LArTrackShowerId/PfoCharacterisationAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  CheatingPfoCharacterisationAlgorithm class
 */
class CheatingPfoCharacterisationAlgorithm : public PfoCharacterisationAlgorithm
{
private:
    bool IsClearTrack(const pandora::ParticleFlowObject *const pPfo) const;
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_PFO_CHARACTERISATION_ALGORITHM_H
