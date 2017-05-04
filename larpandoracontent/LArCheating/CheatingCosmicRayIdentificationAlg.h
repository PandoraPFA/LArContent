/**
 *  @file   larpandoracontent/LArCheating/CheatingCosmicRayIdentificationAlg.h
 * 
 *  @brief  Header file for the cosmic ray identification cheater class.
 * 
 *  $Log: $
 */
#ifndef LAR_CHEATING_COSMIC_RAY_IDENTIFICATION_ALG_H
#define LAR_CHEATING_COSMIC_RAY_IDENTIFICATION_ALG_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  CheatingCosmicRayIdentificationAlg class
 */
class CheatingCosmicRayIdentificationAlg : public pandora::Algorithm
{
private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string     m_inputPfoListName;             ///< The input pfo list name
    std::string     m_outputPfoListName;            ///< The output pfo list name
    std::string     m_inputDaughterPfoListName;     ///< The input daughter pfo list name (if not specified, will assume same as main input list)
    std::string     m_outputDaughterPfoListName;    ///< The output daughter pfo list name (if not specified, will assume same as main output list)
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_COSMIC_RAY_IDENTIFICATION_ALG_H
