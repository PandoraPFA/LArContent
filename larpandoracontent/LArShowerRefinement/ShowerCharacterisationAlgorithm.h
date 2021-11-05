/**
 *  @file   larpandoracontent/LArShowerRefinement/ShowerCharacterisationAlgorithm.h
 *
 *  @brief  Header file for the shower characterisation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_SHOWER_CHARACTERISATION_ALGORITHM_H
#define LAR_SHOWER_CHARACTERISATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  ShowerCharacterisationAlgorithm class
 */
class ShowerCharacterisationAlgorithm : public pandora::Algorithm
{
private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef LAR_SHOWER_CHARACTERISATION_ALGORITHM_H
