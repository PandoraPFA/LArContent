/**
 *  @file   larpandoracontent/LArShowerRefinement/GammaStartRefinementAlgorithm.h
 *
 *  @brief  Header file for the shower characterisation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_GAMMA_START_REFINEMENT_ALGORITHM_H
#define LAR_GAMMA_START_REFINEMENT_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArShowerRefinement/ShowerStartRefinementBaseAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  GammaStartRefinementAlgorithm class
 */
class GammaStartRefinementAlgorithm : public ShowerStartRefinementBaseAlgorithm
{
private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef LAR_GAMMA_START_REFINEMENT_ALGORITHM_H
