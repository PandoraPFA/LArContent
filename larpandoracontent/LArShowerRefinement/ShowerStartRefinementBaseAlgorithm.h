/**
 *  @file   larpandoracontent/LArShowerRefinement/ShowerStartRefinementBaseAlgorithm.h
 *
 *  @brief  Header file for the shower characterisation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_SHOWER_START_REFINEMENT_BASE_ALGORITHM_H
#define LAR_SHOWER_START_REFINEMENT_BASE_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  ShowerStartRefinementBaseAlgorithm class
 */
class ShowerStartRefinementBaseAlgorithm : public pandora::Algorithm
{
private:
    virtual pandora::StatusCode Run();
    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef LAR_SHOWER_START_REFINEMENT_BASE_ALGORITHM_H
