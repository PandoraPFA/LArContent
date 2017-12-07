/**
 *  @file   larpandoracontent/LArMonitoring/PfoValidationAlgorithm.h
 *
 *  @brief  Header file for the pfo validation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_PFO_VALIDATION_ALGORITHM_H
#define LAR_PFO_VALIDATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  PfoValidationAlgorithm class
 */
class PfoValidationAlgorithm: public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    PfoValidationAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string                                m_caloHitListName;          ///< Name of input calo hit list
    std::string                                m_pfoListName;              ///< Name of input pfo list
    LArMCParticleHelper::ValidationParameters  m_parameters;               ///< Parameters used to decide when an MCParticle is reconstructable
};

} // namespace lar_content

#endif // LAR_PFO_VALIDATION_ALGORITHM_H
