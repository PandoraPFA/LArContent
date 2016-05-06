/**
 *  @file   LArContent/LArMonitoring/ParticleMonitoringAlgorithm.h
 *
 *  @brief  Header file for the particle monitoring algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_PARTICLE_MONITORING_ALGORITHM_H
#define LAR_PARTICLE_MONITORING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content
{

/**
 *  @brief  ParticleMonitoringAlgorithm class
 */
class ParticleMonitoringAlgorithm: public pandora::Algorithm
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

    /**
     *  @brief  Default constructor
     */
    ParticleMonitoringAlgorithm();

    /**
     *  @brief  Destructor
     */
    ~ParticleMonitoringAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string     m_caloHitListName;          ///< Name of input calo hit list
    std::string     m_mcParticleListName;       ///< Name of input MC particle list
    std::string     m_pfoListName;              ///< Name of input Pfo list
    std::string     m_fileName;                 ///< Name of output file
    std::string     m_treeName;                 ///< Name of output tree

    bool            m_primaryPfosOnly;          ///< Whether to extract only primary Pfos - top-level Pfos and top-level daughters of top-level neutrinos
    bool            m_collapseToPrimaryPfos;    ///< Whether to collapse hits associated with daughter pfos back to the primary pfo
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ParticleMonitoringAlgorithm::Factory::CreateAlgorithm() const
{
    return new ParticleMonitoringAlgorithm();
}

} // namespace lar_content

#endif // LAR_PARTICLE_MONITORING_ALGORITHM_H
