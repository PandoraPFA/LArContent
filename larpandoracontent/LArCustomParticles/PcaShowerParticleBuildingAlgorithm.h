/**
 *  @file   larpandoracontent/LArCustomParticles/PcaShowerParticleBuildingAlgorithm.h
 *
 *  @brief  Header file for the neutrino event creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_PCA_SHOWER_PARTICLE_BUILDING_ALGORITHM_H
#define LAR_PCA_SHOWER_PARTICLE_BUILDING_ALGORITHM_H 1

#include "larpandoracontent/LArObjects/LArShowerPfo.h"

#include "larpandoracontent/LArCustomParticles/CustomParticleCreationAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  PcaShowerParticleBuildingAlgorithm class
 */
class PcaShowerParticleBuildingAlgorithm : public CustomParticleCreationAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    PcaShowerParticleBuildingAlgorithm();

    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

private:
    void CreatePfo(const pandora::ParticleFlowObject *const pInputPfo, const pandora::ParticleFlowObject *&pOutputPfo) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int m_layerFitHalfWindow; ///<
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *PcaShowerParticleBuildingAlgorithm::Factory::CreateAlgorithm() const
{
    return new PcaShowerParticleBuildingAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_PCA_SHOWER_PARTICLE_BUILDING_ALGORITHM_H
