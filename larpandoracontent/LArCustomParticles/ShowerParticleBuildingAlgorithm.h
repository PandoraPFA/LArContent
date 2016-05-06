/**
 *  @file   LArContent/LArCustomParticles/ShowerParticleBuildingAlgorithm.h
 *
 *  @brief  Header file for the neutrino event creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_SHOWER_PARTICLE_BUILDING_ALGORITHM_H
#define LAR_SHOWER_PARTICLE_BUILDING_ALGORITHM_H 1

#include "larpandoracontent/LArContent/LArObjects/LArShowerPfo.h"

#include "larpandoracontent/LArContent/LArCustomParticles/CustomParticleCreationAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  ShowerParticleBuildingAlgorithm class
 */
class ShowerParticleBuildingAlgorithm : public CustomParticleCreationAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ShowerParticleBuildingAlgorithm();

    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void CreatePfo(const pandora::ParticleFlowObject *const pInputPfo, const pandora::ParticleFlowObject*& pOutputPfo) const;

    bool            m_cosmicMode;             ///<
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ShowerParticleBuildingAlgorithm::Factory::CreateAlgorithm() const
{
    return new ShowerParticleBuildingAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_SHOWER_PARTICLE_BUILDING_ALGORITHM_H
