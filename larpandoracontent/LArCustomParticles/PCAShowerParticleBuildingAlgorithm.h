/**
 *  @file   larpandoracontent/LArCustomParticles/PCAShowerParticleBuildingAlgorithm.h
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
 *  @brief  PCAShowerParticleBuildingAlgorithm class
 */
class PCAShowerParticleBuildingAlgorithm : public CustomParticleCreationAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    PCAShowerParticleBuildingAlgorithm();

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

    typedef std::vector<pandora::CartesianVector> EigenVectors;

    void RunPCA(const pandora::CaloHitList &threeDCaloHitList, pandora::CartesianVector &centroid, pandora::CartesianVector &outputEigenValues, EigenVectors &outputEigenVecs) const;

    pandora::CartesianVector ShowerLength(const pandora::CartesianVector &eigenValues) const;

    float OpeningAngle(const pandora::CartesianVector &principal, const pandora::CartesianVector &secondary, const pandora::CartesianVector &eigenValues) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool            m_cosmicMode;               ///<
    unsigned int    m_layerFitHalfWindow;       ///< 
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *PCAShowerParticleBuildingAlgorithm::Factory::CreateAlgorithm() const
{
    return new PCAShowerParticleBuildingAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_PCA_SHOWER_PARTICLE_BUILDING_ALGORITHM_H
