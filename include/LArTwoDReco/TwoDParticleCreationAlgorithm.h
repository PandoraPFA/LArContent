/**
 *  @file   TwoDParticleCreationAlgorithm.h
 * 
 *  @brief  Header file for the two dimensional particle creation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_TWO_D_PARTICLE_CREATION_ALGORITHM_H
#define LAR_TWO_D_PARTICLE_CREATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  TwoDParticleCreationAlgorithm class
 */
class TwoDParticleCreationAlgorithm : public pandora::Algorithm
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

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Create pfos for provided clusters
     * 
     *  @param  pClusterList address of the cluster list
     */
    pandora::StatusCode CreatePFOs(const pandora::ClusterList *const pClusterList) const;

    std::string     m_inputClusterListNameU;                 ///< The input cluster list name for the U view
    std::string     m_inputClusterListNameV;                 ///< The input cluster list name for the V view
    std::string     m_inputClusterListNameW;                 ///< The input cluster list name for the W view

    std::string     m_outputPfoListName;                    ///< The output pfo list name
    unsigned int    m_minHitsInCluster;                     ///< Min number of hits for clusters to form pfos
    float           m_minClusterEnergy;                     ///< Min energy for clusters to form pfos
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *TwoDParticleCreationAlgorithm::Factory::CreateAlgorithm() const
{
    return new TwoDParticleCreationAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_TWO_D_PARTICLE_CREATION_ALGORITHM_H
