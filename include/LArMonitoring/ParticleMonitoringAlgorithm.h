/**
 *  @file   LArContent/include/LArMonitoring/ParticleMonitoringAlgorithm.h
 * 
 *  @brief  Header file for the particle monitoring algorithm.
 * 
 *  $Log: $
 */
#ifndef LAR_PARTICLE_MONITORING_ALGORITHM_H
#define LAR_PARTICLE_MONITORING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
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

    typedef std::map<pandora::Uid, int> ContributionMap;
    typedef std::map<pandora::Uid, pandora::Uid> UidRelationMap;
    typedef std::map<pandora::Uid, const pandora::MCParticle*> UidToMCParticleMap;
    typedef std::map<pandora::Uid, const pandora::ParticleFlowObject*> UidToPfoMap;

    /**
     *  @brief  Figure out the primary ancestor of each MC particle and store in a map
     * 
     *  @param  pMCParticleList
     *  @param  uidToPrimaryMap
     *  @param  uidToMCParticleMap
     */
    void GetMCParticleMaps(const pandora::MCParticleList *const pMCParticleList, UidRelationMap &uidToPrimaryMap,
        UidToMCParticleMap &uidToMCParticleMap) const;

    /**
     *  @brief  Find pfo containing most hits from each primary MC particle, and store properties in a map indexed by particle Uid
     * 
     *  @param  pCaloHitList
     *  @param  pfoList
     *  @param  uidToPrimaryMap
     *  @param  uidToPfoMap
     *  @param  contributionMap
     */
    void GetMCParticleToPfoMatches(const pandora::CaloHitList *const pCaloHitList, const pandora::PfoList &pfoList, 
        const UidRelationMap &uidToPrimaryMap, UidToPfoMap &uidToPfoMap, ContributionMap &contributionMap) const;

    /**
     *  @brief  Work out the total number of hits associated with each primary MC particle
     * 
     *  @param  pCaloHitList
     *  @param  uidToPrimaryMap
     *  @param  contributionMap
     */
    void GetMCParticleToCaloHitMatches(const pandora::CaloHitList *const pCaloHitList, const UidRelationMap &uidToPrimaryMap,
        ContributionMap &contributionMap) const;

    /**
     *  @brief  Get 3D position in best agreement with provided reference point, using combinations of u,v and w inner/outer layer centroids
     * 
     *  @param  clusterList
     *  @param  referencePoint
     *  @param  useInnerLayer
     */
    pandora::CartesianVector GetSpacePoint(const pandora::ClusterList &clusterList, const pandora::CartesianVector &referencePoint,
        const bool useInnerLayer) const;

    /**
     *  @brief  Whether a mc particle is a final-state particle from a neutrino (or antineutrino) interaction
     * 
     *  @param  pMCParticle
     * 
     *  @return boolean
     */
    bool IsNeutrinoInduced(const pandora::MCParticle *const pMCParticle) const;

    /**
     *  @brief  Whether a mc particle is a neutrino or (antineutrino)
     * 
     *  @param  pMCParticle
     * 
     *  @return boolean
     */
    bool IsNeutrino(const pandora::MCParticle *const pMCParticle) const;

    /**
     *  @brief  Get primary neutrino or antineutrino
     * 
     *  @param  pMCParticle
     * 
     *  @return pdg code of neutrino (or zero, otherwise)
     */
    int GetPrimaryNeutrino(const pandora::MCParticle *const pMCParticle) const;


    std::string     m_caloHitListName;                  ///< 
    std::string     m_mcParticleListName;               ///< 
    std::string     m_pfoListName;                      ///< 
    std::string     m_fileName;                         ///< 
    std::string     m_treeName;                         ///< 
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ParticleMonitoringAlgorithm::Factory::CreateAlgorithm() const
{
    return new ParticleMonitoringAlgorithm();
}

} // namespace lar

#endif // LAR_PARTICLE_MONITORING_ALGORITHM_H
