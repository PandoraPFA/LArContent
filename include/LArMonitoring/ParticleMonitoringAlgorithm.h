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

    typedef std::map<const pandora::MCParticle*, const pandora::MCParticle*> MCRelationMap;
    typedef std::map<const pandora::MCParticle*, const pandora::ParticleFlowObject*> MCToPfoMap;
    typedef std::map<const pandora::ParticleFlowObject*, const pandora::MCParticle*> PfoToMCMap;
    typedef std::map<const pandora::CaloHit*, const pandora::MCParticle*> CaloHitToMCMap;
    typedef std::map<const pandora::CaloHit*, const pandora::ParticleFlowObject*> CaloHitToPfoMap;
    typedef std::map<const pandora::ParticleFlowObject*, int> PfoContributionMap;
    typedef std::map<const pandora::MCParticle*, int> MCContributionMap;

    /**
     *  @brief  Create map between each MC particle and its primary parent MC particle
     *
     *  @param  pMCParticleList the input MC particle list
     *  @param  mcPrimaryMap the output mapping between particles and their parents
     */
    void GetMCParticleMaps(const pandora::MCParticleList *const pMCParticleList, MCRelationMap &mcPrimaryMap) const;

    /**
     *  @brief  Match calo hits to their parent particles
     *
     *  @param  pCaloHitList the input list of calo hits
     *  @param  mcPrimaryMap the input map between particles and their primaries
     *  @param  mcHitMap the output mapping between calo hits and their main MC particles
     *  @param  mcContributionMap the output mapping between MC particles and their number of associated hits
     */
    void GetMCParticleToCaloHitMatches(const pandora::CaloHitList *const pCaloHitList, const MCRelationMap &mcPrimaryMap,
        CaloHitToMCMap &mcHitMap, MCContributionMap &mcContributionMap) const;

    /**
     *  @brief  Match calo hits to their parent Pfos
     *
     *  @param  pCaloHitList the input list of calo hits
     *  @param  pfoList the input list of Pfos
     *  @param  pfoHitMap the output mapping between calo hits and their parent Pfos
     *  @param  pfoContributionMap the output mapping between Pfos and their number of associated hits
     */
    void GetPfoToCaloHitMatches(const pandora::CaloHitList *const pCaloHitList, const pandora::PfoList &pfoList,
        CaloHitToPfoMap &pfoHitMap, PfoContributionMap &pfoContributionMap) const;

    /**
     *  @brief  Match MC particle and Pfos
     *
     *  @param  pCaloHitList the input list of calo hits
     *  @param  pfoList the input list of Pfos
     *  @param  mcHitMap the input mapping between calo hits and their parent Pfos
     *  @param  outputPrimaryMap the output mapping between MC particles to their best matched Pfo
     *  @param  outputContributionMap the output mapping between MC particles and their number of matched hits
     */
    void GetMCParticleToPfoMatches(const pandora::CaloHitList *const pCaloHitList, const pandora::PfoList &pfoList,
        const CaloHitToMCMap &mcHitMap, MCToPfoMap &outputPrimaryMap, MCContributionMap &outputContributionMap) const;

    /**
     *  @brief  Collect up all calo hits associated with a Pfo and its daughters
     *
     *  @param  pParentPfo the input Pfo
     *  @param  caloHitList the output calo hit list
     */
    void CollectCaloHits(const pandora::ParticleFlowObject *const pParentPfo, pandora::CaloHitList &caloHitList) const;

    /**
     *  @brief  Collect up all daughter clusters of a Pfo
     *
     *  @param  pParentPfo the input Pfo
     *  @param  clusterList the output cluster list
     */
    void CollectDaughterClusters(const pandora::ParticleFlowObject *const pParentPfo, pandora::ClusterList &clusterList) const;

    std::string     m_caloHitListName;                  ///< Name of input calo hit list
    std::string     m_mcParticleListName;               ///< Name of input MC particle list
    std::string     m_pfoListName;                      ///< Name of input Pfo list
    std::string     m_fileName;                         ///< Name of output file
    std::string     m_treeName;                         ///< Name of output tree

    bool            m_useDaughterPfos;                  ///< Whether to include daughter Pfos in performance metrics
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ParticleMonitoringAlgorithm::Factory::CreateAlgorithm() const
{
    return new ParticleMonitoringAlgorithm();
}

} // namespace lar

#endif // LAR_PARTICLE_MONITORING_ALGORITHM_H
