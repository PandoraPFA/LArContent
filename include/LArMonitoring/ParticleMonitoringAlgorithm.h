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

#include "LArHelpers/LArMCParticleHelper.h"

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

    typedef std::unordered_map<const pandora::MCParticle*, const pandora::MCParticle*> MCRelationMap;
    typedef std::unordered_map<const pandora::MCParticle*, const pandora::ParticleFlowObject*> MCToPfoMap;
    typedef std::unordered_map<const pandora::ParticleFlowObject*, const pandora::MCParticle*> PfoToMCMap;
    typedef std::unordered_map<const pandora::CaloHit*, const pandora::MCParticle*> CaloHitToMCMap;
    typedef std::unordered_map<const pandora::CaloHit*, const pandora::ParticleFlowObject*> CaloHitToPfoMap;
    typedef std::unordered_map<const pandora::ParticleFlowObject*, pandora::CaloHitList> PfoContributionMap;
    typedef std::unordered_map<const pandora::MCParticle*, pandora::CaloHitList> MCContributionMap;

    /**
     *  @brief  Get neutrino MC particles from an input MC particle list
     */
    void GetTrueNeutrinos(const pandora::MCParticleList *pMCParticleList, pandora::MCParticleList &trueNeutrinos) const;

    /**
     *  @brief  Get neutrino pfos from an input pfo list
     */
    void GetRecoNeutrinos(const pandora::PfoList *pPfoList, pandora::PfoList &recoNeutrinos) const;

    /**
     *  @brief  Modify a pfo list, recursively removing top-level neutrinos and replacing them with their daughter pfos
     *
     *  @param  pfoList the pfo list, which may be modified
     */
    void ExtractNeutrinoDaughters(pandora::PfoList &pfoList) const;

    /**
     *  @brief  Create map from true neutrinos to reconstructed neutrinos
     *
     *  @param  pCaloHitList the input list of calo hits
     *  @param  recoNeutrinos the input list of reconstructed neutrinos
     *  @param  mcHitMap the input mapping between calo hits and their main MC particles
     *  @param  outputPrimaryMap the ouput mapping between from true to reconstructed neutrinos
     */
    void GetNeutrinoMatches(const pandora::CaloHitList *const pCaloHitList, const pandora::PfoList &recoNeutrinos,
        const CaloHitToMCMap &mcHitMap, MCToPfoMap &outputNeutrinoMap) const;

    /**
     *  @brief  Match calo hits to their parent particles
     *
     *  @param  pCaloHitList the input list of calo hits
     *  @param  mcPrimaryMap the input map between particles and their primaries
     *  @param  mcHitMap the output mapping between calo hits and their main MC particles
     *  @param  mcContributionMap the output mapping between MC particles and their number of associated hits
     */
    void GetMCParticleToCaloHitMatches(const pandora::CaloHitList *const pCaloHitList, const LArMCParticleHelper::MCRelationMap &mcPrimaryMap,
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
     *  @param  matchedPrimaryMap the output mapping between MC particles to their best matched Pfo
     *  @param  matchedContributionMap the output mapping between MC particles and their number of matched hits
     *  @param  fullContributionMap the output mapping between MC particles and any pfo hits
     */
    void GetMCParticleToPfoMatches(const pandora::CaloHitList *const pCaloHitList, const pandora::PfoList &pfoList,
        const CaloHitToMCMap &mcHitMap, MCToPfoMap &matchedPrimaryMap, MCContributionMap &matchedContributionMap,
        MCContributionMap &fullContributionMap) const;

    /**
     *  @brief  Collect up all calo hits associated with a Pfo and its daughters
     *
     *  @param  pParentPfo the input Pfo
     *  @param  caloHitList the output calo hit list
     */
    void CollectCaloHits(const pandora::ParticleFlowObject *const pParentPfo, pandora::CaloHitList &caloHitList) const;


    void CollectCaloHits(const pandora::PfoList &pfoList, pandora::CaloHitList &caloHitList) const;


    /**
     *  @brief  Collect up all daughter clusters of a Pfo
     *
     *  @param  pParentPfo the input Pfo
     *  @param  clusterList the output cluster list
     */
    void CollectDaughterClusters(const pandora::ParticleFlowObject *const pParentPfo, pandora::ClusterList &clusterList) const;

    /**
     *  @brief  Count the number of calo hits, in a provided list, of a specified type
     *
     *  @param  hitType the hit type
     *  @param  caloHitList the calo hit list
     *
     *  @return the number of calo hits of the specified type
     */
    unsigned int CountHitsByType(const pandora::HitType hitType, const pandora::CaloHitList &caloHitList) const;

    std::string     m_caloHitListName;          ///< Name of input calo hit list
    std::string     m_mcParticleListName;       ///< Name of input MC particle list
    std::string     m_pfoListName;              ///< Name of input Pfo list
    std::string     m_fileName;                 ///< Name of output file
    std::string     m_treeName;                 ///< Name of output tree

    bool            m_useDaughterPfos;          ///< Whether to include daughter pfos in performance metrics
    bool            m_extractNeutrinoDaughters; ///< Whether to treat each neutrino pfo daughter as a standalone top-level pfo
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ParticleMonitoringAlgorithm::Factory::CreateAlgorithm() const
{
    return new ParticleMonitoringAlgorithm();
}

} // namespace lar_content

#endif // LAR_PARTICLE_MONITORING_ALGORITHM_H
