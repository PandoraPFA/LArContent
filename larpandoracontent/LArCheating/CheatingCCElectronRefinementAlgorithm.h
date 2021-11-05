/**
 *  @file   larpandoracontent/LArCheating/CheatingCCElectronRefinementAlgorithm.h
 *
 *  @brief  Header file for the cheating vertex creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_CC_ELECTRON_REFINEMENT_ALGORITHM_H
#define LAR_CHEATING_CC_ELECTRON_REFINEMENT_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content
{

/**
 *  @brief  CheatingCCElectronRefinementAlgorithm::Algorithm class
 */
class CheatingCCElectronRefinementAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CheatingCCElectronRefinementAlgorithm();

    typedef std::unordered_map<const pandora::CaloHit*, const pandora::Cluster*> HitToClusterMap;
    typedef std::unordered_map<const pandora::Cluster*, const pandora::ParticleFlowObject*> ClusterToPfoMap;
    typedef std::unordered_map<const pandora::MCParticle*, pandora::PfoList> MCParticleToPfoMap;

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void FillElectronHitMap(const pandora::CaloHitList *const pCaloHitList, LArMCParticleHelper::MCContributionMap &mcElectronToHitMap) const;

    const pandora::ParticleFlowObject *FindBestPfoMatch(const pandora::MCParticle *const pMCElectron, const pandora::PfoList *const pShowerPfoList, 
        LArMCParticleHelper::MCContributionMap &mcElectronToHitMap) const;

    void FillOwnershipMaps(const pandora::PfoList *const pShowerPfoList);

    void RefineElectronPfos(const pandora::ParticleFlowObject *const pElectronPfo, const pandora::MCParticle *const pMCElectron, 
        LArMCParticleHelper::MCContributionMap &mcElectronToHitMap);

    float FindConeLength(const pandora::ParticleFlowObject *const pPfo, const pandora::CartesianVector &mcVertex, 
        const pandora::CartesianVector &mcDirection, const pandora::HitType hitType) const;

    bool DoesPassCut(const pandora::CaloHit *const pCaloHit, const pandora::CartesianVector &mcVertex, const pandora::CartesianVector &mcDirection, 
        const float coneLength) const;

    std::string m_mcParticleListName;
    std::string m_showerPfoListName;
    std::string m_trackPfoListName;
    std::string m_caloHitListName;
    pandora::StringVector m_clusterListNames; 

    bool m_completenessMode;
    float m_thresholdCompleteness;
    float m_thresholdPurity;
    float m_maxOpeningAngle;

    HitToClusterMap m_hitToClusterMapU;
    HitToClusterMap m_hitToClusterMapV;
    HitToClusterMap m_hitToClusterMapW;
    ClusterToPfoMap m_clusterToPfoMapU;
    ClusterToPfoMap m_clusterToPfoMapV;
    ClusterToPfoMap m_clusterToPfoMapW;
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_CC_ELECTRON_REFINEMENT_ALGORITHM_H
