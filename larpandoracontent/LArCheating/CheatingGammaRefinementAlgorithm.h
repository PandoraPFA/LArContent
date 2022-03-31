/**
 *  @file   larpandoracontent/LArCheating/CheatingGammaRefinementAlgorithm.h
 *
 *  @brief  Header file for the cheating vertex creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_GAMMA_REFINEMENT_ALGORITHM_H
#define LAR_CHEATING_GAMMA_REFINEMENT_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include <exception>
#include <stdexcept>

namespace lar_content
{

/**
 *  @brief  CheatingGammaRefinementAlgorithm::Algorithm class
 */
class CheatingGammaRefinementAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CheatingGammaRefinementAlgorithm();

    typedef std::unordered_map<const pandora::MCParticle*, pandora::CaloHitList> MCParticleToHitListMap;

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void FillGammaHitMap(const pandora::CaloHitList *const pCaloHitList, LArMCParticleHelper::MCContributionMap &gammaHitMap) const;

    const pandora::MCParticle *FindMatchedGamma(const pandora::ParticleFlowObject *const pPfo, const LArMCParticleHelper::MCContributionMap &gammaHitMap) const;

    void FindContaminantMCParticles(const pandora::MCParticle *const pMCGamma, const pandora::ParticleFlowObject *const pGammaPfo, pandora::MCParticleList &contaminantMCParticleList) const;

    void FilterContaminantMCParticleList(const pandora::MCParticle *const pMCGamma, const pandora::ParticleFlowObject *const pGammaPfo, pandora::MCParticleList &mcContaminantList) const;

    bool IsShower(const pandora::MCParticle *const pMCParticle) const;

    bool AreParticlesSeparated(const pandora::MCParticle *const pMCParticle1, const pandora::MCParticle *const pMCParticle2) const;

    void GetDistinguishableChildren(const pandora::MCParticle *const pMCContaminant, const pandora::MCParticle *const pMCGamma, 
        pandora::MCParticleList &mcHierarchyList) const;

    void RemoveContaminantHits(const pandora::MCParticle *const pMCGamma, const pandora::ParticleFlowObject *const pGammaPfo, const pandora::MCParticleList &mcContaminantList,  
        const LArMCParticleHelper::MCContributionMap &mcParticleHitMap) const;

    const pandora::Cluster *CreateCluster(const pandora::MCParticle *const pMCParticle, const pandora::CaloHitList &caloHitList, const pandora::HitType hitType) const;

    void handle_eptr(std::exception_ptr eptr);

    std::string m_caloHitListName;
    std::string m_showerPfoListName;
    float m_maxImpactT;
    float m_minOpeningAngle;
    float m_minGammaCompleteness;
    bool m_removeHierarchy;
    bool m_truncateMode;
    float m_creationCompletenessThreshold;
    bool m_isActive;
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_CC_ELECTRON_REFINEMENT_ALGORITHM_H
