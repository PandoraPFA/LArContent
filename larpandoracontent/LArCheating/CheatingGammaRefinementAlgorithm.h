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

    bool IsMatchedToGamma(const pandora::ParticleFlowObject *const pPfo) const;

    void FindContaminantMCParticles(const pandora::ParticleFlowObject *const pGammaPfo, pandora::MCParticleList &contaminantMCParticleList) const;

    void FilterContaminantMCParticleList(const pandora::ParticleFlowObject *const pGammaPfo, pandora::MCParticleList &mcContaminantList) const;

    bool IsShower(const pandora::MCParticle *const pMCParticle) const;

    bool AreShowersSeparated(const pandora::MCParticle *const pMCParticle1, const pandora::MCParticle *const pMCParticle2) const;

    void RemoveContaminantHits(const pandora::ParticleFlowObject *const pGammaPfo, const pandora::MCParticleList &mcContaminantList) const;

    const pandora::Cluster *CreateCluster(const pandora::MCParticle *const pMCParticle, const pandora::CaloHitList &caloHitList, const pandora::HitType hitType) const;

    std::string m_showerPfoListName;
    bool m_truncateMode;
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_CC_ELECTRON_REFINEMENT_ALGORITHM_H
