/**
 *  @file   
 *
 *  @brief  
 *
 *  $Log: $
 */
#ifndef LAR_ROLL_UP_H
#define LAR_ROLL_UP_H 1

#include "Pandora/AlgorithmHeaders.h"

#include <memory>

namespace lar_content
{

class IRollUpPolicy
{
public:
    virtual ~IRollUpPolicy() = default;

    virtual bool ShouldRollUpMC(const pandora::MCParticle *const pMC) const = 0;

    virtual bool ShouldRollUpCaloHit(const pandora::CaloHit *const pCaloHit, const pandora::MCParticle *const pRolledUpMainMC) const = 0;
};

class RollUpEMPolicy : public IRollUpPolicy
{
public:
    bool ShouldRollUpMC(const pandora::MCParticle *const pMC) const override;
    bool ShouldRollUpCaloHit(const pandora::CaloHit *const pCaloHit, const pandora::MCParticle *const pRolledUpMainMC) const override;
};

class RollUpEMAndDeltaRayPolicy : public IRollUpPolicy
{
public:
    bool ShouldRollUpMC(const pandora::MCParticle *const pMC) const override;
    bool ShouldRollUpCaloHit(const pandora::CaloHit *const pCaloHit, const pandora::MCParticle *const pRolledUpMainMC) const override;
};

class RollUpEMAndAmbiguousDeltaRayHitsPolicy : public IRollUpPolicy
{
public:
    RollUpEMAndAmbiguousDeltaRayHitsPolicy(
        const float deltaRayParentWeightThreshold, const std::map<pandora::HitType, float> deltaRayLengthThresholds);

    bool ShouldRollUpMC(const pandora::MCParticle *const pMC) const override;
    bool ShouldRollUpCaloHit(const pandora::CaloHit *const pCaloHit, const pandora::MCParticle *const pRolledUpMainMC) const override;

private:
    bool CausesShower(const pandora::MCParticle *const pMC, int nDescendantElectrons) const;

    float m_deltaRayParentWeightThreshold;
    std::map<pandora::HitType, float> m_deltaRayLengthThresholdsSquared;
};

class RollUpper
{
public:
    RollUpper(std::unique_ptr<IRollUpPolicy> policy);

    const pandora::MCParticle *RollUpMC(const pandora::MCParticle *const pMC);

    const pandora::MCParticle *RollUpCaloHit(const pandora::CaloHit *const pCaloHit);

private:
    std::unique_ptr<IRollUpPolicy> m_policy;
    std::unordered_map<const pandora::MCParticle*, const pandora::MCParticle*> m_mcCache;
    std::unordered_map<const pandora::CaloHit*, const pandora::MCParticle*> m_caloHitCache;
};

} // namespace lar_content

#endif // #ifndef LAR_ROLL_UP_H
