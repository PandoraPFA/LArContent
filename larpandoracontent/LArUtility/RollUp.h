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

    virtual const pandora::MCParticle *GetRollUpTargetMC(const pandora::MCParticle *const pMC) const = 0;

    virtual bool ShouldFoldCaloHit(const pandora::CaloHit *const pCaloHit, const pandora::MCParticle *const pRolledUpMainMC) const = 0;
};

class RollUpNullPolicy : public IRollUpPolicy
{
public:
    const pandora::MCParticle *GetRollUpTargetMC(const pandora::MCParticle *const pMC) const override;
    bool ShouldFoldCaloHit(const pandora::CaloHit *const pCaloHit, const pandora::MCParticle *const pRolledUpMainMC) const override;
};

class RollUpEMPolicy : public IRollUpPolicy
{
public:
    const pandora::MCParticle *GetRollUpTargetMC(const pandora::MCParticle *const pMC) const override;
    bool ShouldFoldCaloHit(const pandora::CaloHit *const pCaloHit, const pandora::MCParticle *const pRolledUpMainMC) const override;
};

class RollUpEMAndDeltaRayPolicy : public RollUpEMPolicy
{
public:
    const pandora::MCParticle *GetRollUpTargetMC(const pandora::MCParticle *const pMC) const override;
};

class RollUpEMAndAmbiguousDeltaRayHitsPolicy : public RollUpEMPolicy
{
public:
    RollUpEMAndAmbiguousDeltaRayHitsPolicy(
        const float deltaRayParentWeightThreshold, const std::map<pandora::HitType, float> deltaRayLengthThresholds);

    bool ShouldFoldCaloHit(const pandora::CaloHit *const pCaloHit, const pandora::MCParticle *const pRolledUpMainMC) const override;

private:
    float m_deltaRayParentWeightThreshold;
    std::map<pandora::HitType, float> m_deltaRayLengthThresholdsSquared;
};

class RollUpEMWithComptonFilterAndAmbiguousDeltaRayHitsPolicy : public RollUpEMAndAmbiguousDeltaRayHitsPolicy
{
public:
    using RollUpEMAndAmbiguousDeltaRayHitsPolicy::RollUpEMAndAmbiguousDeltaRayHitsPolicy;

    const pandora::MCParticle *GetRollUpTargetMC(const pandora::MCParticle *const pMC) const override;
};

class RollUpper
{
public:
    RollUpper();

    RollUpper(std::unique_ptr<IRollUpPolicy> policy);

    void Reset();

    const pandora::MCParticle *RollUpMC(const pandora::MCParticle *const pMC) const;

    const pandora::MCParticle *RollUpCaloHit(const pandora::CaloHit *const pCaloHit) const;

private:
    std::unique_ptr<IRollUpPolicy> m_policy;
    mutable std::unordered_map<const pandora::MCParticle*, const pandora::MCParticle*> m_mcCache;
    mutable std::unordered_map<const pandora::CaloHit*, const pandora::MCParticle*> m_caloHitCache;
};

} // namespace lar_content

#endif // #ifndef LAR_ROLL_UP_H
