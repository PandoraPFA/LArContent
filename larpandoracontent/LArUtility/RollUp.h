/**
 *  @file   larpandoracontent/LArUtility/RollUp.h
 *
 *  @brief  Header file for classes related to roll-up of EM activity. 
 *
 *  $Log: $
 */
#ifndef LAR_ROLL_UP_H
#define LAR_ROLL_UP_H 1

#include "Pandora/AlgorithmHeaders.h"

#include <memory>

namespace lar_content
{

/**
 *  @brief  The IRollUpPolicy class
 *
 *  Interface for policies that define how MC particles and calo hits are rolled up. Roll-up maps a given MC particle to an ancestor,
 *  collapsing EM shower or delta ray sub-trees into a single representative particle.
 */
class IRollUpPolicy
{
public:
    virtual ~IRollUpPolicy() = default;

    /**
     *  @brief  Get the MC particle a given MC particle should be rolled-up to
     *
     *  @param[in]  pMC  the input MC particle
     *
     *  @return  the target MC particle (may be the input MC particle itself if no roll-up applies)
     */
    virtual const pandora::MCParticle *GetRollUpTargetMC(const pandora::MCParticle *const pMC) const = 0;

    /**
     *  @brief  Whether a calo hit should be folded into the parent of its rolled-up main MC particle
     *
     *  @param[in]  pCaloHit         the calo hit
     *  @param[in]  pRolledUpMainMC  the main MC particle after MC particle level roll-up
     *
     *  @return  whether the hit should be further folded to the parent of hit's rolled-up MC particle
     */
    virtual bool ShouldFurtherRollUpCaloHit(
        const pandora::CaloHit *const pCaloHit, const pandora::MCParticle *const pRolledUpMainMC) const = 0;
};

/**
 *  @brief  The RollUpNullPolicy class
 *
 *  A no-op policy, no MC particles are rolled-up and no calo hits are folded.
 */
class RollUpNullPolicy : public IRollUpPolicy
{
public:
    const pandora::MCParticle *GetRollUpTargetMC(const pandora::MCParticle *const pMC) const override;
    bool ShouldFurtherRollUpCaloHit(
        const pandora::CaloHit *const pCaloHit, const pandora::MCParticle *const pRolledUpMainMC) const override;
};

/**
 *  @brief  The RollUpEMPolicy class
 *
 *  Roll up EM particles to the highest ancestor that is still EM.
 */
class RollUpEMPolicy : public IRollUpPolicy
{
public:
    const pandora::MCParticle *GetRollUpTargetMC(const pandora::MCParticle *const pMC) const override;
    bool ShouldFurtherRollUpCaloHit(
        const pandora::CaloHit *const pCaloHit, const pandora::MCParticle *const pRolledUpMainMC) const override;
};

/**
 *  @brief  The RollUpEMAndDeltaRayPolicy class
 *
 *  Extends RollUpEMPolicy to also roll up delta ray electrons to their parent track particle.
 */
class RollUpEMAndDeltaRayPolicy : public RollUpEMPolicy
{
public:
    const pandora::MCParticle *GetRollUpTargetMC(const pandora::MCParticle *const pMC) const override;
};

/**
 *  @brief  The RollUpEMAndAmbiguousDeltaRayHitsPolicy class
 *
 *  Extends RollUpEMPolicy to handle ambiguous delta ray hits at the calo hit level, rolling up parts of delta rays that
 *  significantly overlap with the parent track in an attempt to stop discontinuities in the parent track.
 */
class RollUpEMAndAmbiguousDeltaRayHitsPolicy : public RollUpEMPolicy
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param[in]  deltaRayParentWeightThreshold  minimum weight contribution of the parent track to a hit for that hit to be folded
     *  @param[in]  deltaRayLengthThresholds       per-view length thresholds below which a non-showering delta ray is folded
     *                                             to its parent track
     */
    RollUpEMAndAmbiguousDeltaRayHitsPolicy(
        const float deltaRayParentWeightThreshold, const std::map<pandora::HitType, float> deltaRayLengthThresholds);

    bool ShouldFurtherRollUpCaloHit(
        const pandora::CaloHit *const pCaloHit, const pandora::MCParticle *const pRolledUpMainMC) const override;

private:
    float m_deltaRayParentWeightThreshold;                               ///< Minimum parent weight for a hit to be folded to the parent track
    std::map<pandora::HitType, float> m_deltaRayLengthThresholdsSquared; ///< Squared per-view length thresholds for short delta ray folding
};

/**
 *  @brief  The RollUpEMWithComptonFilterAndAmbiguousDeltaRayHitsPolicy class
 *
 *  A concisely named class that extends RollUpEMAndAmbiguousDeltaRayHitsPolicy with a Compton scatter filter.
 *  Pure Compton scatter chains (photon → electron(s) with succeeding Brems) are blocked from roll-up,
 *  preventing potential associations between distant and diffuse hits.
 */
class RollUpEMWithComptonFilterAndAmbiguousDeltaRayHitsPolicy : public RollUpEMAndAmbiguousDeltaRayHitsPolicy
{
public:
    using RollUpEMAndAmbiguousDeltaRayHitsPolicy::RollUpEMAndAmbiguousDeltaRayHitsPolicy;

    const pandora::MCParticle *GetRollUpTargetMC(const pandora::MCParticle *const pMC) const override;
};

/**
 *  @brief  The RollUpper class
 *
 *  Applies a configurable roll-up policy to MC particles and calo hits using per-event caching to save compute.
 *  Construct with a IRollUpPolicy to control the roll-up logic. Reset() must be called between events to clear the caches.
 */
class RollUpper
{
public:
    /**
     *  @brief  Default constructor, uses the no roll-up policy
     */
    RollUpper();

    /**
     *  @brief  Constructor
     *
     *  @param[in]  policy  the roll-up policy to apply
     */
    RollUpper(std::unique_ptr<IRollUpPolicy> policy);

    /**
     *  @brief  Clear the per-event MC particle and calo hit caches
     */
    void Reset() const;

    /**
     *  @brief  Get the MC particle that the input MC particle rolls up to under the current policy
     *
     *  @param[in]  pMC  the input MC particle
     *
     *  @return  the rolled-up MC particle (may be the input MC Particle itself)
     */
    const pandora::MCParticle *RollUpMC(const pandora::MCParticle *const pMC) const;

    /**
     *  @brief  Get the main rolled-up MC particle for a calo hit
     *
     *  Aggregates the hit's MC particle weight map after rolling up each contributing MC particle, selects the highest-weight result,
     *  then applies any hit-level folding defined by the policy.
     *
     *  @param[in]  pCaloHit  the input calo hit
     *
     *  @return  the main rolled-up MC particle for the hit (may be null if the hit's weight map is empty)
     */
    const pandora::MCParticle *RollUpCaloHit(const pandora::CaloHit *const pCaloHit) const;

private:
    std::unique_ptr<IRollUpPolicy> m_policy;                                                          ///< The roll-up policy
    mutable std::unordered_map<const pandora::MCParticle *, const pandora::MCParticle *> m_mcCache;   ///< Per-event cache of MC roll-up results
    mutable std::unordered_map<const pandora::CaloHit *, const pandora::MCParticle *> m_caloHitCache; ///< Per-event cache of calo hit roll-up results
};

} // namespace lar_content

#endif // #ifndef LAR_ROLL_UP_H
