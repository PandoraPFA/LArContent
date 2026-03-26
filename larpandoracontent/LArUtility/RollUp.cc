/**
 *  @file   
 *
 *  @brief  
 *
 *  $Log: $
 */

#include "larpandoracontent/LArUtility/RollUp.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

using namespace pandora;

namespace lar_content
{

RollUpper::RollUpper(std::unique_ptr<IRollUpPolicy> policy) :
    m_policy{std::move(policy)}
{
}

const MCParticle *RollUpper::RollUpMC(const MCParticle *const pMC)
{
    if (m_mcCache.find(pMC) != m_mcCache.end())
        return m_mcCache.at(pMC);

    const MCParticle *pLeadingMC{pMC};
    while (!pLeadingMC->IsRootParticle())
    {
        if (!m_policy->ShouldRollUpMC(pLeadingMC))
            break;
        pLeadingMC = LArMCParticleHelper::GetNextParentMCParticle(pLeadingMC);
    }

    m_mcCache.insert({pMC, pLeadingMC});

    return pLeadingMC;
}

const MCParticle *RollUpper::RollUpCaloHit(const CaloHit *const pCaloHit)
{
    if (m_caloHitCache.find(pCaloHit) != m_caloHitCache.end())
        return m_caloHitCache.at(pCaloHit);

    const MCParticleWeightMap &weightMap{pCaloHit->GetMCParticleWeightMap()};
    MCParticleWeightMap rollUpWeightMap;
    for (const auto &[pMC, weight] : weightMap)
        rollUpWeightMap[this->RollUpMC(pMC)] += weight;

    const MCParticle *pMainMC{nullptr};
    float maxWeight{0.f};
    for (const auto &[pMC, weight] : rollUpWeightMap)
    {
        if (weight > maxWeight)
        {
            pMainMC = pMC;
            maxWeight = weight;
        }
        else if (weight == maxWeight) // tie-breaker (very unlikely)
        {
            if (LArMCParticleHelper::SortByMomentum(pMC, pMainMC))
                pMainMC = pMC;
        }
    }

    if (pMainMC && m_policy->ShouldRollUpCaloHit(pCaloHit, pMainMC))
        pMainMC = LArMCParticleHelper::GetNextParentMCParticle(pMainMC);

    m_caloHitCache.insert({pCaloHit, pMainMC});

    return pMainMC;
}

bool RollUpEMPolicy::ShouldRollUpMC(const MCParticle *const pMC) const
{
    if (pMC->IsRootParticle() || !LArMCParticleHelper::IsEM(pMC))
        return false;

    const MCParticle *const pParentMC{LArMCParticleHelper::GetNextParentMCParticle(pMC)};

    if (LArMCParticleHelper::IsEM(pParentMC))
        return true;

    return false;
}

bool RollUpEMPolicy::ShouldRollUpCaloHit(
    [[maybe_unused]] const CaloHit *const pCaloHit, [[maybe_unused]] const MCParticle *const pRolledUpMainMC) const
{
    return false;
}

bool RollUpEMAndDeltaRayPolicy::ShouldRollUpMC(const MCParticle *const pMC) const
{
    if (pMC->IsRootParticle() || !LArMCParticleHelper::IsEM(pMC))
        return false;

    const MCParticle *const pParentMC{LArMCParticleHelper::GetNextParentMCParticle(pMC)};

    if (LArMCParticleHelper::IsEM(pParentMC))
        return true;

    // is electron and parent is not EM and parent is charged -> is delta ray -> roll-up
    const bool parentIsCharged{PdgTable::GetParticleCharge(pParentMC->GetParticleId()) != 0};
    if (pMC->GetParticleId() == E_MINUS && parentIsCharged)
        return true;

    return false;
}

bool RollUpEMAndDeltaRayPolicy::ShouldRollUpCaloHit(
    [[maybe_unused]] const CaloHit *const pCaloHit, [[maybe_unused]] const MCParticle *const pRolledUpMainMC) const
{
    return false;
}

RollUpEMAndAmbiguousDeltaRayHitsPolicy::RollUpEMAndAmbiguousDeltaRayHitsPolicy(
    const float deltaRayParentWeightThreshold, const std::map<HitType, float> deltaRayLengthThresholds) :
    m_deltaRayParentWeightThreshold{deltaRayParentWeightThreshold}
{
    for (const auto &[view, threshold] : deltaRayLengthThresholds)
        m_deltaRayLengthThresholdsSquared[view] = threshold * threshold;
}

bool RollUpEMAndAmbiguousDeltaRayHitsPolicy::ShouldRollUpMC(const MCParticle *const pMC) const
{
    if (pMC->IsRootParticle() || !LArMCParticleHelper::IsEM(pMC))
        return false;

    const MCParticle *const pParentMC{LArMCParticleHelper::GetNextParentMCParticle(pMC)};

    if (LArMCParticleHelper::IsEM(pParentMC))
        return true;

    return false;
}

bool RollUpEMAndAmbiguousDeltaRayHitsPolicy::ShouldRollUpCaloHit(
    const CaloHit *const pCaloHit, const MCParticle *const pRolledUpMainMC) const
{
    // Can't be a delta ray? -> don't roll-up this hit
    if (pRolledUpMainMC->IsRootParticle() || pRolledUpMainMC->GetParticleId() != E_MINUS)
        return false;

    const MCParticle *const pParentMC{LArMCParticleHelper::GetNextParentMCParticle(pRolledUpMainMC)};

    // Is not a delta ray? -> don't roll-up this hit
    const bool parentIsCharged{PdgTable::GetParticleCharge(pParentMC->GetParticleId()) != 0};
    if (LArMCParticleHelper::IsEM(pParentMC) || !parentIsCharged)
        return false;

    // Delta ray starts a shower and is short? -> roll-up this hit
    const float lengthThresholdSquared{m_deltaRayLengthThresholdsSquared.at(pCaloHit->GetHitType())};
    if (!this->CausesShower(pRolledUpMainMC, 0) &&
        (pRolledUpMainMC->GetVertex() - pRolledUpMainMC->GetEndpoint()).GetMagnitudeSquared() < lengthThresholdSquared)
        return true;

    // Delta ray hit that has a significant contribution from the parent -> roll-up this hit
    float parentWeight{std::numeric_limits<float>::lowest()};
    const MCParticleWeightMap &weightMap{pCaloHit->GetMCParticleWeightMap()};
    for (const auto &[pContributingMC, weight] : weightMap)
    {
        if (pContributingMC == pParentMC)
        {
            parentWeight = weight;
            break;
        }
    }
    if (parentWeight > m_deltaRayParentWeightThreshold)
        return true;

    return false;
}

bool RollUpEMAndAmbiguousDeltaRayHitsPolicy::CausesShower(const MCParticle *const pMC, int nDescendantElectrons) const
{
    if (nDescendantElectrons > 1)
        return true;

    if (std::abs(pMC->GetParticleId()) == E_MINUS)
        nDescendantElectrons++; // Including the parent particle, ie. the first in the recursion, as a descendent

    for (const MCParticle *pChildMC : pMC->GetDaughterList())
    {
        if (this->CausesShower(pChildMC, nDescendantElectrons))
            return true;
    }

    return false;
}

} // namespace lar_content
