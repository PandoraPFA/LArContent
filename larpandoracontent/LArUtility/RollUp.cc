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

void RollUpper::Reset()
{
    m_mcCache.clear();
    m_caloHitCache.clear();
}

const MCParticle *RollUpper::RollUpMC(const MCParticle *const pMC)
{
    if (m_mcCache.find(pMC) == m_mcCache.end())
        m_mcCache.insert({pMC, m_policy->GetRollUpTargetMC(pMC)});

    return m_mcCache.at(pMC);
}

const MCParticle *RollUpper::RollUpCaloHit(const CaloHit *const pCaloHit)
{
    if (m_caloHitCache.find(pCaloHit) == m_caloHitCache.end())
    {
        const MCParticleWeightMap &weightMap{pCaloHit->GetMCParticleWeightMap()};
        MCParticleWeightMap rolledUpWeightMap;
        for (const auto &[pMC, weight] : weightMap)
            rolledUpWeightMap[this->RollUpMC(pMC)] += weight;

        const MCParticle *pMainMC{nullptr};
        float maxWeight{0.f};
        for (const auto &[pMC, weight] : rolledUpWeightMap)
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

        if (pMainMC && m_policy->ShouldFoldCaloHit(pCaloHit, pMainMC))
            pMainMC = LArMCParticleHelper::GetNextParentMCParticle(pMainMC);

        m_caloHitCache.insert({pCaloHit, pMainMC});
    }

    return m_caloHitCache.at(pCaloHit);
}

const MCParticle *RollUpEMPolicy::GetRollUpTargetMC(const MCParticle *const pMC) const
{
    if (!LArMCParticleHelper::IsEM(pMC))
        return pMC;

    const MCParticle *pLeadingMC{pMC};
    while (!pLeadingMC->IsRootParticle())
    {
        const MCParticle *const pParentMC{LArMCParticleHelper::GetNextParentMCParticle(pLeadingMC)};
        if (!LArMCParticleHelper::IsEM(pParentMC))
            break;
        pLeadingMC = pParentMC;
    }

    return pLeadingMC;
}

bool RollUpEMPolicy::ShouldFoldCaloHit(
    [[maybe_unused]] const CaloHit *const pCaloHit, [[maybe_unused]] const MCParticle *const pRolledUpMainMC) const
{
    return false;
}

const MCParticle *RollUpEMAndDeltaRayPolicy::GetRollUpTargetMC(const MCParticle *const pMC) const
{
    if (!LArMCParticleHelper::IsEM(pMC))
        return pMC;

    const MCParticle *pLeadingMC{pMC};
    while (!pLeadingMC->IsRootParticle())
    {
        const MCParticle *const pParentMC{LArMCParticleHelper::GetNextParentMCParticle(pLeadingMC)};
        if (!LArMCParticleHelper::IsEM(pParentMC))
        {
            // Is electron and parent is a charged track -> jump up the chain by one to roll-up delta ray
            const bool parentIsCharged{PdgTable::GetParticleCharge(pParentMC->GetParticleId()) != 0};
            if (pLeadingMC->GetParticleId() == E_MINUS && parentIsCharged)
                pLeadingMC = pParentMC;
            break;
        }
        pLeadingMC = pParentMC;
    }

    return pLeadingMC;
}

RollUpEMAndAmbiguousDeltaRayHitsPolicy::RollUpEMAndAmbiguousDeltaRayHitsPolicy(
    const float deltaRayParentWeightThreshold, const std::map<HitType, float> deltaRayLengthThresholds) :
    m_deltaRayParentWeightThreshold{deltaRayParentWeightThreshold}
{
    for (const auto &[view, threshold] : deltaRayLengthThresholds)
        m_deltaRayLengthThresholdsSquared[view] = threshold * threshold;
}

bool RollUpEMAndAmbiguousDeltaRayHitsPolicy::ShouldFoldCaloHit(
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

    // Delta ray doesn't shower and is short? -> roll-up this hit
    const float lengthThresholdSquared{m_deltaRayLengthThresholdsSquared.at(pCaloHit->GetHitType())};
    if (!LArMCParticleHelper::CausesShower(pRolledUpMainMC, 0) &&
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

const MCParticle *RollUpEMWithComptonFilterAndAmbiguousDeltaRayHitsPolicy::GetRollUpTargetMC(const MCParticle *const pMC) const
{
    if (!LArMCParticleHelper::IsEM(pMC))
        return pMC;

    bool hasAncestorElectron{false};
    const MCParticle *pLeadingMC{pMC};
    while (!pLeadingMC->IsRootParticle())
    {
        const MCParticle *const pParentMC{LArMCParticleHelper::GetNextParentMCParticle(pLeadingMC)};
        if (!LArMCParticleHelper::IsEM(pParentMC))
            break;
        if (std::abs(pParentMC->GetParticleId()) == E_MINUS)
            hasAncestorElectron = true;
        pLeadingMC = pParentMC;
    }

    // Don't roll-up chains of EM that consist only of compton scatters, prevents distant diffuse hits sharing the same true mc
    if (!hasAncestorElectron && std::abs(pLeadingMC->GetParticleId()) != E_MINUS && !LArMCParticleHelper::CausesShower(pLeadingMC, 0))
        return pMC;

    return pLeadingMC;
}

} // namespace lar_content
