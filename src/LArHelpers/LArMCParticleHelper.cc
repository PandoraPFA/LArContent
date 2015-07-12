/**
 *  @file   LArContent/src/LArHelpers/LArMCParticleHelper.cc
 *
 *  @brief  Implementation of the lar monte carlo particle helper class.
 *
 *  $Log: $
 */

#include "Objects/MCParticle.h"
#include "Objects/CaloHit.h"
#include "Objects/Cluster.h"

#include "Pandora/PdgTable.h"

#include "LArHelpers/LArMCParticleHelper.h"

#include <cstdlib>

namespace lar_content
{

using namespace pandora;

bool LArMCParticleHelper::IsNeutrinoFinalState(const MCParticle *const pMCParticle)
{
    return ((pMCParticle->GetParentList().size() == 1) && (LArMCParticleHelper::IsNeutrino(*(pMCParticle->GetParentList().begin()))));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsNeutrinoInduced(const MCParticle *const pMCParticle)
{
    return (LArMCParticleHelper::GetParentNeutrinoId(pMCParticle) != 0);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsNeutrino(const MCParticle *const pMCParticle)
{
    const int absoluteParticleId(std::abs(pMCParticle->GetParticleId()));

    if ((NU_E == absoluteParticleId) || (NU_MU == absoluteParticleId) || (NU_TAU == absoluteParticleId))
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsVisible(const MCParticle *const pMCParticle)
{
    const int absoluteParticleId(std::abs(pMCParticle->GetParticleId()));

    if ((E_MINUS == absoluteParticleId) || (MU_MINUS == absoluteParticleId) ||
        (PI_PLUS == absoluteParticleId) || (K_PLUS == absoluteParticleId) ||
        (SIGMA_MINUS == absoluteParticleId) || (SIGMA_PLUS == absoluteParticleId) || (HYPERON_MINUS == absoluteParticleId) ||
        (PROTON == absoluteParticleId) || (PHOTON == absoluteParticleId))
        return true;

    // TODO: What about ions or neutrons?

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const MCParticle *LArMCParticleHelper::GetParentMCParticle(const MCParticle *const pMCParticle)
{
    const MCParticle *pParentMCParticle = pMCParticle;

    while (pParentMCParticle->GetParentList().empty() == false)
    {
        pParentMCParticle = *(pParentMCParticle->GetParentList().begin());
    }

    return pParentMCParticle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const MCParticle *LArMCParticleHelper::GetPrimaryMCParticle(const MCParticle *const pMCParticle)
{
    // Navigate upward through MC daughter/parent links - collect this particle and all its parents
    MCParticleVector mcVector;

    const MCParticle *pParentMCParticle = pMCParticle;
    mcVector.push_back(pParentMCParticle);

    while (!pParentMCParticle->GetParentList().empty())
    {
        pParentMCParticle = *(pParentMCParticle->GetParentList().begin());
        mcVector.push_back(pParentMCParticle);
    }

    // Navigate downward through MC parent/daughter links - return the first long-lived charged particle
    for (MCParticleVector::const_reverse_iterator iter = mcVector.rbegin(), iterEnd = mcVector.rend(); iter != iterEnd; ++iter)
    {
        const MCParticle *const pNextParticle = *iter;

        if (LArMCParticleHelper::IsVisible(pNextParticle))
            return pNextParticle;
    }

    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const MCParticle *LArMCParticleHelper::GetParentNeutrino(const MCParticle *const pMCParticle)
{
    const MCParticle *const pParentMCParticle = LArMCParticleHelper::GetParentMCParticle(pMCParticle);

    if(!LArMCParticleHelper::IsNeutrino(pParentMCParticle))
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return pParentMCParticle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

int LArMCParticleHelper::GetParentNeutrinoId(const MCParticle *const pMCParticle)
{
    try
    {
        const MCParticle *const pParentMCParticle = LArMCParticleHelper::GetParentNeutrino(pMCParticle);
        return pParentMCParticle->GetParticleId();
    }
    catch (const StatusCodeException &)
    {
        return 0;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsNeutrinoInduced(const Cluster *const pCluster, const float minWeight)
{
    return (LArMCParticleHelper::GetNeutrinoWeight(pCluster) > minWeight);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsNeutrinoInduced(const CaloHit *const pCaloHit, const float minWeight)
{
    return (LArMCParticleHelper::GetNeutrinoWeight(pCaloHit) > minWeight);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArMCParticleHelper::GetNeutrinoWeight(const Cluster *const pCluster)
{
    float neutrinoWeight(0.f);
    float totalWeight(0.f);

    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            const CaloHit *const pCaloHit = *hitIter;

            try
            {
                // note: order is important here
                neutrinoWeight += LArMCParticleHelper::GetNeutrinoWeight(pCaloHit);
                totalWeight += 1.f;
            }
            catch (const StatusCodeException &)
            {
            }
        }
    }

    if (totalWeight > std::numeric_limits<float>::epsilon())
        return (neutrinoWeight/totalWeight);

    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArMCParticleHelper::GetNeutrinoWeight(const CaloHit *const pCaloHit)
{
    const MCParticleWeightMap &weights = pCaloHit->GetMCParticleWeightMap();

    if (weights.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    float neutrinoWeight(0.f);
    float totalWeight(0.f);

    for (MCParticleWeightMap::const_iterator iter = weights.begin(), iterEnd = weights.end(); iter != iterEnd; ++iter)
    {
        const MCParticle *const pMCParticle = iter->first;
        const float weight = iter->second;

        if (LArMCParticleHelper::IsNeutrinoInduced(pMCParticle))
            neutrinoWeight += weight;

        totalWeight += weight;
    }

    if (totalWeight > std::numeric_limits<float>::epsilon())
        return (neutrinoWeight/totalWeight);

    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsPrimary(const pandora::MCParticle *const pMCParticle)
{
    try
    {
        const MCParticle *const pPrimaryMCParticle = LArMCParticleHelper::GetPrimaryMCParticle(pMCParticle);
        return (pPrimaryMCParticle == pMCParticle);
    }
    catch (const StatusCodeException &)
    {
        return 0;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetMCPrimaryMap(const MCParticleList *const pMCParticleList, MCRelationMap &mcPrimaryMap)
{
    for (MCParticleList::const_iterator iter = pMCParticleList->begin(), iterEnd = pMCParticleList->end(); iter != iterEnd; ++iter)
    {
        const MCParticle *const pMCParticle = *iter;

        try
        {
            const MCParticle *const pPrimaryMCParticle = LArMCParticleHelper::GetPrimaryMCParticle(pMCParticle);
            mcPrimaryMap[pMCParticle] = pPrimaryMCParticle;
        }
        catch (const StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetPrimaryMCParticleList(const MCParticleList *const pMCParticleList, MCParticleVector &mcPrimaryVector)
{
    for (MCParticleList::const_iterator iter = pMCParticleList->begin(), iterEnd = pMCParticleList->end(); iter != iterEnd; ++iter)
    {
        const MCParticle *const pMCParticle = *iter;

        if (LArMCParticleHelper::IsPrimary(pMCParticle))
            mcPrimaryVector.push_back(pMCParticle);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetNeutrinoMCParticleList(const MCParticleList *const pMCParticleList, MCParticleVector &mcNeutrinoVector)
{
    for (MCParticleList::const_iterator iter = pMCParticleList->begin(), iterEnd = pMCParticleList->end(); iter != iterEnd; ++iter)
    {
        const MCParticle *const pMCParticle = *iter;

        if (LArMCParticleHelper::IsNeutrino(pMCParticle) && pMCParticle->GetParentList().empty())
            mcNeutrinoVector.push_back(pMCParticle);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::SortBySource(const MCParticle *const pLhs, const MCParticle *const pRhs)
{
    // Put neutrino-induced particles first
    const int parentLhs(LArMCParticleHelper::GetParentNeutrinoId(pLhs));
    const int parentRhs(LArMCParticleHelper::GetParentNeutrinoId(pRhs));

    if (parentLhs != parentRhs)
        return (parentLhs > parentRhs);

    return LArMCParticleHelper::SortByMomentum(pLhs, pRhs);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::SortByMomentum(const MCParticle *const pLhs, const MCParticle *const pRhs)
{
    // Sort by momentum (prefer higher momentum)
    const float momentumLhs(pLhs->GetMomentum().GetMagnitudeSquared());
    const float momentumRhs(pRhs->GetMomentum().GetMagnitudeSquared());

    if (momentumLhs != momentumRhs)
        return (momentumLhs > momentumRhs);

    // Sort by energy (prefer lighter particles)
    if (pLhs->GetEnergy() != pRhs->GetEnergy())
        return (pLhs->GetEnergy() < pRhs->GetEnergy());

    // Sort by PDG code (prefer smaller numbers)
    if (pLhs->GetParticleId() != pRhs->GetParticleId())
        return (pLhs->GetParticleId() < pRhs->GetParticleId());

    // Sort by vertex position (tie-breaker)
    const float positionLhs(pLhs->GetVertex().GetMagnitudeSquared());
    const float positionRhs(pRhs->GetVertex().GetMagnitudeSquared());

    return (positionLhs < positionRhs);
}

} // namespace lar_content
