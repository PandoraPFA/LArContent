/**
 *  @file   larpandoracontent/LArHelpers/LArMCParticleHelper.cc
 *
 *  @brief  Implementation of the lar monte carlo particle helper class.
 *
 *  $Log: $
 */

#include "Helpers/MCParticleHelper.h"

#include "Objects/CaloHit.h"
#include "Objects/Cluster.h"
#include "Objects/MCParticle.h"

#include "Pandora/PdgTable.h"
#include "Pandora/StatusCodes.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include <algorithm>
#include <cstdlib>

namespace lar_content
{

using namespace pandora;

LArMCParticleHelper::ValidationParameters::ValidationParameters() :
    m_minPrimaryGoodHits(15),
    m_minHitsForGoodView(5),
    m_minPrimaryGoodViews(2),
    m_selectInputHits(true),
    m_maxPhotonPropagation(2.5f),
    m_minHitSharingFraction(0.9f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsNeutrinoFinalState(const MCParticle *const pMCParticle)
{
    return (LArMCParticleHelper::IsNeutrinoInduced(pMCParticle) && LArMCParticleHelper::IsPrimary(pMCParticle));
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
        (PROTON == absoluteParticleId) || (PHOTON == absoluteParticleId) ||
        (NEUTRON == absoluteParticleId))
        return true;

    // TODO: What about ions or neutrons? Neutrons currently included - they are parents of what would otherwise be large numbers of primary photons
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetTrueNeutrinos(const MCParticleList *const pMCParticleList, MCParticleVector &trueNeutrinos)
{
    if (!pMCParticleList)
        return;

    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        if (pMCParticle->GetParentList().empty() && LArMCParticleHelper::IsNeutrino(pMCParticle))
            trueNeutrinos.push_back(pMCParticle);
    }

    std::sort(trueNeutrinos.begin(), trueNeutrinos.end(), LArMCParticleHelper::SortByMomentum);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const MCParticle *LArMCParticleHelper::GetParentMCParticle(const MCParticle *const pMCParticle)
{
    const MCParticle *pParentMCParticle = pMCParticle;

    while (pParentMCParticle->GetParentList().empty() == false)
    {
        if (1 != pParentMCParticle->GetParentList().size())
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

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
        if (1 != pParentMCParticle->GetParentList().size())
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

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

    if (!LArMCParticleHelper::IsNeutrino(pParentMCParticle))
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

bool LArMCParticleHelper::IsNeutrinoInduced(const Cluster *const pCluster, const float minFraction)
{
    return (LArMCParticleHelper::GetNeutrinoFraction(pCluster) > minFraction);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsNeutrinoInduced(const CaloHit *const pCaloHit, const float minFraction)
{
    return (LArMCParticleHelper::GetNeutrinoFraction(pCaloHit) > minFraction);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
float LArMCParticleHelper::GetNeutrinoFraction(const T *const pT)
{
    float neutrinoWeight(0.f), totalWeight(0.f);
    LArMCParticleHelper::GetNeutrinoWeight(pT, neutrinoWeight, totalWeight);

    if (totalWeight > std::numeric_limits<float>::epsilon())
        return (neutrinoWeight / totalWeight);

    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <>
void LArMCParticleHelper::GetNeutrinoWeight(const CaloHit *const pCaloHit, float &neutrinoWeight, float &totalWeight)
{
    neutrinoWeight = 0.f; totalWeight = 0.f;
    const MCParticleWeightMap &hitMCParticleWeightMap(pCaloHit->GetMCParticleWeightMap());

    if (hitMCParticleWeightMap.empty())
        return;

    MCParticleList mcParticleList;
    for (const auto &mapEntry : hitMCParticleWeightMap) mcParticleList.push_back(mapEntry.first);
    mcParticleList.sort(LArMCParticleHelper::SortByMomentum);

    for (const MCParticle *const pMCParticle : mcParticleList)
    {
        const float weight(hitMCParticleWeightMap.at(pMCParticle));

        if (LArMCParticleHelper::IsNeutrinoInduced(pMCParticle))
            neutrinoWeight += weight;

        totalWeight += weight;
    }

    // ATTN normalise arbitrary input weights at this point
    if (totalWeight > std::numeric_limits<float>::epsilon())
    {
        neutrinoWeight *= 1.f / totalWeight;
        totalWeight = 1.f;
    }
    else
    {
        neutrinoWeight = 0.f;
        totalWeight = 0.f;
    }
}

template <>
void LArMCParticleHelper::GetNeutrinoWeight(const Cluster *const pCluster, float &neutrinoWeight, float &totalWeight)
{
    const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));

    if ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    CaloHitList caloHitList;
    pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);
    LArMCParticleHelper::GetNeutrinoWeight(&caloHitList, neutrinoWeight, totalWeight);
}

template <>
void LArMCParticleHelper::GetNeutrinoWeight(const ParticleFlowObject *const pPfo, float &neutrinoWeight, float &totalWeight)
{
    ClusterList twoDClusters;
    LArPfoHelper::GetTwoDClusterList(pPfo, twoDClusters);
    LArMCParticleHelper::GetNeutrinoWeight(&twoDClusters, neutrinoWeight, totalWeight);
}

template <typename T>
void LArMCParticleHelper::GetNeutrinoWeight(const T *const pT, float &neutrinoWeight, float &totalWeight)
{
    neutrinoWeight = 0.f; totalWeight = 0.f;

    for (const auto *const pValueT : *pT)
    {
        float thisNeutrinoWeight = 0.f, thisTotalWeight = 0.f;
        LArMCParticleHelper::GetNeutrinoWeight(pValueT, thisNeutrinoWeight, thisTotalWeight);
        neutrinoWeight += thisNeutrinoWeight;
        totalWeight += thisTotalWeight;
    }
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
    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
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
    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        if (LArMCParticleHelper::IsPrimary(pMCParticle))
            mcPrimaryVector.push_back(pMCParticle);
    }

    std::sort(mcPrimaryVector.begin(), mcPrimaryVector.end(), LArMCParticleHelper::SortByMomentum);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetNeutrinoMCParticleList(const MCParticleList *const pMCParticleList, MCParticleVector &mcNeutrinoVector)
{
    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        if (LArMCParticleHelper::IsNeutrino(pMCParticle) && pMCParticle->GetParentList().empty())
            mcNeutrinoVector.push_back(pMCParticle);
    }

    std::sort(mcNeutrinoVector.begin(), mcNeutrinoVector.end(), LArMCParticleHelper::SortByMomentum);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::SortByMomentum(const MCParticle *const pLhs, const MCParticle *const pRhs)
{
    // Sort by momentum (prefer higher momentum)
    const float momentumLhs(pLhs->GetMomentum().GetMagnitudeSquared());
    const float momentumRhs(pRhs->GetMomentum().GetMagnitudeSquared());

    if (std::fabs(momentumLhs - momentumRhs) > std::numeric_limits<float>::epsilon())
        return (momentumLhs > momentumRhs);

    // Sort by energy (prefer lighter particles)
    if (std::fabs(pLhs->GetEnergy() - pRhs->GetEnergy()) > std::numeric_limits<float>::epsilon())
        return (pLhs->GetEnergy() < pRhs->GetEnergy());

    // Sort by PDG code (prefer smaller numbers)
    if (pLhs->GetParticleId() != pRhs->GetParticleId())
        return (pLhs->GetParticleId() < pRhs->GetParticleId());

    // Sort by vertex position (tie-breaker)
    const float positionLhs(pLhs->GetVertex().GetMagnitudeSquared());
    const float positionRhs(pRhs->GetVertex().GetMagnitudeSquared());

    return (positionLhs < positionRhs);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetMCParticleToCaloHitMatches(const CaloHitList *const pCaloHitList, const MCRelationMap &mcToPrimaryMCMap,
    CaloHitToMCMap &hitToMCMap, MCContributionMap &mcToTrueHitListMap)
{
    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        try
        {
            const MCParticle *const pHitParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
            const MCParticle *pPrimaryParticle(pHitParticle);

            // ATTN Do not map back to mc primaries if mc to primary mc map not provided
            if (!mcToPrimaryMCMap.empty())
            {
                MCRelationMap::const_iterator mcIter = mcToPrimaryMCMap.find(pHitParticle);

                if (mcToPrimaryMCMap.end() == mcIter)
                    continue;

                pPrimaryParticle = mcIter->second;
            }

            mcToTrueHitListMap[pPrimaryParticle].push_back(pCaloHit);
            hitToMCMap[pCaloHit] = pPrimaryParticle;
        }
        catch (StatusCodeException &statusCodeException)
        {
            if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                throw statusCodeException;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::SelectTrueNeutrinos(const MCParticleList *const pAllMCParticleList, MCParticleVector &selectedMCNeutrinoVector)
{
    MCParticleVector allMCNeutrinoVector;
    LArMCParticleHelper::GetNeutrinoMCParticleList(pAllMCParticleList, allMCNeutrinoVector);

    for (const MCParticle *const pMCNeutrino : allMCNeutrinoVector)
    {
        // ATTN Now demand that input mc neutrinos LArMCParticles, with addition of interaction type
        const LArMCParticle *const pLArMCNeutrino(dynamic_cast<const LArMCParticle*>(pMCNeutrino));

        if (pLArMCNeutrino && (0 != pLArMCNeutrino->GetNuanceCode()))
            selectedMCNeutrinoVector.push_back(pMCNeutrino);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::SelectReconstructableMCParticles(const MCParticleList *pMCParticleList, const CaloHitList *pCaloHitList, const ValidationParameters &parameters,
    std::function<bool(const MCParticle *const)> fCriteria, MCContributionMap &selectedMCParticlesToGoodHitsMap)
{
    // Obtain map: [mc particle -> primary mc particle]
    LArMCParticleHelper::MCRelationMap mcToPrimaryMCMap;
    LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcToPrimaryMCMap);

    // Remove non-reconstructable hits, e.g. those downstream of a neutron
    CaloHitList selectedCaloHitList;
    LArMCParticleHelper::SelectCaloHits(pCaloHitList, mcToPrimaryMCMap, selectedCaloHitList, parameters.m_selectInputHits, parameters.m_maxPhotonPropagation);

    // Remove shared hits where target particle deposits below threshold energy fraction
    CaloHitList goodCaloHitList;
    LArMCParticleHelper::SelectGoodCaloHits(&selectedCaloHitList, mcToPrimaryMCMap, goodCaloHitList, parameters.m_selectInputHits, parameters.m_minHitSharingFraction);

    // Obtain maps: [good hit -> primary mc particle], [primary mc particle -> list of good hits]
    CaloHitToMCMap goodHitToPrimaryMCMap;
    MCContributionMap mcToGoodTrueHitListMap;
    LArMCParticleHelper::GetMCParticleToCaloHitMatches(&goodCaloHitList, mcToPrimaryMCMap, goodHitToPrimaryMCMap, mcToGoodTrueHitListMap);

    // Obtain vector: primary mc particles
    MCParticleVector mcPrimaryVector;
    LArMCParticleHelper::GetPrimaryMCParticleList(pMCParticleList, mcPrimaryVector);

    // Select MCParticles matching criteria
    MCParticleVector candidateTargets;
    LArMCParticleHelper::SelectParticlesMatchingCriteria(mcPrimaryVector, fCriteria, candidateTargets);

    // Ensure the MCParticles have enough hits to be reconstructed
    LArMCParticleHelper::SelectParticlesByHitCount(candidateTargets, mcToGoodTrueHitListMap, parameters, selectedMCParticlesToGoodHitsMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsBeamNeutrinoFinalState(const MCParticle *const pMCParticle)
{
    try
    {
        const MCParticle *const pMCNeutrino(LArMCParticleHelper::GetParentNeutrino(pMCParticle));
        const int nuance(LArMCParticleHelper::GetNuanceCode(pMCNeutrino));
        return (LArMCParticleHelper::IsNeutrinoFinalState(pMCParticle) && nuance != 0 && nuance != 2000 && nuance != 3000);
    }
    catch (const StatusCodeException &) {}

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsBeamParticle(const MCParticle *const pMCParticle)
{
    const int nuance(LArMCParticleHelper::GetNuanceCode(pMCParticle));
    return (LArMCParticleHelper::IsPrimary(pMCParticle) && nuance == 2000);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsCosmicRay(const MCParticle *const pMCParticle)
{
    const int nuance(LArMCParticleHelper::GetNuanceCode(pMCParticle));
    return (LArMCParticleHelper::IsPrimary(pMCParticle) && ((nuance == 0 && !LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCParticle)) || nuance == 3000));
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int LArMCParticleHelper::GetNuanceCode(const MCParticle *const pMCParticle)
{
    const LArMCParticle *const pLArMCParticle(dynamic_cast<const LArMCParticle*>(pMCParticle));
    if (!pLArMCParticle)
    {
        std::cout << "LArMCParticleHelper::GetNuanceCode - Error: Can't cast to LArMCParticle" << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);
    }

    return pLArMCParticle->GetNuanceCode();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(const PfoList &pfoList, const MCContributionMap &selectedMCParticleToGoodHitsMap,
    PfoContributionMap &pfoToReconstructable2DHitsMap)
{
    LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(pfoList, MCContributionMapVector({selectedMCParticleToGoodHitsMap}), pfoToReconstructable2DHitsMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(const PfoList &pfoList, const MCContributionMapVector &selectedMCParticleToGoodHitsMaps,
    PfoContributionMap &pfoToReconstructable2DHitsMap)
{
    for (const ParticleFlowObject *const pPfo : pfoList)
    {
        CaloHitList pfoHitList;
        LArMCParticleHelper::CollectReconstructable2DHits(pPfo, selectedMCParticleToGoodHitsMaps, pfoHitList);

        if (!pfoToReconstructable2DHitsMap.insert(PfoContributionMap::value_type(pPfo, pfoHitList)).second)
            throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(const PfoContributionMap &pfoToReconstructable2DHitsMap, const MCContributionMapVector &selectedMCParticleToGoodHitsMaps,
    PfoToMCParticleHitSharingMap &pfoToMCParticleHitSharingMap, MCParticleToPfoHitSharingMap &mcParticleToPfoHitSharingMap)
{
    PfoVector sortedPfos;
    for (const auto &mapEntry : pfoToReconstructable2DHitsMap) sortedPfos.push_back(mapEntry.first);
    std::sort(sortedPfos.begin(), sortedPfos.end(), LArPfoHelper::SortByNHits);

    for (const ParticleFlowObject *const pPfo : sortedPfos)
    {
        for (const MCContributionMap &mcParticleToGoodHitsMap : selectedMCParticleToGoodHitsMaps)
        {
            MCParticleVector sortedMCParticles;
            for (const auto &mapEntry : mcParticleToGoodHitsMap) sortedMCParticles.push_back(mapEntry.first);
            std::sort(sortedMCParticles.begin(), sortedMCParticles.end(), PointerLessThan<MCParticle>());

            for (const MCParticle *const pMCParticle : sortedMCParticles)
            {
                // Add map entries for this Pfo & MCParticle if required
                if (pfoToMCParticleHitSharingMap.find(pPfo) == pfoToMCParticleHitSharingMap.end())
                    if (!pfoToMCParticleHitSharingMap.insert(PfoToMCParticleHitSharingMap::value_type(pPfo, std::vector<MCParticleIntPair>())).second)
                        throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT); // ATTN maybe overkill

                if (mcParticleToPfoHitSharingMap.find(pMCParticle) == mcParticleToPfoHitSharingMap.end())
                    if (!mcParticleToPfoHitSharingMap.insert(MCParticleToPfoHitSharingMap::value_type(pMCParticle, std::vector<PfoIntPair>())).second)
                        throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);

                // Check this Pfo & MCParticle pairing hasn't already been checked
                std::vector<MCParticleIntPair> &mcHitPairs(pfoToMCParticleHitSharingMap.at(pPfo));
                std::vector<PfoIntPair> &pfoHitPairs(mcParticleToPfoHitSharingMap.at(pMCParticle));

                if (std::any_of(mcHitPairs.begin(), mcHitPairs.end(), [&] (const MCParticleIntPair &pair) { return (pair.first == pMCParticle); }))
                    throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);

                if (std::any_of(pfoHitPairs.begin(), pfoHitPairs.end(), [&] (const PfoIntPair &pair) { return (pair.first == pPfo); }))
                    throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);

                // Add the number of shared hits to the maps
                const unsigned int nSharedHits(LArMCParticleHelper::CountSharedHits(pfoToReconstructable2DHitsMap.at(pPfo), mcParticleToGoodHitsMap.at(pMCParticle)));
                mcHitPairs.push_back(MCParticleIntPair(pMCParticle, nSharedHits));
                pfoHitPairs.push_back(PfoIntPair(pPfo, nSharedHits));
            }
        }
    }
}

// private
//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::CollectReconstructable2DHits(const ParticleFlowObject *const pPfo, const MCContributionMapVector &selectedMCParticleToGoodHitsMaps,
    pandora::CaloHitList &reconstructableCaloHitList2D)
{
    // Collect all 2D calo hits
    CaloHitList caloHitList2D;
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_U, caloHitList2D);
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_V, caloHitList2D);
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, caloHitList2D);

    // Filter for only reconstructable hits
    for (const CaloHit *const pCaloHit : caloHitList2D)
    {
        bool isTargetHit(false);
        for (const MCContributionMap &mcParticleToGoodHitsMap : selectedMCParticleToGoodHitsMaps)
        {
            // ATTN This map is unordered, but this does not impact search for specific target hit
            for (const MCContributionMap::value_type &mapEntry : mcParticleToGoodHitsMap)
            {
                if (std::find(mapEntry.second.begin(), mapEntry.second.end(), pCaloHit) != mapEntry.second.end())
                {
                    isTargetHit = true;
                    break;
                }
            }
            if (isTargetHit) break;
        }

        if (isTargetHit)
            reconstructableCaloHitList2D.push_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::SelectCaloHits(const CaloHitList *const pCaloHitList, const LArMCParticleHelper::MCRelationMap &mcToPrimaryMCMap,
    CaloHitList &selectedCaloHitList, const bool selectInputHits, const float maxPhotonPropagation)
{
    if (!selectInputHits)
    {
        selectedCaloHitList.insert(selectedCaloHitList.end(), pCaloHitList->begin(), pCaloHitList->end());
        return;
    }

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        try
        {
            const MCParticle *const pHitParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));

            LArMCParticleHelper::MCRelationMap::const_iterator mcIter = mcToPrimaryMCMap.find(pHitParticle);

            if (mcToPrimaryMCMap.end() == mcIter)
                continue;

            const MCParticle *const pPrimaryParticle = mcIter->second;

            if (PassMCParticleChecks(pPrimaryParticle, pPrimaryParticle, pHitParticle, maxPhotonPropagation))
                selectedCaloHitList.push_back(pCaloHit);
        }
        catch (const StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::SelectGoodCaloHits(const CaloHitList *const pSelectedCaloHitList, const LArMCParticleHelper::MCRelationMap &mcToPrimaryMCMap,
    CaloHitList &selectedGoodCaloHitList, const bool selectInputHits, const float minHitSharingFraction)
{
    if (!selectInputHits)
    {
        selectedGoodCaloHitList.insert(selectedGoodCaloHitList.end(), pSelectedCaloHitList->begin(), pSelectedCaloHitList->end());
        return;
    }

    for (const CaloHit *const pCaloHit : *pSelectedCaloHitList)
    {
        MCParticleVector mcParticleVector;
        for (const auto &mapEntry : pCaloHit->GetMCParticleWeightMap()) mcParticleVector.push_back(mapEntry.first);
        std::sort(mcParticleVector.begin(), mcParticleVector.end(), PointerLessThan<MCParticle>());

        MCParticleWeightMap primaryWeightMap;

        for (const MCParticle *const pMCParticle : mcParticleVector)
        {
            const float weight(pCaloHit->GetMCParticleWeightMap().at(pMCParticle));
            LArMCParticleHelper::MCRelationMap::const_iterator mcIter = mcToPrimaryMCMap.find(pMCParticle);

            if (mcToPrimaryMCMap.end() != mcIter)
                primaryWeightMap[mcIter->second] += weight;
        }

        MCParticleVector mcPrimaryVector;
        for (const auto &mapEntry : primaryWeightMap) mcPrimaryVector.push_back(mapEntry.first);
        std::sort(mcPrimaryVector.begin(), mcPrimaryVector.end(), PointerLessThan<MCParticle>());

        const MCParticle *pBestPrimaryParticle(nullptr);
        float bestPrimaryWeight(0.f), primaryWeightSum(0.f);

        for (const MCParticle *const pPrimaryMCParticle : mcPrimaryVector)
        {
            const float primaryWeight(primaryWeightMap.at(pPrimaryMCParticle));
            primaryWeightSum += primaryWeight;

            if (primaryWeight > bestPrimaryWeight)
            {
                bestPrimaryWeight = primaryWeight;
                pBestPrimaryParticle = pPrimaryMCParticle;
            }
        }

        if (!pBestPrimaryParticle || (primaryWeightSum < std::numeric_limits<float>::epsilon()) || ((bestPrimaryWeight / primaryWeightSum) < minHitSharingFraction))
            continue;

        selectedGoodCaloHitList.push_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::SelectParticlesMatchingCriteria(const MCParticleVector &inputMCParticles, std::function<bool(const MCParticle *const)> fCriteria,
    MCParticleVector &selectedParticles)
{
    for (const MCParticle *const pMCParticle : inputMCParticles)
    {
        if (fCriteria(pMCParticle))
            selectedParticles.push_back(pMCParticle);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::SelectParticlesByHitCount(const MCParticleVector &candidateTargets, const MCContributionMap &mcToGoodTrueHitListMap,
    const ValidationParameters &parameters, MCContributionMap &selectedMCParticlesToGoodHitsMap)
{
    // Apply restrictions on the number of good hits associated with the MCParticles
    for (const MCParticle * const pMCTarget : candidateTargets)
    {
        MCContributionMap::const_iterator goodTrueHitsIter = mcToGoodTrueHitListMap.find(pMCTarget);
        if (mcToGoodTrueHitListMap.end() == goodTrueHitsIter)
            continue;

        const CaloHitList &caloHitList(goodTrueHitsIter->second);
        if (caloHitList.size() < parameters.m_minPrimaryGoodHits)
            continue;

        unsigned int nGoodViews(0);
        if (LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, caloHitList) >= parameters.m_minHitsForGoodView)
            ++nGoodViews;

        if (LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, caloHitList) >= parameters.m_minHitsForGoodView)
            ++nGoodViews;

        if (LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, caloHitList) >= parameters.m_minHitsForGoodView)
            ++nGoodViews;

        if (nGoodViews < parameters.m_minPrimaryGoodViews)
            continue;

        if (!selectedMCParticlesToGoodHitsMap.insert(MCContributionMap::value_type(goodTrueHitsIter->first, goodTrueHitsIter->second)).second)
            throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::PassMCParticleChecks(const MCParticle *const pOriginalPrimary, const MCParticle *const pThisMCParticle,
    const MCParticle *const pHitMCParticle, const float maxPhotonPropagation)
{
    if (NEUTRON == std::abs(pThisMCParticle->GetParticleId()))
        return false;

    if ((PHOTON == pThisMCParticle->GetParticleId()) && (PHOTON != pOriginalPrimary->GetParticleId()) && (E_MINUS != std::abs(pOriginalPrimary->GetParticleId())))
    {
        if ((pThisMCParticle->GetEndpoint() - pThisMCParticle->GetVertex()).GetMagnitude() > maxPhotonPropagation)
            return false;
    }

    if (pThisMCParticle == pHitMCParticle)
        return true;

    for (const MCParticle *const pDaughterMCParticle : pThisMCParticle->GetDaughterList())
    {
        if (PassMCParticleChecks(pOriginalPrimary, pDaughterMCParticle, pHitMCParticle, maxPhotonPropagation))
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int LArMCParticleHelper::CountSharedHits(const CaloHitList &hitListA, const CaloHitList &hitListB)
{
    unsigned int nSharedHits(0);
    for (const CaloHit *const pCaloHit : hitListA)
    {
        if (std::find(hitListB.begin(), hitListB.end(), pCaloHit) != hitListB.end())
            nSharedHits++;
    }

    return nSharedHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template float LArMCParticleHelper::GetNeutrinoFraction(const CaloHit *const);
template float LArMCParticleHelper::GetNeutrinoFraction(const Cluster *const);
template float LArMCParticleHelper::GetNeutrinoFraction(const ParticleFlowObject *const);
template float LArMCParticleHelper::GetNeutrinoFraction(const CaloHitList *const);
template float LArMCParticleHelper::GetNeutrinoFraction(const ClusterList *const);
template float LArMCParticleHelper::GetNeutrinoFraction(const PfoList *const);

template void LArMCParticleHelper::GetNeutrinoWeight(const CaloHitList *const, float &, float &);
template void LArMCParticleHelper::GetNeutrinoWeight(const ClusterList *const, float &, float &);
template void LArMCParticleHelper::GetNeutrinoWeight(const PfoList *const, float &, float &);

} // namespace lar_content
