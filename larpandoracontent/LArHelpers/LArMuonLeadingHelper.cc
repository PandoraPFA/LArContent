/**
 *  @file   larpandoracontent/LArHelpers/LArMuonLeadingHelper.cc
 *
 *  @brief  Implementation of the lar delta ray helper class.
 *
 *  $Log: $
 */

#include "Helpers/MCParticleHelper.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMuonLeadingHelper.h"

#include "Pandora/PdgTable.h"

#include "Objects/CaloHit.h"
#include "Objects/ParticleFlowObject.h"

namespace lar_content
{

using namespace pandora;

LArMuonLeadingHelper::ValidationParameters::ValidationParameters() :
    LArMCParticleHelper::PrimaryParameters(), m_maxBremsstrahlungSeparation(2.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

MCProcess LArMuonLeadingHelper::GetLeadingProcess(const MCParticle *const pMCParticle)
{
    if (LArMCParticleHelper::GetHierarchyTier(pMCParticle) == 0)
        return MC_PROC_UNKNOWN;

    const MCParticle *const pLeadingParticle(LArMuonLeadingHelper::GetLeadingParticle(pMCParticle));
    const LArMCParticle *const pLArMCParticle(dynamic_cast<const LArMCParticle *>(pLeadingParticle));

    if (!pLArMCParticle)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    return pLArMCParticle->GetProcess();
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMuonLeadingHelper::IsDeltaRay(const MCParticle *const pMCParticle)
{
    return (LArMuonLeadingHelper::GetLeadingProcess(pMCParticle) == MC_PROC_MU_IONI);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMuonLeadingHelper::IsMichel(const MCParticle *const pMCParticle)
{
    return (LArMuonLeadingHelper::GetLeadingProcess(pMCParticle) == MC_PROC_DECAY);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMuonLeadingHelper::IsMuonLeading(const MCParticle *const pMCParticle)
{
    const MCParticle *const pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pMCParticle));

    return ((LArMCParticleHelper::GetHierarchyTier(pMCParticle) == 1) && (LArMCParticleHelper::IsCosmicRay(pParentMCParticle)));
}

//------------------------------------------------------------------------------------------------------------------------------------------

const MCParticle *LArMuonLeadingHelper::GetLeadingParticle(const MCParticle *const pMCParticle)
{
    if (!pMCParticle)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    if (LArMCParticleHelper::GetHierarchyTier(pMCParticle) == 0)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const MCParticle *pParentMCParticle = pMCParticle;

    while (LArMCParticleHelper::GetHierarchyTier(pParentMCParticle) != 1)
    {
        const MCParticleList &parentList(pParentMCParticle->GetParentList());

        if (parentList.size() != 1)
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        pParentMCParticle = parentList.front();
    }

    return pParentMCParticle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMuonLeadingHelper::GetMCToLeadingMap(const MCParticleList *const pMCParticleList, LArMCParticleHelper::MCRelationMap &mcToLeadingMap)
{
    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        const MCParticle *const pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pMCParticle));

        if (!LArMCParticleHelper::IsCosmicRay(pParentMCParticle))
            continue;

        // For the CRs: fold hits to themselves, for the DRs: fold hits to the leading MCParticle
        if (pMCParticle == pParentMCParticle)
        {
            mcToLeadingMap[pMCParticle] = pMCParticle;
        }
        else
        {
            const MCParticle *const pLeadingMCParticle(LArMuonLeadingHelper::GetLeadingParticle(pMCParticle));
            mcToLeadingMap[pMCParticle] = pLeadingMCParticle;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMuonLeadingHelper::SelectReconstructableLeadingParticles(const MCParticleList *pMCParticleList, const CaloHitList *pCaloHitList,
    const ValidationParameters &parameters, const CaloHitList &recoMuonHitList, LArMCParticleHelper::MCContributionMap &selectedMCParticlesToHitsMap)
{
    // Obtain hierarchy folding map:
    LArMCParticleHelper::MCRelationMap mcToLeadingMCMap;
    LArMuonLeadingHelper::GetMCToLeadingMap(pMCParticleList, mcToLeadingMCMap);

    // Select reconstructable hits, e.g. remove delta ray hits 'stolen' by the cosmic rays
    CaloHitList selectedCaloHitList;
    LeadingMCParticleToPostBremsstrahlungHitList leadingMCParticleToPostBremsstrahlungHitList;
    LArMuonLeadingHelper::SelectCaloHits(pCaloHitList, mcToLeadingMCMap, selectedCaloHitList, parameters.m_selectInputHits,
        parameters.m_minHitSharingFraction, recoMuonHitList, leadingMCParticleToPostBremsstrahlungHitList);

    // Obtain maps: [hit -> leading MCParticle], [leading MCParticle -> list of hits]
    LArMCParticleHelper::CaloHitToMCMap trueHitToLeadingMCMap;
    LArMCParticleHelper::MCContributionMap leadingMCToTrueHitListMap;
    LArMCParticleHelper::GetMCParticleToCaloHitMatches(&selectedCaloHitList, mcToLeadingMCMap, trueHitToLeadingMCMap, leadingMCToTrueHitListMap);

    // Add in close post bremsstrahlung hits
    LArMuonLeadingHelper::AddInPostBremsstrahlungHits(
        leadingMCParticleToPostBremsstrahlungHitList, parameters.m_maxBremsstrahlungSeparation, leadingMCToTrueHitListMap);

    // Obtain vector: all mc particles
    MCParticleVector leadingMCVector;
    LArMuonLeadingHelper::SelectLeadingMCParticles(pMCParticleList, leadingMCVector);

    // Ensure the MCParticles have enough "good" hits to be reconstructed
    LArMCParticleHelper::SelectParticlesByHitCount(leadingMCVector, leadingMCToTrueHitListMap, mcToLeadingMCMap, parameters, selectedMCParticlesToHitsMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMuonLeadingHelper::SelectCaloHits(const CaloHitList *const pCaloHitList, const LArMCParticleHelper::MCRelationMap &mcToTargetMCMap,
    CaloHitList &selectedCaloHitList, const bool selectInputHits, const float minHitSharingFraction, const CaloHitList &recoMuonHitList,
    LeadingMCParticleToPostBremsstrahlungHitList &leadingMCParticleToPostBremsstrahlungHitList)
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

            if (mcToTargetMCMap.find(pHitParticle) == mcToTargetMCMap.end())
                continue;

            // Remove delta ray hits that have been 'stolen' by the muon
            if (!LArMCParticleHelper::IsCosmicRay(mcToTargetMCMap.at(pHitParticle)))
            {
                if (std::find(recoMuonHitList.begin(), recoMuonHitList.end(), pCaloHit) != recoMuonHitList.end())
                    continue;
            }

            MCParticleVector mcParticleContributionVector;
            for (const auto &mapEntry : pCaloHit->GetMCParticleWeightMap())
                mcParticleContributionVector.push_back(mapEntry.first);

            std::sort(mcParticleContributionVector.begin(), mcParticleContributionVector.end(), PointerLessThan<MCParticle>());

            MCParticleWeightMap targetWeightMap;
            for (const MCParticle *const pMCParticle : mcParticleContributionVector)
            {
                const float weight(pCaloHit->GetMCParticleWeightMap().at(pMCParticle));
                LArMCParticleHelper::MCRelationMap::const_iterator mcIter = mcToTargetMCMap.find(pMCParticle);

                if (mcToTargetMCMap.end() != mcIter)
                    targetWeightMap[mcIter->second] += weight;
            }

            MCParticleVector mcTargetContributionVector;
            for (const auto &mapEntry : targetWeightMap)
                mcTargetContributionVector.push_back(mapEntry.first);
            std::sort(mcTargetContributionVector.begin(), mcTargetContributionVector.end(), PointerLessThan<MCParticle>());

            float bestTargetWeight(0.f), targetWeightSum(0.f);

            for (const MCParticle *const pTargetMCParticle : mcTargetContributionVector)
            {
                const float targetWeight(targetWeightMap.at(pTargetMCParticle));
                targetWeightSum += targetWeight;

                if (targetWeight > bestTargetWeight)
                {
                    bestTargetWeight = targetWeight;
                }
            }

            if ((targetWeightSum < std::numeric_limits<float>::epsilon()) || ((bestTargetWeight / targetWeightSum) < minHitSharingFraction))
                continue;

            // Remove and record post bremsstrahlung hits
            if (LArMuonLeadingHelper::RejectBremsstrahlungHits(pCaloHit, leadingMCParticleToPostBremsstrahlungHitList))
                continue;

            selectedCaloHitList.push_back(pCaloHit);
        }
        catch (const StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMuonLeadingHelper::RejectBremsstrahlungHits(
    const CaloHit *const pCaloHit, LeadingMCParticleToPostBremsstrahlungHitList &leadingMCParticleToPostBremsstrahlungHitList)
{
    const MCParticle *const pHitMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));

    MCParticleList ancestorMCParticleList;
    LArMCParticleHelper::GetAllAncestorMCParticles(pHitMCParticle, ancestorMCParticleList);

    bool isPostBremsstrahlung(false);
    const MCParticle *leadingMCParticle(nullptr);

    for (const MCParticle *const pAncestorMCParticle : ancestorMCParticleList)
    {
        if (LArMCParticleHelper::GetHierarchyTier(pAncestorMCParticle) == 1)
        {
            if (LArMuonLeadingHelper::IsMuonLeading(pAncestorMCParticle))
                leadingMCParticle = pAncestorMCParticle;
        }

        if (pAncestorMCParticle->GetParticleId() == PHOTON)
            isPostBremsstrahlung = true;
    }

    if (isPostBremsstrahlung && leadingMCParticle)
    {
        leadingMCParticleToPostBremsstrahlungHitList[leadingMCParticle].push_back(pCaloHit);
        return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMuonLeadingHelper::AddInPostBremsstrahlungHits(const LeadingMCParticleToPostBremsstrahlungHitList &leadingMCParticleToPostBremsstrahlungHitList,
    const float maxBremsstrahlungSeparation, LArMCParticleHelper::MCContributionMap &leadingMCToTrueHitListMap)
{
    MCParticleVector leadingMCParticleVector;
    for (auto &entry : leadingMCParticleToPostBremsstrahlungHitList)
        leadingMCParticleVector.push_back(entry.first);

    for (const MCParticle *const pLeadingMCParticle : leadingMCParticleVector)
    {
        // Do not add in hits for which there are no main particle hits
        if (leadingMCToTrueHitListMap.find(pLeadingMCParticle) == leadingMCToTrueHitListMap.end())
            continue;

        LArMuonLeadingHelper::AddInPostBremsstrahlungHits(pLeadingMCParticle, leadingMCParticleToPostBremsstrahlungHitList,
            maxBremsstrahlungSeparation, leadingMCToTrueHitListMap, TPC_VIEW_U);
        LArMuonLeadingHelper::AddInPostBremsstrahlungHits(pLeadingMCParticle, leadingMCParticleToPostBremsstrahlungHitList,
            maxBremsstrahlungSeparation, leadingMCToTrueHitListMap, TPC_VIEW_V);
        LArMuonLeadingHelper::AddInPostBremsstrahlungHits(pLeadingMCParticle, leadingMCParticleToPostBremsstrahlungHitList,
            maxBremsstrahlungSeparation, leadingMCToTrueHitListMap, TPC_VIEW_W);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMuonLeadingHelper::AddInPostBremsstrahlungHits(const MCParticle *const pLeadingMCParticle,
    const LeadingMCParticleToPostBremsstrahlungHitList &leadingMCParticleToPostBremsstrahlungHitList,
    const float maxBremsstrahlungSeparation, LArMCParticleHelper::MCContributionMap &leadingMCToTrueHitListMap, const HitType tpcView)
{
    CaloHitList leadingViewHitList;
    for (const CaloHit *const pCaloHit : leadingMCToTrueHitListMap.at(pLeadingMCParticle))
    {
        if (pCaloHit->GetHitType() == tpcView)
            leadingViewHitList.push_back(pCaloHit);
    }

    if (leadingViewHitList.empty())
        return;

    CaloHitList postBremsstrahlungViewHitList;
    for (const CaloHit *const pCaloHit : leadingMCParticleToPostBremsstrahlungHitList.at(pLeadingMCParticle))
    {
        if (pCaloHit->GetHitType() == tpcView)
            postBremsstrahlungViewHitList.push_back(pCaloHit);
    }

    if (postBremsstrahlungViewHitList.empty())
        return;

    bool hitsAdded(false);

    do
    {
        hitsAdded = false;

        for (const CaloHit *const pPostBremsstrahlungHit : postBremsstrahlungViewHitList)
        {
            if (std::find(leadingViewHitList.begin(), leadingViewHitList.end(), pPostBremsstrahlungHit) != leadingViewHitList.end())
                continue;

            const float separationDistance(LArClusterHelper::GetClosestDistance(pPostBremsstrahlungHit->GetPositionVector(), leadingViewHitList));

            if (separationDistance < maxBremsstrahlungSeparation)
            {
                leadingViewHitList.push_back(pPostBremsstrahlungHit);
                hitsAdded = true;
                break;
            }
        }
    } while (hitsAdded);

    CaloHitList &leadingHitList(leadingMCToTrueHitListMap.at(pLeadingMCParticle));
    for (const CaloHit *const pCaloHit : leadingViewHitList)
    {
        if (std::find(leadingHitList.begin(), leadingHitList.end(), pCaloHit) == leadingHitList.end())
            leadingHitList.push_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMuonLeadingHelper::SelectLeadingMCParticles(const MCParticleList *pMCParticleList, MCParticleVector &selectedParticles)
{
    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        const MCParticle *const pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pMCParticle));

        if (!LArMCParticleHelper::IsCosmicRay(pParentMCParticle))
            continue;

        if (pMCParticle == pParentMCParticle)
        {
            selectedParticles.push_back(pMCParticle);
        }
        else
        {
            if (LArMuonLeadingHelper::IsMuonLeading(pMCParticle))
                selectedParticles.push_back(pMCParticle);
        }
    }

    std::sort(selectedParticles.begin(), selectedParticles.end(), LArMCParticleHelper::SortByMomentum);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMuonLeadingHelper::GetPfoMatchContamination(const MCParticle *const pLeadingParticle, const CaloHitList &matchedPfoHitList,
    CaloHitList &parentTrackHits, CaloHitList &otherTrackHits, CaloHitList &otherShowerHits)
{
    const MCParticle *const pParentCosmicRay(LArMCParticleHelper::GetParentMCParticle(pLeadingParticle));

    for (const CaloHit *const pCaloHit : matchedPfoHitList)
    {
        const MCParticle *const pHitParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));

        if (LArMCParticleHelper::IsCosmicRay(pHitParticle))
        {
            (pHitParticle == pParentCosmicRay) ? parentTrackHits.push_back(pCaloHit) : otherTrackHits.push_back(pCaloHit);
        }
        else
        {
            const MCParticle *const pHitLeadingParticle(LArMuonLeadingHelper::GetLeadingParticle(pHitParticle));

            if (pHitLeadingParticle != pLeadingParticle)
                otherShowerHits.push_back(pCaloHit);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMuonLeadingHelper::GetMuonPfoContaminationContribution(
    const CaloHitList &cosmicRayPfoHitList, const CaloHitList &leadingMCHitList, CaloHitList &leadingHitsInParentCosmicRay)
{
    for (const CaloHit *const pCaloHit : cosmicRayPfoHitList)
    {
        if (std::find(leadingMCHitList.begin(), leadingMCHitList.end(), pCaloHit) != leadingMCHitList.end())
            leadingHitsInParentCosmicRay.push_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArMuonLeadingHelper::GetClosestDistance(const Cluster *const pCluster, const CartesianPointVector &cartesianPointVector)
{
    float closestDistance(std::numeric_limits<float>::max());

    CaloHitList caloHitList;
    pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const float distance(LArMuonLeadingHelper::GetClosestDistance(pCaloHit, cartesianPointVector));

        if (distance < closestDistance)
            closestDistance = distance;
    }

    return closestDistance;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArMuonLeadingHelper::GetClosestDistance(const CaloHit *const pCaloHit, const CartesianPointVector &cartesianPointVector)
{
    float shortestDistanceSquared(std::numeric_limits<float>::max());
    const CartesianVector referencePoint(pCaloHit->GetPositionVector());

    for (const CartesianVector &testPosition : cartesianPointVector)
    {
        const float separationSquared((testPosition - referencePoint).GetMagnitudeSquared());

        if (separationSquared < shortestDistanceSquared)
            shortestDistanceSquared = separationSquared;
    }

    return std::sqrt(shortestDistanceSquared);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LArMuonLeadingHelper::GetClosestPosition(const CartesianVector &referencePoint, const CartesianPointVector &cartesianPointVector,
    const Cluster *const pCluster, const float maxDistanceToCluster, const float maxDistanceToReferencePoint, CartesianVector &closestPosition)
{
    bool found(false);
    float shortestDistanceSquared(std::numeric_limits<float>::max());

    for (const CartesianVector &testPosition : cartesianPointVector)
    {
        if (LArClusterHelper::GetClosestDistance(testPosition, pCluster) > maxDistanceToCluster)
            continue;

        const float separationSquared((testPosition - referencePoint).GetMagnitude());

        if (separationSquared > maxDistanceToReferencePoint)
            continue;

        if (separationSquared < shortestDistanceSquared)
        {
            shortestDistanceSquared = separationSquared;
            closestPosition = testPosition;
            found = true;
        }
    }

    return found ? STATUS_CODE_SUCCESS : STATUS_CODE_NOT_FOUND;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMuonLeadingHelper::GetClosestPositions(const CartesianPointVector &cartesianPointVector1, const Cluster *const pCluster2,
    CartesianVector &outputPosition1, CartesianVector &outputPosition2)
{
    bool distanceFound(false);
    float minDistanceSquared(std::numeric_limits<float>::max());

    CartesianVector closestPosition1(0.f, 0.f, 0.f);
    CartesianVector closestPosition2(0.f, 0.f, 0.f);

    CaloHitList caloHitList2;
    pCluster2->GetOrderedCaloHitList().FillCaloHitList(caloHitList2);

    for (const CartesianVector &positionVector1 : cartesianPointVector1)
    {
        for (const CaloHit *const pCaloHit : caloHitList2)
        {
            const CartesianVector &positionVector2(pCaloHit->GetPositionVector());

            const float distanceSquared((positionVector1 - positionVector2).GetMagnitudeSquared());

            if (distanceSquared < minDistanceSquared)
            {
                minDistanceSquared = distanceSquared;
                closestPosition1 = positionVector1;
                closestPosition2 = positionVector2;
                distanceFound = true;
            }
        }
    }

    if (!distanceFound)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    outputPosition1 = closestPosition1;
    outputPosition2 = closestPosition2;
}

} // namespace lar_content
