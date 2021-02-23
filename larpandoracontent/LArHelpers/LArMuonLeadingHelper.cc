/**
 *  @file   larpandoracontent/LArHelpers/LArMuonLeadingHelper.cc
 *
 *  @brief  Implementation of the lar delta ray helper class.
 *
 *  $Log: $
 */

#include "Helpers/MCParticleHelper.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArMuonLeadingHelper.h"

#include "Pandora/PdgTable.h"

#include "PandoraMonitoringApi.h"

#include "Objects/ParticleFlowObject.h"
#include "Objects/CaloHit.h"

namespace lar_content
{

using namespace pandora;

LArMuonLeadingHelper::ValidationParameters::ValidationParameters() :
    LArMCParticleHelper::PrimaryParameters(),
    m_maxBremsstrahlungSeparation(2.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------    

bool LArMuonLeadingHelper::IsDeltaRay(const MCParticle *const pMCParticle)
{
    const MCParticleList parentList(pMCParticle->GetParentList());

    if(parentList.empty())
        return false;

    if (1 != parentList.size())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    if (!LArMCParticleHelper::IsCosmicRay(parentList.front()))
        return false;

    const LArMCParticle *const pLArMCParticle(dynamic_cast<const LArMCParticle*>(pMCParticle));
    
    return pLArMCParticle->GetIsDR();
}

//------------------------------------------------------------------------------------------------------------------------------------------    

bool LArMuonLeadingHelper::IsMichel(const MCParticle *const pMCParticle)
{
    if (std::abs(pMCParticle->GetParticleId()) != 11)
        return false;
    
    const MCParticleList parentList(pMCParticle->GetParentList());

    if(parentList.empty())
        return false;
    
    if (1 != parentList.size())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    
    if (!LArMCParticleHelper::IsCosmicRay(parentList.front()))
        return false;    
    
    const LArMCParticle *const pLArMCParticle(dynamic_cast<const LArMCParticle*>(pMCParticle));

    return pLArMCParticle->GetIsDecay();
}

//------------------------------------------------------------------------------------------------------------------------------------------    

bool LArMuonLeadingHelper::IsLeading(const MCParticle *const pMCParticle)
{
    const MCParticle *const pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pMCParticle));
    
    return ((LArMCParticleHelper::GetHierarchyTier(pMCParticle) == 1) && (LArMCParticleHelper::IsCosmicRay(pParentMCParticle)));
}     

//------------------------------------------------------------------------------------------------------------------------------------------

const MCParticle *LArMuonLeadingHelper::GetLeadingParticle(const MCParticle *const pMCParticle)
{
    // Navigate upward through MC daughter/parent links - collect this particle and all its parents
    MCParticleVector mcParticleVector;

    const MCParticle *pParentMCParticle = pMCParticle;
    mcParticleVector.push_back(pParentMCParticle);

    while (!pParentMCParticle->GetParentList().empty())
    {
        if (1 != pParentMCParticle->GetParentList().size())
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        pParentMCParticle = *(pParentMCParticle->GetParentList().begin());
        mcParticleVector.push_back(pParentMCParticle);
    }

    return *(mcParticleVector.end() - 2);
}    

//------------------------------------------------------------------------------------------------------------------------------------------    

void LArMuonLeadingHelper::GetMCToLeadingMap(const MCParticleList *const pMCParticleList, LArMCParticleHelper::MCRelationMap &mcToLeadingMap)
{
    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        const MCParticle *const pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pMCParticle));

        if (!LArMCParticleHelper::IsCosmicRay(pParentMCParticle))
            continue;

        // For the CRs: fold hits to themselves, for the DRs: fold hits to the leading non-muon MC
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

void LArMuonLeadingHelper::SelectReconstructableLeadingParticles(const MCParticleList *pMCParticleList, const CaloHitList *pCaloHitList, const ValidationParameters &parameters,
    const CaloHitList &recoMuonHitList, LArMCParticleHelper::MCContributionMap &selectedMCParticlesToHitsMap, const Pandora &pandora)
{
    // Obtain map: [mc particle -> leading muon child mc]
    LArMCParticleHelper::MCRelationMap mcToLeadingMCMap;
    LArMuonLeadingHelper::GetMCToLeadingMap(pMCParticleList, mcToLeadingMCMap);

    // Select reconstructable hits, e.g. remove those downstream of a neutron
    // Unless selectInputHits == false
    CaloHitList selectedCaloHitList;
    LeadingMCParticleToPostPhotonHitLists leadingMCParticleToPostPhotonHitLists;
    LArMuonLeadingHelper::SelectCaloHits(pCaloHitList, mcToLeadingMCMap, selectedCaloHitList, parameters.m_selectInputHits, parameters.m_minHitSharingFraction,
        recoMuonHitList, leadingMCParticleToPostPhotonHitLists);

    // Obtain maps: [hit -> leading muon child mc], [leading muon child mc -> list of hits]
    LArMCParticleHelper::CaloHitToMCMap trueHitToLeadingMCMap;
    LArMCParticleHelper::MCContributionMap leadingMCToTrueHitListMap;
    LArMCParticleHelper::GetMCParticleToCaloHitMatches(&selectedCaloHitList, mcToLeadingMCMap, trueHitToLeadingMCMap, leadingMCToTrueHitListMap);

    // Add back in post bremsstrahlung hits that are close to the non-muon leading
    LArMuonLeadingHelper::AddInReconstructablePostPhotonHits(leadingMCParticleToPostPhotonHitLists, parameters.m_maxBremsstrahlungSeparation, leadingMCToTrueHitListMap, pandora);

    // Obtain vector: all mc particles
    MCParticleVector leadingMCVector;
    LArMuonLeadingHelper::SelectLeadingMCParticles(pMCParticleList, leadingMCVector);
    
    // Ensure the MCParticles have enough "good" hits to be reconstructed
    LArMCParticleHelper::SelectParticlesByHitCount(leadingMCVector, leadingMCToTrueHitListMap, mcToLeadingMCMap, parameters, selectedMCParticlesToHitsMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMuonLeadingHelper::SelectCaloHits(const CaloHitList *const pCaloHitList, const LArMCParticleHelper::MCRelationMap &mcToTargetMCMap,
    CaloHitList &selectedCaloHitList, const bool selectInputHits, const float minHitSharingFraction, const CaloHitList &recoMuonHitList,
    LeadingMCParticleToPostPhotonHitLists &leadingMCParticleToPostPhotonHitLists)
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

            // Remove hits from metrics that have been taken by the muon
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
            for (const auto &mapEntry : targetWeightMap) mcTargetContributionVector.push_back(mapEntry.first);
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
            
            if(RejectBremsstrahlungHits(pCaloHit, leadingMCParticleToPostPhotonHitLists))
                continue;

            selectedCaloHitList.push_back(pCaloHit);
        }
        catch (const StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMuonLeadingHelper::RejectBremsstrahlungHits(const CaloHit *const pCaloHit, LeadingMCParticleToPostPhotonHitLists &leadingMCParticleToPostPhotonHitLists)
{
    const MCParticle *const pHitMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
    
    MCParticleList ancestorMCParticleList;
    LArMCParticleHelper::GetAllAncestorMCParticles(pHitMCParticle, ancestorMCParticleList);

    int highestTier(0);
    const MCParticle *leadingMCParticle(nullptr), *highestTierPhoton(nullptr);

    for (const MCParticle *const pAncestorMCParticle : ancestorMCParticleList)
    {
        if (LArMCParticleHelper::GetHierarchyTier(pAncestorMCParticle) == 1)
        {
            if (LArMuonLeadingHelper::IsLeading(pAncestorMCParticle))
                leadingMCParticle = pAncestorMCParticle;
        }

        if (pAncestorMCParticle->GetParticleId() == PHOTON)
        {
            if (LArMCParticleHelper::GetHierarchyTier(pAncestorMCParticle) > highestTier)
            {
                highestTier = LArMCParticleHelper::GetHierarchyTier(pAncestorMCParticle);
                highestTierPhoton = pAncestorMCParticle;
            }
        }
    }
   
    if (leadingMCParticle && highestTierPhoton)
    {
        leadingMCParticleToPostPhotonHitLists[leadingMCParticle][pHitMCParticle].push_back(pCaloHit);        
        return true;
    }
    else
    {
        return false;
    }
}    

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMuonLeadingHelper::AddInReconstructablePostPhotonHits(const LeadingMCParticleToPostPhotonHitLists &leadingMCParticleToPostPhotonHitLists, const float maxBremsstrahlungSeparation,
    LArMCParticleHelper::MCContributionMap &leadingMCToTrueHitListMap, const Pandora &pandora)
{
    MCParticleVector leadingMCParticleVector;
    for (auto &entry : leadingMCParticleToPostPhotonHitLists)
        leadingMCParticleVector.push_back(entry.first);
    std::sort(leadingMCParticleVector.begin(), leadingMCParticleVector.end(), LArMCParticleHelper::SortByMomentum);

    for (const MCParticle *const pLeadingMCParticle : leadingMCParticleVector)
    {
        // DO NOT ADD IN HITS FOR WHICH THERE IS NO MAIN PARTICLE HITS
        if(leadingMCToTrueHitListMap.find(pLeadingMCParticle) == leadingMCToTrueHitListMap.end())
        {
            ////////////////////////////////
            /*
            std::cout << "ISOBEL - CASE WHERE THERE IS NO INITIAL DR TO ADD HITS ON TO" << std::endl;
            for (auto &entry : leadingMCParticleToPostPhotonHitLists.at(pLeadingMCParticle))
            {
                for (const CaloHit *const pCaloHit : entry.second)
                {
                    const CartesianVector &pos(pCaloHit->GetPositionVector());
                    PandoraMonitoringApi::AddMarkerToVisualization(pandora, &pos, "POSITION", BLACK, 2);
                }
            }
            */
            ////////////////////////////////            
            continue;
        }

        LArMuonLeadingHelper::AddHits(pLeadingMCParticle, leadingMCParticleToPostPhotonHitLists, maxBremsstrahlungSeparation, leadingMCToTrueHitListMap, TPC_VIEW_U, pandora);
        LArMuonLeadingHelper::AddHits(pLeadingMCParticle, leadingMCParticleToPostPhotonHitLists, maxBremsstrahlungSeparation, leadingMCToTrueHitListMap, TPC_VIEW_V, pandora);
        LArMuonLeadingHelper::AddHits(pLeadingMCParticle, leadingMCParticleToPostPhotonHitLists, maxBremsstrahlungSeparation, leadingMCToTrueHitListMap, TPC_VIEW_W, pandora);        
    }    
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMuonLeadingHelper::AddHits(const MCParticle *const pLeadingMCParticle, const LeadingMCParticleToPostPhotonHitLists &leadingMCParticleToPostPhotonHitLists,
    const float maxBremsstrahlungSeparation, LArMCParticleHelper::MCContributionMap &leadingMCToTrueHitListMap, const HitType &tpcView, const Pandora &/*pandora*/)
{
    CaloHitList leadingHitList;
    for (const CaloHit *const pCaloHit : leadingMCToTrueHitListMap.at(pLeadingMCParticle))
    {
        if (pCaloHit->GetHitType() == tpcView)
            leadingHitList.push_back(pCaloHit);
    }

    if (leadingHitList.empty())
        return;
    
    CaloHitList postBremsstrahlungHits;
    for (auto &entry : leadingMCParticleToPostPhotonHitLists.at(pLeadingMCParticle))
    {
        for (const CaloHit *const pCaloHit : entry.second)
        {
            if (pCaloHit->GetHitType() == tpcView)
                postBremsstrahlungHits.push_back(pCaloHit);
        }
    }

    if (postBremsstrahlungHits.empty())
        return;

    //////////////////////////////
    /*
    for (const CaloHit *const pCaloHit : leadingHitList)
    {
        const CartesianVector shiftedPosition(pCaloHit->GetPositionVector().GetX() - pCaloHit->GetX0(), pCaloHit->GetPositionVector().GetY(), pCaloHit->GetPositionVector().GetZ());
        PandoraMonitoringApi::AddMarkerToVisualization(pandora, &shiftedPosition, "Main Delta Ray Hit", BLACK, 2);
    }

    PandoraMonitoringApi::ViewEvent(pandora);
    for (const CaloHit *const pCaloHit : postBremsstrahlungHits)
    {
        const CartesianVector shiftedPosition(pCaloHit->GetPositionVector().GetX() - pCaloHit->GetX0(), pCaloHit->GetPositionVector().GetY(), pCaloHit->GetPositionVector().GetZ());
        PandoraMonitoringApi::AddMarkerToVisualization(pandora, &shiftedPosition, "Post Bremsstrahlung Hit", BLUE, 2);
    }
        
    PandoraMonitoringApi::ViewEvent(pandora);
    */
    //////////////////////////////    


    //std::cout << "maxBremsstrahlungSeparation: " << maxBremsstrahlungSeparation << std::endl;
    
    bool hitsAdded(true);
    while (hitsAdded)
    {
        hitsAdded = false;

        for (const CaloHit *const pPostBremsstrahlungHit : postBremsstrahlungHits)
        {
            if (std::find(leadingHitList.begin(), leadingHitList.end(), pPostBremsstrahlungHit) != leadingHitList.end())
                continue;
            
            const float separationDistance(LArClusterHelper::GetClosestDistance(pPostBremsstrahlungHit->GetPositionVector(), leadingHitList));

            if (separationDistance < maxBremsstrahlungSeparation)
            {
                leadingHitList.push_back(pPostBremsstrahlungHit);
                hitsAdded = true;
                break;
            }
        }
    }

    CaloHitList &hits(leadingMCToTrueHitListMap.at(pLeadingMCParticle));
    for (const CaloHit *const pCaloHit : leadingHitList)
    {
        if (std::find(hits.begin(), hits.end(), pCaloHit) == hits.end())
            hits.push_back(pCaloHit);
    }  
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
            if (pHitParticle == pParentCosmicRay)
            {
                parentTrackHits.push_back(pCaloHit);
            }
            else
            {
                otherTrackHits.push_back(pCaloHit);
            }
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

void LArMuonLeadingHelper::GetMuonPfoContaminationContribution(const CaloHitList &cosmicRayPfoHitList, const CaloHitList &leadingMCHitList,
    CaloHitList &leadingHitsInParentCosmicRay)
{
    for (const CaloHit *const pCaloHit : cosmicRayPfoHitList)
    {
        if (std::find(leadingMCHitList.begin(), leadingMCHitList.end(), pCaloHit) != leadingMCHitList.end())
            leadingHitsInParentCosmicRay.push_back(pCaloHit);
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
            if (LArMuonLeadingHelper::IsLeading(pMCParticle))
                selectedParticles.push_back(pMCParticle);
        }
    }

    std::sort(selectedParticles.begin(), selectedParticles.end(), LArMCParticleHelper::SortByMomentum);
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

float LArMuonLeadingHelper::GetClosestDistance(const CaloHit *const pReferenceCaloHit, const CaloHitList &caloHitList)
{
    CartesianPointVector cartesianPointVector;
    
    for (const CaloHit *const pCaloHit : caloHitList)
        cartesianPointVector.push_back(pCaloHit->GetPositionVector());

    return LArMuonLeadingHelper::GetClosestDistance(pReferenceCaloHit, cartesianPointVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMuonLeadingHelper::GetClosestPositions(const CartesianPointVector &pCluster1, const Cluster *const pCluster2, CartesianVector &outputPosition1,
    CartesianVector &outputPosition2)
{
    bool distanceFound(false);
    float minDistanceSquared(std::numeric_limits<float>::max());

    CartesianVector closestPosition1(0.f, 0.f, 0.f);
    CartesianVector closestPosition2(0.f, 0.f, 0.f);

    CaloHitList caloHitList2;
    pCluster2->GetOrderedCaloHitList().FillCaloHitList(caloHitList2);

    for (const CartesianVector &positionVector1 : pCluster1)
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

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArMuonLeadingHelper::GetClosestPosition(const CartesianVector &referencePoint, const CartesianPointVector &cartesianPointVector,
    const Cluster *const pCluster)
{
    CartesianVector closestPoint(0.f, 0.f, 0.f);
    float shortestDistanceSquared(std::numeric_limits<float>::max());

    for (const CartesianVector &testPosition : cartesianPointVector)
    {
        if (LArClusterHelper::GetClosestDistance(testPosition, pCluster) > 0.5f)
            continue;

        const float separationSquared((testPosition - referencePoint).GetMagnitude());

        if (separationSquared > 5.f)
            continue;

        if (separationSquared < shortestDistanceSquared)
        {
            shortestDistanceSquared = separationSquared;
            closestPoint = testPosition;
        }
    }

    return closestPoint;
}

} // namespace lar_content

