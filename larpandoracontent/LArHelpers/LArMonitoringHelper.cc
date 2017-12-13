/**
 *  @file   larpandoracontent/LArHelpers/LArMonitoringHelper.cc
 *
 *  @brief  Implementation of the lar monitoring helper class.
 *
 *  $Log: $
 */

#include "Pandora/PdgTable.h"

#include "Objects/CaloHit.h"
#include "Objects/Cluster.h"
#include "Objects/MCParticle.h"
#include "Objects/ParticleFlowObject.h"

#include "Helpers/MCParticleHelper.h"

#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArFormattingHelper.h"

#include <algorithm>

using namespace pandora;

namespace lar_content
{

void LArMonitoringHelper::ExtractTargetPfos(const PfoList &inputPfoList, const bool primaryPfosOnly, PfoList &outputPfoList)
{
    if (primaryPfosOnly)
    {
        for (const ParticleFlowObject *const pPfo : inputPfoList)
        {
            if (LArPfoHelper::IsFinalState(pPfo))
            {
                outputPfoList.push_back(pPfo);
            }
            else if (pPfo->GetParentPfoList().empty() && LArPfoHelper::IsNeutrino(pPfo))
            {
                outputPfoList.insert(outputPfoList.end(), pPfo->GetDaughterPfoList().begin(), pPfo->GetDaughterPfoList().end());
            }
        }
    }
    else
    {
        LArPfoHelper::GetAllDownstreamPfos(inputPfoList, outputPfoList);

        for (const ParticleFlowObject *const pPfo : inputPfoList)
        {
            if (LArPfoHelper::IsNeutrino(pPfo))
            {
                PfoList::iterator eraseIter(std::find(outputPfoList.begin(), outputPfoList.end(), pPfo));

                if (outputPfoList.end() != eraseIter)
                    outputPfoList.erase(eraseIter);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMonitoringHelper::GetNeutrinoMatches(const CaloHitList *const pCaloHitList, const PfoList &recoNeutrinos,
    const LArMCParticleHelper::CaloHitToMCMap &hitToPrimaryMCMap, LArMCParticleHelper::MCToPfoMap &outputPrimaryMap)
{
    const CaloHitSet caloHitSet(pCaloHitList->begin(), pCaloHitList->end());

    for (const ParticleFlowObject *const pNeutrinoPfo : recoNeutrinos)
    {
        if (!LArPfoHelper::IsNeutrino(pNeutrinoPfo))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        PfoList pfoList;
        LArPfoHelper::GetAllDownstreamPfos(pNeutrinoPfo, pfoList);

        CaloHitList clusterHits;
        LArMonitoringHelper::CollectCaloHits(pfoList, clusterHits);

        LArMCParticleHelper::MCContributionMap inputContributionMap;

        for (const CaloHit *const pCaloHit : clusterHits)
        {
            if (!caloHitSet.count(pCaloHit))
                continue;

            LArMCParticleHelper::CaloHitToMCMap::const_iterator mcIter = hitToPrimaryMCMap.find(pCaloHit);

            if (mcIter == hitToPrimaryMCMap.end())
                continue;

            const MCParticle *const pFinalStateParticle = mcIter->second;
            const MCParticle *const pNeutrinoParticle = LArMCParticleHelper::GetParentMCParticle(pFinalStateParticle);

            if (!LArMCParticleHelper::IsNeutrino(pNeutrinoParticle))
                continue;

            inputContributionMap[pNeutrinoParticle].push_back(pCaloHit);
        }

        CaloHitList biggestHitList;
        const MCParticle *biggestContributor(nullptr);

        MCParticleList mcParticleList;
        for (const auto &mapEntry : inputContributionMap) mcParticleList.push_back(mapEntry.first);
        mcParticleList.sort(LArMCParticleHelper::SortByMomentum);

        for (const MCParticle *const pMCParticle : mcParticleList)
        {
            const CaloHitList &caloHitList(inputContributionMap.at(pMCParticle));

            if (caloHitList.size() > biggestHitList.size())
            {
                biggestHitList = caloHitList;
                biggestContributor = pMCParticle;
            }
        }

        if (biggestContributor)
            outputPrimaryMap[biggestContributor] = pNeutrinoPfo;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMonitoringHelper::GetMCParticleToCaloHitMatches(const CaloHitList *const pCaloHitList, const LArMCParticleHelper::MCRelationMap &mcToPrimaryMCMap,
    LArMCParticleHelper::CaloHitToMCMap &hitToPrimaryMCMap, LArMCParticleHelper::MCContributionMap &mcToTrueHitListMap)
{
    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        try
        {
            const MCParticle *const pHitParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));

            LArMCParticleHelper::MCRelationMap::const_iterator mcIter = mcToPrimaryMCMap.find(pHitParticle);

            if (mcToPrimaryMCMap.end() == mcIter)
                continue;

            const MCParticle *const pPrimaryParticle = mcIter->second;
            mcToTrueHitListMap[pPrimaryParticle].push_back(pCaloHit);
            hitToPrimaryMCMap[pCaloHit] = pPrimaryParticle;
        }
        catch (StatusCodeException &statusCodeException)
        {
            if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                throw statusCodeException;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMonitoringHelper::GetPfoToCaloHitMatches(const CaloHitList *const pCaloHitList, const PfoList &pfoList, const bool collapseToPrimaryPfos,
    LArMCParticleHelper::CaloHitToPfoMap &hitToPfoMap, LArMCParticleHelper::PfoContributionMap &pfoToHitListMap)
{
    const CaloHitSet caloHitSet(pCaloHitList->begin(), pCaloHitList->end());

    for (const ParticleFlowObject *const pPfo : pfoList)
    {
        ClusterList clusterList;

        if (!collapseToPrimaryPfos)
        {
            LArPfoHelper::GetTwoDClusterList(pPfo, clusterList);
        }
        else
        {
            if (!LArPfoHelper::IsFinalState(pPfo))
                continue;

            PfoList downstreamPfoList;
            LArPfoHelper::GetAllDownstreamPfos(pPfo, downstreamPfoList);

            for (PfoList::const_iterator dIter = downstreamPfoList.begin(), dIterEnd = downstreamPfoList.end(); dIter != dIterEnd; ++dIter)
                LArPfoHelper::GetTwoDClusterList(*dIter, clusterList);
        }

        CaloHitList pfoHitList;

        for (const Cluster *const pCluster : clusterList)
        {
            CaloHitList clusterHits;
            pCluster->GetOrderedCaloHitList().FillCaloHitList(clusterHits);
            clusterHits.insert(clusterHits.end(), pCluster->GetIsolatedCaloHitList().begin(), pCluster->GetIsolatedCaloHitList().end());

            for (const CaloHit *const pCaloHit : clusterHits)
            {
                if (TPC_3D == pCaloHit->GetHitType())
                    throw StatusCodeException(STATUS_CODE_FAILURE);

                if (!caloHitSet.count(pCaloHit))
                    continue;

                hitToPfoMap[pCaloHit] = pPfo;
                pfoHitList.push_back(pCaloHit);
            }
        }

        CaloHitList &caloHitListInMap(pfoToHitListMap[pPfo]);
        caloHitListInMap.insert(caloHitListInMap.end(), pfoHitList.begin(), pfoHitList.end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMonitoringHelper::GetMCParticleToPfoMatches(const CaloHitList *const pCaloHitList, const LArMCParticleHelper::PfoContributionMap &pfoToHitListMap,
    const LArMCParticleHelper::CaloHitToMCMap &hitToPrimaryMCMap, LArMCParticleHelper::MCToPfoMap &mcToBestPfoMap, LArMCParticleHelper::MCContributionMap &mcToBestPfoHitsMap,
    LArMCParticleHelper::MCToPfoMatchingMap &mcToFullPfoMatchingMap)
{
    const CaloHitSet caloHitSet(pCaloHitList->begin(), pCaloHitList->end());

    PfoList pfoList;
    for (const auto &mapEntry : pfoToHitListMap) pfoList.push_back(mapEntry.first);
    pfoList.sort(LArPfoHelper::SortByNHits);

    for (const ParticleFlowObject *const pPfo : pfoList)
    {
        const CaloHitList &pfoHits(pfoToHitListMap.at(pPfo));

        for (const CaloHit *const pCaloHit : pfoHits)
        {
            if (TPC_3D == pCaloHit->GetHitType())
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

            if (!caloHitSet.count(pCaloHit))
                continue;

            LArMCParticleHelper::CaloHitToMCMap::const_iterator mcIter = hitToPrimaryMCMap.find(pCaloHit);

            if (mcIter == hitToPrimaryMCMap.end())
                continue;

            const MCParticle *const pPrimaryParticle = mcIter->second;
            mcToFullPfoMatchingMap[pPrimaryParticle][pPfo].push_back(pCaloHit);
        }
    }

    MCParticleList mcParticleList;
    for (const auto &mapEntry : mcToFullPfoMatchingMap) mcParticleList.push_back(mapEntry.first);
    mcParticleList.sort(LArMCParticleHelper::SortByMomentum);

    for (const MCParticle *const pPrimaryParticle : mcParticleList)
    {
        const LArMCParticleHelper::PfoContributionMap &pfoContributionMap(mcToFullPfoMatchingMap.at(pPrimaryParticle));
        const ParticleFlowObject *pBestMatchPfo(nullptr);
        CaloHitList bestMatchCaloHitList;

        PfoList matchedPfoList;
        for (const auto &mapEntry : pfoContributionMap) matchedPfoList.push_back(mapEntry.first);
        matchedPfoList.sort(LArPfoHelper::SortByNHits);

        for (const ParticleFlowObject *const pPfo : matchedPfoList)
        {
            const CaloHitList &caloHitList(pfoContributionMap.at(pPfo));

            if (caloHitList.size() > bestMatchCaloHitList.size())
            {
                pBestMatchPfo = pPfo;
                bestMatchCaloHitList = caloHitList;
            }
        }

        if (pBestMatchPfo)
        {
            mcToBestPfoMap[pPrimaryParticle] = pBestMatchPfo;
            mcToBestPfoHitsMap[pPrimaryParticle] = bestMatchCaloHitList;
        }
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void LArMonitoringHelper::CollectCaloHits(const ParticleFlowObject *const pParentPfo, CaloHitList &caloHitList)
{
    ClusterList clusterList;
    LArPfoHelper::GetTwoDClusterList(pParentPfo, clusterList);

    for (const Cluster *const pCluster : clusterList)
    {
        if (TPC_3D == LArClusterHelper::GetClusterHitType(pCluster))
            throw StatusCodeException(STATUS_CODE_FAILURE);

        pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);
        caloHitList.insert(caloHitList.end(), pCluster->GetIsolatedCaloHitList().begin(), pCluster->GetIsolatedCaloHitList().end());
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void LArMonitoringHelper::CollectCaloHits(const PfoList &pfoList, CaloHitList &caloHitList)
{
    for (const ParticleFlowObject *const pPfo : pfoList)
    {
        LArMonitoringHelper::CollectCaloHits(pPfo, caloHitList);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int LArMonitoringHelper::CountHitsByType(const HitType hitType, const CaloHitList &caloHitList)
{
    unsigned int nHitsOfSpecifiedType(0);

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        if (hitType == pCaloHit->GetHitType())
            ++nHitsOfSpecifiedType;
    }

    return nHitsOfSpecifiedType;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMonitoringHelper::GetOrderedMCParticleVector(const LArMCParticleHelper::LArMCParticleHelper::MCContributionMapVector &selectedMCParticleToGoodHitsMaps, MCParticleVector &orderedMCParticleVector)
{
    for (const LArMCParticleHelper::MCContributionMap &mcParticleToGoodHitsMap : selectedMCParticleToGoodHitsMaps)
    {
        if (mcParticleToGoodHitsMap.empty())
            continue;

        // Copy map contents to vector it can be sorted
        std::vector<LArMCParticleHelper::MCParticleCaloHitListPair> mcParticleToGoodHitsVect;
        std::copy(mcParticleToGoodHitsMap.begin(), mcParticleToGoodHitsMap.end(), std::back_inserter(mcParticleToGoodHitsVect));

        // Sort by number of hits descending
        std::sort(mcParticleToGoodHitsVect.begin(), mcParticleToGoodHitsVect.end(), [] (const LArMCParticleHelper::MCParticleCaloHitListPair &a, const LArMCParticleHelper::MCParticleCaloHitListPair &b) -> bool
        {
            if (a.second.size() != b.second.size())
                return (a.second.size() > b.second.size());

            // ATTN default to normal MCParticle sorting to avoid tie-breakers
            return (a.first < b.first);
        });

        for (const LArMCParticleHelper::MCParticleCaloHitListPair &mcParticleCaloHitPair : mcParticleToGoodHitsVect)
            orderedMCParticleVector.push_back(mcParticleCaloHitPair.first);
    }

    // Check that all elements of the vector are unique
    const unsigned int nMCParticles(orderedMCParticleVector.size());
    if (std::distance(orderedMCParticleVector.begin(), std::unique(orderedMCParticleVector.begin(), orderedMCParticleVector.end())) != nMCParticles)
        throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMonitoringHelper::GetOrderedPfoVector(const LArMCParticleHelper::PfoContributionMap &pfoToReconstructable2DHitsMap, pandora::PfoVector &orderedPfoVector)
{
    // Copy map contents to vector it can be sorted
    std::vector<LArMCParticleHelper::PfoCaloHitListPair> pfoToReconstructable2DHitsVect;
    std::copy(pfoToReconstructable2DHitsMap.begin(), pfoToReconstructable2DHitsMap.end(), std::back_inserter(pfoToReconstructable2DHitsVect));

    // Sort by number of hits descending putting neutrino final states first
    std::sort(pfoToReconstructable2DHitsVect.begin(), pfoToReconstructable2DHitsVect.end(), [] (const LArMCParticleHelper::PfoCaloHitListPair &a, const LArMCParticleHelper::PfoCaloHitListPair &b) -> bool
    {
        bool isANuFinalState(LArPfoHelper::IsNeutrinoFinalState(a.first));
        bool isBNuFinalState(LArPfoHelper::IsNeutrinoFinalState(b.first));

        if (isANuFinalState != isBNuFinalState)
            return isANuFinalState;

        if (a.second.size() != b.second.size())
            return (a.second.size() > b.second.size());

        // ATTN fall back on using all hits as a tie-breaker
        return LArPfoHelper::SortByNHits(a.first, b.first);
    });

    for (const LArMCParticleHelper::PfoCaloHitListPair &pfoCaloHitPair : pfoToReconstructable2DHitsVect)
        orderedPfoVector.push_back(pfoCaloHitPair.first);

    // Check that all elements of the vector are unique
    const unsigned int nPfos(orderedPfoVector.size());
    if (std::distance(orderedPfoVector.begin(), std::unique(orderedPfoVector.begin(), orderedPfoVector.end())) != nPfos)
        throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMonitoringHelper::PrintMCParticleTable(const LArMCParticleHelper::MCContributionMap &selectedMCParticleToGoodHitsMap, const MCParticleVector &orderedMCParticleVector)
{
    if (selectedMCParticleToGoodHitsMap.empty())
    {
        std::cout << "No MCParticles supplied." << std::endl;
        return;
    }

    LArFormattingHelper::Table table({"ID", "NUANCE", "TYPE", "", "E", "dist", "", "nGoodHits", "U", "V", "W"});

    unsigned int usedParticleCount(0);
    for (unsigned int id = 0; id < orderedMCParticleVector.size(); ++id)
    {
        const MCParticle *const pMCParticle(orderedMCParticleVector.at(id));

        LArMCParticleHelper::MCContributionMap::const_iterator it = selectedMCParticleToGoodHitsMap.find(pMCParticle);
        if (selectedMCParticleToGoodHitsMap.end() == it)
            continue;  // ATTN MCParticles in selectedMCParticleToGoodHitsMap may be a subset of orderedMCParticleVector

        table.AddElement(id);
        table.AddElement(LArMCParticleHelper::GetNuanceCode(pMCParticle));
        table.AddElement(PdgTable::GetParticleName(pMCParticle->GetParticleId()));

        table.AddElement(pMCParticle->GetEnergy());
        table.AddElement((pMCParticle->GetEndpoint() - pMCParticle->GetVertex()).GetMagnitude());

        table.AddElement(it->second.size());
        table.AddElement(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, it->second));
        table.AddElement(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, it->second));
        table.AddElement(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, it->second));

        usedParticleCount++;
    }

    // Check every MCParticle in selectedMCParticleToGoodHitsMap has been printed
    if (usedParticleCount != selectedMCParticleToGoodHitsMap.size())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    table.Print();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMonitoringHelper::PrintPfoTable(const LArMCParticleHelper::PfoContributionMap &pfoToReconstructable2DHitsMap, const PfoVector &orderedPfoVector)
{
    if (pfoToReconstructable2DHitsMap.empty())
    {
        std::cout << "No Pfos supplied." << std::endl;
        return;
    }

    LArFormattingHelper::Table table({"ID", "PID", "Is Nu FS", "", "nHits", "U", "V", "W", "", "nGoodHits", "U", "V", "W"});

    for (unsigned int id = 0; id < orderedPfoVector.size(); ++id)
    {
        const ParticleFlowObject *const pPfo(orderedPfoVector.at(id));

        LArMCParticleHelper::PfoContributionMap::const_iterator it = pfoToReconstructable2DHitsMap.find(pPfo);
        if (pfoToReconstructable2DHitsMap.end() == it)
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);

        table.AddElement(id);
        table.AddElement(pPfo->GetParticleId());
        table.AddElement(LArPfoHelper::IsNeutrinoFinalState(pPfo));

        CaloHitList all2DCaloHits;
        LArMonitoringHelper::CollectCaloHits(pPfo, all2DCaloHits);

        table.AddElement(all2DCaloHits.size());
        table.AddElement(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, all2DCaloHits));
        table.AddElement(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, all2DCaloHits));
        table.AddElement(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, all2DCaloHits));

        table.AddElement(it->second.size());
        table.AddElement(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, it->second));
        table.AddElement(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, it->second));
        table.AddElement(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, it->second));
    }

    table.Print();
}

//------------------------------------------------------------------------------------------------------------------------------------------

/*
void LArMonitoringHelper::PrintMatchingTable(const PfoVector &orderedPfoVector, const MCParticleVector &orderedMCParticleVector, const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcParticleToPfoHitSharingMap)
{
    if (orderedPfoVector.empty())
    {
        std::cout << "No Pfos supplied." << std::endl;
        return;
    }

    if (orderedMCParticleVector.empty())
    {
        std::cout << "No MCParticles supplied." << std::endl;
        return;
    }

    std::vector<std::string> tableHeaders({"MCParticle",""});
    for (unsigned int pfoId = 0; pfoId < orderedPfoVector.size(); ++pfoId)
        tableHeaders.push_back(std::to_string(pfoId));

    LArFormattingHelper::Table table(tableHeaders);

    for (unsigned int mcParticleId = 0; mcParticleId < orderedMCParticleVector.size(); ++mcParticleId)
    {
        LArMCParticleHelper::MCParticleToPfoHitSharingMap::const_iterator it = mcParticleToPfoHitSharingMap.find(orderedMCParticleVector.at(mcParticleId));
        if (it == mcParticleToPfoHitSharingMap.end())
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);

        table.AddElement(mcParticleId);
        for (unsigned int pfoId = 0; pfoId < orderedPfoVector.size(); ++pfoId)
        {
            bool foundPfo(false);
            unsigned int nSharedHits(std::numeric_limits<unsigned int>::max());

            for (const LArMCParticleHelper::PfoIntPair &pair : it->second)
            {
                if (pair.first == orderedPfoVector.at(pfoId))
                {
                    nSharedHits = pair.second;
                    foundPfo = true;
                    break;
                }
            }

            if (!foundPfo)
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);

            table.AddElement(nSharedHits);
        }
    }

    table.Print();
}
*/

} // namespace lar_content
