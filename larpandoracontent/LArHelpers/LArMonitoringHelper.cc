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

void LArMonitoringHelper::GetOrderedMCParticleVector(const LArMCParticleHelper::MCContributionMapVector &selectedMCParticleToGoodHitsMaps, MCParticleVector &orderedMCParticleVector)
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
            return (*(a.first) < *(b.first));
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
        LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_U, all2DCaloHits);
        LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_V, all2DCaloHits);
        LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, all2DCaloHits);

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
