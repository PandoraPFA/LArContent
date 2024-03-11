/**
 *  @file   larpandoracontent/LArHelpers/LArMonitoringHelper.cc
 *
 *  @brief  Implementation of the lar monitoring helper class.
 *
 *  $Log: $
 */

#include "Pandora/PdgTable.h"

#include "Objects/CaloHit.h"
#include "Objects/MCParticle.h"
#include "Objects/ParticleFlowObject.h"

#include "Helpers/MCParticleHelper.h"

#include "larpandoracontent/LArHelpers/LArFormattingHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

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

void LArMonitoringHelper::GetOrderedMCParticleVector(
    const LArMCParticleHelper::MCContributionMapVector &selectedMCParticleToGoodHitsMaps, MCParticleVector &orderedMCParticleVector)
{
    for (const LArMCParticleHelper::MCContributionMap &mcParticleToGoodHitsMap : selectedMCParticleToGoodHitsMaps)
    {
        if (mcParticleToGoodHitsMap.empty())
            continue;

        // Copy map contents to vector it can be sorted
        std::vector<LArMCParticleHelper::MCParticleCaloHitListPair> mcParticleToGoodHitsVect;
        std::copy(mcParticleToGoodHitsMap.begin(), mcParticleToGoodHitsMap.end(), std::back_inserter(mcParticleToGoodHitsVect));

        // Sort by number of hits descending
        std::sort(mcParticleToGoodHitsVect.begin(), mcParticleToGoodHitsVect.end(),
            [](const LArMCParticleHelper::MCParticleCaloHitListPair &a, const LArMCParticleHelper::MCParticleCaloHitListPair &b) -> bool
            {
                // Neutrinos, then beam particles, then cosmic rays
                const bool isANuFinalState(LArMCParticleHelper::IsBeamNeutrinoFinalState(a.first)),
                    isBNuFinalState(LArMCParticleHelper::IsBeamNeutrinoFinalState(b.first));

                if (isANuFinalState != isBNuFinalState)
                    return isANuFinalState;

                const bool isABeamParticle(LArMCParticleHelper::IsBeamParticle(a.first)),
                    isBBeamParticle(LArMCParticleHelper::IsBeamParticle(b.first));

                if (isABeamParticle != isBBeamParticle)
                    return isABeamParticle;

                // Then sort by numbers of true hits
                if (a.second.size() != b.second.size())
                    return (a.second.size() > b.second.size());

                // Default to normal MCParticle sorting
                return LArMCParticleHelper::SortByMomentum(a.first, b.first);
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
    std::sort(pfoToReconstructable2DHitsVect.begin(), pfoToReconstructable2DHitsVect.end(),
        [](const LArMCParticleHelper::PfoCaloHitListPair &a, const LArMCParticleHelper::PfoCaloHitListPair &b) -> bool
        {
            // Neutrinos before cosmic rays
            const bool isANuFinalState(LArPfoHelper::IsNeutrinoFinalState(a.first)), isBNuFinalState(LArPfoHelper::IsNeutrinoFinalState(b.first));

            if (isANuFinalState != isBNuFinalState)
                return isANuFinalState;

            if (a.second.size() != b.second.size())
                return (a.second.size() > b.second.size());

            // Default to normal pfo sorting
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

void LArMonitoringHelper::PrintMCParticleTable(
    const LArMCParticleHelper::MCContributionMap &selectedMCParticleToGoodHitsMap, const MCParticleVector &orderedMCParticleVector)
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
            continue; // ATTN MCParticles in selectedMCParticleToGoodHitsMap may be a subset of orderedMCParticleVector

        // ATTN enumerate from 1 to match event validation algorithm
        table.AddElement(id + 1);
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

    LArFormattingHelper::Table table({"ID", "PID", "Is Nu FS", "", "nGoodHits", "U", "V", "W"});

    for (unsigned int id = 0; id < orderedPfoVector.size(); ++id)
    {
        const ParticleFlowObject *const pPfo(orderedPfoVector.at(id));

        LArMCParticleHelper::PfoContributionMap::const_iterator it = pfoToReconstructable2DHitsMap.find(pPfo);
        if (pfoToReconstructable2DHitsMap.end() == it)
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);

        // ATTN enumerate from 1 to match event validation algorithm
        table.AddElement(id + 1);
        table.AddElement(pPfo->GetParticleId());
        table.AddElement(LArPfoHelper::IsNeutrinoFinalState(pPfo));

        table.AddElement(it->second.size());
        table.AddElement(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, it->second));
        table.AddElement(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, it->second));
        table.AddElement(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, it->second));
    }

    table.Print();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMonitoringHelper::PrintMatchingTable(const PfoVector &orderedPfoVector, const MCParticleVector &orderedMCParticleVector,
    const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcParticleToPfoHitSharingMap, const unsigned int nMatches)
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

    // Get the maximum number of MCParticle to Pfos matches that need to be shown
    unsigned int maxMatches(0);
    for (const auto &entry : mcParticleToPfoHitSharingMap)
        maxMatches = std::max(static_cast<unsigned int>(entry.second.size()), maxMatches);

    const bool showOthersColumn(maxMatches > nMatches);
    const unsigned int nMatchesToShow(std::min(maxMatches, nMatches));

    // Set up the table headers
    std::vector<std::string> tableHeaders({"MCParticle", ""});
    for (unsigned int i = 0; i < nMatchesToShow; ++i)
    {
        tableHeaders.push_back("");
        tableHeaders.push_back("Pfo");
        tableHeaders.push_back("nSharedHits");
    }

    if (showOthersColumn)
    {
        tableHeaders.push_back("");
        tableHeaders.push_back("");
        tableHeaders.push_back("nOtherPfos");
        tableHeaders.push_back("nSharedHits");
    }

    LArFormattingHelper::Table table(tableHeaders);

    // Make a new row for each MCParticle
    for (unsigned int mcParticleId = 0; mcParticleId < orderedMCParticleVector.size(); ++mcParticleId)
    {
        const MCParticle *const pMCParticle(orderedMCParticleVector.at(mcParticleId));
        LArMCParticleHelper::MCParticleToPfoHitSharingMap::const_iterator it = mcParticleToPfoHitSharingMap.find(pMCParticle);
        LArMCParticleHelper::PfoToSharedHitsVector pfoToSharedHitsVector;

        if (it != mcParticleToPfoHitSharingMap.end())
            pfoToSharedHitsVector = it->second;

        const LArFormattingHelper::Color mcCol(LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCParticle)
                ? LArFormattingHelper::LIGHT_GREEN
                : (LArMCParticleHelper::IsBeamParticle(pMCParticle) ? LArFormattingHelper::LIGHT_BLUE : LArFormattingHelper::LIGHT_RED));

        // ATTN enumerate from 1 to match event validation algorithm
        table.AddElement(mcParticleId + 1, LArFormattingHelper::REGULAR, mcCol);

        // Get the matched Pfos
        unsigned int nPfosShown(0);
        unsigned int nOtherHits(0);
        for (const auto &pfoNSharedHitsPair : pfoToSharedHitsVector)
        {
            for (unsigned int pfoId = 0; pfoId < orderedPfoVector.size(); ++pfoId)
            {
                if (pfoNSharedHitsPair.first != orderedPfoVector.at(pfoId))
                    continue;

                if (nPfosShown < nMatchesToShow)
                {
                    // ATTN enumerate from 1 to match event validation algorithm
                    const LArFormattingHelper::Color pfoCol(
                        LArPfoHelper::IsNeutrinoFinalState(pfoNSharedHitsPair.first) ? LArFormattingHelper::LIGHT_GREEN : LArFormattingHelper::LIGHT_RED);
                    table.AddElement(pfoId + 1, LArFormattingHelper::REGULAR, pfoCol);
                    table.AddElement(pfoNSharedHitsPair.second.size(), LArFormattingHelper::REGULAR, pfoCol);
                    nPfosShown++;
                }
                else
                {
                    nOtherHits += pfoNSharedHitsPair.second.size();
                }
                break;
            }
        }

        // Pad the rest of the row with empty entries
        for (unsigned int i = 0; i < nMatchesToShow - nPfosShown; ++i)
        {
            table.AddElement("");
            table.AddElement("");
        }

        // Print any remaining matches
        if (!showOthersColumn)
            continue;

        if (nOtherHits != 0)
        {
            table.AddElement(pfoToSharedHitsVector.size() - nPfosShown);
            table.AddElement(nOtherHits);
        }
        else
        {
            table.AddElement("");
            table.AddElement("");
        }
    }
    table.Print();
}

} // namespace lar_content
