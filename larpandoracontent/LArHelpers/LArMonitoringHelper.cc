/**
 *  @file   larpandoracontent/LArHelpers/LArMonitoringHelper.cc
 *
 *  @brief  Implementation of the lar monitoring helper class.
 *
 *  $Log: $
 */

#include "Objects/CaloHit.h"
#include "Objects/Cluster.h"
#include "Objects/MCParticle.h"
#include "Objects/ParticleFlowObject.h"

#include "Helpers/MCParticleHelper.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

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
    const CaloHitToMCMap &hitToPrimaryMCMap, MCToPfoMap &outputPrimaryMap)
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

        MCContributionMap inputContributionMap;

        for (const CaloHit *const pCaloHit : clusterHits)
        {
            if (!caloHitSet.count(pCaloHit))
                continue;

            CaloHitToMCMap::const_iterator mcIter = hitToPrimaryMCMap.find(pCaloHit);

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

void LArMonitoringHelper::GetMCParticleToCaloHitMatches(const CaloHitList *const pCaloHitList, const MCRelationMap &mcToPrimaryMCMap,
    CaloHitToMCMap &hitToPrimaryMCMap, MCContributionMap &mcToTrueHitListMap)
{
    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        try
        {
            const MCParticle *const pHitParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));

            MCRelationMap::const_iterator mcIter = mcToPrimaryMCMap.find(pHitParticle);

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
    CaloHitToPfoMap &hitToPfoMap, PfoContributionMap &pfoToHitListMap)
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

void LArMonitoringHelper::GetMCParticleToPfoMatches(const CaloHitList *const pCaloHitList, const PfoContributionMap &pfoToHitListMap,
    const CaloHitToMCMap &hitToPrimaryMCMap, MCToPfoMap &mcToBestPfoMap, MCContributionMap &mcToBestPfoHitsMap,
    MCToPfoMatchingMap &mcToFullPfoMatchingMap)
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

            CaloHitToMCMap::const_iterator mcIter = hitToPrimaryMCMap.find(pCaloHit);

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
        const PfoContributionMap &pfoContributionMap(mcToFullPfoMatchingMap.at(pPrimaryParticle));
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

} // namespace lar_content
