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
        for (PfoList::const_iterator iter = inputPfoList.begin(), iterEnd = inputPfoList.end(); iter != iterEnd; ++iter)
        {
            const ParticleFlowObject *const pPfo(*iter);

            if (LArPfoHelper::IsFinalState(pPfo))
            {
                outputPfoList.insert(pPfo);
            }
            else if (pPfo->GetParentPfoList().empty() && LArPfoHelper::IsNeutrino(pPfo))
            {
                outputPfoList.insert(pPfo->GetDaughterPfoList().begin(), pPfo->GetDaughterPfoList().end());
            }
        }
    }
    else
    {
        LArPfoHelper::GetAllDownstreamPfos(inputPfoList, outputPfoList);

        for (PfoList::const_iterator iter = inputPfoList.begin(), iterEnd = inputPfoList.end(); iter != iterEnd; ++iter)
        {
            const ParticleFlowObject *const pPfo(*iter);

            if (LArPfoHelper::IsNeutrino(pPfo))
                outputPfoList.erase(pPfo);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMonitoringHelper::GetNeutrinoMatches(const CaloHitList *const pCaloHitList, const PfoList &recoNeutrinos,
    const CaloHitToMCMap &hitToPrimaryMCMap, MCToPfoMap &outputPrimaryMap)
{
    for (PfoList::const_iterator pIter = recoNeutrinos.begin(), pIterEnd = recoNeutrinos.end(); pIter != pIterEnd; ++pIter)
    {
        const ParticleFlowObject *const pNeutrinoPfo = *pIter;

        if (!LArPfoHelper::IsNeutrino(pNeutrinoPfo))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        PfoList pfoList;
        LArPfoHelper::GetAllDownstreamPfos(pNeutrinoPfo, pfoList);

        CaloHitList clusterHits;
        LArMonitoringHelper::CollectCaloHits(pfoList, clusterHits);

        MCContributionMap inputContributionMap;

        for (CaloHitList::const_iterator hIter = clusterHits.begin(), hIterEnd = clusterHits.end(); hIter != hIterEnd; ++hIter)
        {
            const CaloHit *const pCaloHit = *hIter;

            if (pCaloHitList->count(pCaloHit) == 0)
                continue;

            CaloHitToMCMap::const_iterator mcIter = hitToPrimaryMCMap.find(pCaloHit);

            if (mcIter == hitToPrimaryMCMap.end())
                continue;

            const MCParticle *const pFinalStateParticle = mcIter->second;
            const MCParticle *const pNeutrinoParticle = LArMCParticleHelper::GetParentMCParticle(pFinalStateParticle);

            if (!LArMCParticleHelper::IsNeutrino(pNeutrinoParticle))
                continue;

            inputContributionMap[pNeutrinoParticle].insert(pCaloHit);
        }

        CaloHitList biggestHitList;
        const MCParticle *biggestContributor = NULL;

        for (MCContributionMap::const_iterator cIter = inputContributionMap.begin(), cIterEnd = inputContributionMap.end(); cIter != cIterEnd; ++cIter)
        {
            if (cIter->second.size() > biggestHitList.size())
            {
                biggestHitList = cIter->second;
                biggestContributor = cIter->first;
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
    for (CaloHitList::const_iterator iter = pCaloHitList->begin(), iterEnd = pCaloHitList->end(); iter != iterEnd; ++iter)
    {
        try
        {
            const CaloHit *const pCaloHit = *iter;
            const MCParticle *const pHitParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));

            MCRelationMap::const_iterator mcIter = mcToPrimaryMCMap.find(pHitParticle);

            if (mcToPrimaryMCMap.end() == mcIter)
                continue;

            const MCParticle *const pPrimaryParticle = mcIter->second;
            mcToTrueHitListMap[pPrimaryParticle].insert(pCaloHit);
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
    for (PfoList::const_iterator pIter = pfoList.begin(), pIterEnd = pfoList.end(); pIter != pIterEnd; ++pIter)
    {
        const ParticleFlowObject *const pPfo = *pIter;
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

        for (ClusterList::const_iterator cIter = clusterList.begin(), cIterEnd = clusterList.end(); cIter != cIterEnd; ++cIter)
        {
            const Cluster *const pCluster = *cIter;

            CaloHitList clusterHits;
            pCluster->GetOrderedCaloHitList().GetCaloHitList(clusterHits);
            clusterHits.insert(pCluster->GetIsolatedCaloHitList().begin(), pCluster->GetIsolatedCaloHitList().end());

            for (CaloHitList::const_iterator hIter = clusterHits.begin(), hIterEnd = clusterHits.end(); hIter != hIterEnd; ++hIter)
            {
                const CaloHit *const pCaloHit = *hIter;

                if (TPC_3D == pCaloHit->GetHitType())
                    throw StatusCodeException(STATUS_CODE_FAILURE);

                if (pCaloHitList->count(pCaloHit) == 0)
                    continue;

                hitToPfoMap[pCaloHit] = pPfo;
                pfoHitList.insert(pCaloHit);
            }
        }

        pfoToHitListMap[pPfo] = pfoHitList;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMonitoringHelper::GetMCParticleToPfoMatches(const CaloHitList *const pCaloHitList, const PfoContributionMap &pfoToHitListMap,
    const CaloHitToMCMap &hitToPrimaryMCMap, MCToPfoMap &mcToBestPfoMap, MCContributionMap &mcToBestPfoHitsMap,
    MCToPfoMatchingMap &mcToFullPfoMatchingMap)
{
    for (PfoContributionMap::const_iterator iter = pfoToHitListMap.begin(), iterEnd = pfoToHitListMap.end(); iter != iterEnd; ++iter)
    {
        const ParticleFlowObject *const pPfo(iter->first);
        const CaloHitList &pfoHits(iter->second);

        for (CaloHitList::const_iterator hIter = pfoHits.begin(), hIterEnd = pfoHits.end(); hIter != hIterEnd; ++hIter)
        {
            const CaloHit *const pCaloHit = *hIter;

            if (TPC_3D == pCaloHit->GetHitType())
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

            if (pCaloHitList->count(pCaloHit) == 0)
                continue;

            CaloHitToMCMap::const_iterator mcIter = hitToPrimaryMCMap.find(pCaloHit);

            if (mcIter == hitToPrimaryMCMap.end())
                continue;

            const MCParticle *const pPrimaryParticle = mcIter->second;
            mcToFullPfoMatchingMap[pPrimaryParticle][pPfo].insert(pCaloHit);
        }
    }

    for (MCToPfoMatchingMap::const_iterator iter = mcToFullPfoMatchingMap.begin(), iterEnd = mcToFullPfoMatchingMap.end(); iter != iterEnd; ++iter)
    {
        const MCParticle *const pPrimaryParticle = iter->first;
        const PfoContributionMap &pfoContributionMap = iter->second;

        const ParticleFlowObject *pBestMatchPfo(NULL);
        CaloHitList bestMatchCaloHitList;

        for (PfoContributionMap::const_iterator pIter = pfoContributionMap.begin(), pIterEnd = pfoContributionMap.end(); pIter != pIterEnd; ++pIter)
        {
            if (pIter->second.size() > bestMatchCaloHitList.size())
            {
                pBestMatchPfo = pIter->first;
                bestMatchCaloHitList = pIter->second;
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

    for (ClusterList::const_iterator cIter = clusterList.begin(), cIterEnd = clusterList.end(); cIter != cIterEnd; ++cIter)
    {
        const Cluster *const pCluster = *cIter;

        if (TPC_3D == LArClusterHelper::GetClusterHitType(pCluster))
            throw StatusCodeException(STATUS_CODE_FAILURE);

        pCluster->GetOrderedCaloHitList().GetCaloHitList(caloHitList);
        caloHitList.insert(pCluster->GetIsolatedCaloHitList().begin(), pCluster->GetIsolatedCaloHitList().end());
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void LArMonitoringHelper::CollectCaloHits(const PfoList &pfoList, CaloHitList &caloHitList)
{
    for (PfoList::const_iterator pIter = pfoList.begin(), pIterEnd = pfoList.end(); pIter != pIterEnd; ++pIter)
    {
        const ParticleFlowObject *const pPfo = *pIter;
        LArMonitoringHelper::CollectCaloHits(pPfo, caloHitList);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int LArMonitoringHelper::CountHitsByType(const HitType hitType, const CaloHitList &caloHitList)
{
    unsigned int nHitsOfSpecifiedType(0);

    for (CaloHitList::const_iterator iter = caloHitList.begin(), iterEnd = caloHitList.end(); iter != iterEnd; ++iter)
    {
        if (hitType == (*iter)->GetHitType())
            ++nHitsOfSpecifiedType;
    }

    return nHitsOfSpecifiedType;
}

} // namespace lar_content
