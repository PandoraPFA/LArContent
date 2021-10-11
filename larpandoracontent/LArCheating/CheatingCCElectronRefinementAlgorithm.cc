/**
 *  @file   larpandoracontent/LArCheating/CheatingCCElectronRefinementAlgorithm.cc
 *
 *  @brief  Implementation of the cheating vertex creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArCheating/CheatingCCElectronRefinementAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CheatingCCElectronRefinementAlgorithm::CheatingCCElectronRefinementAlgorithm() : 
    m_maxOpeningAngle(5.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingCCElectronRefinementAlgorithm::Run()
{
    //PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_DEFAULT, -1.f, 1.f, 1.f);

    PfoList allPfoList; ClusterList clusterList;
    const MCParticleList *pMCParticleList(nullptr); const PfoList *pShowerPfoList(nullptr); const CaloHitList *pCaloHitList(nullptr);

    if (this->FillLists(pMCParticleList, pShowerPfoList, allPfoList, pCaloHitList, clusterList) != STATUS_CODE_SUCCESS)
        return STATUS_CODE_SUCCESS;

    if (!LArMCParticleHelper::IsCCNuEvent(pMCParticleList, 12))
    {
        //std::cout << "This is not a nue CC event" << std::endl;
        return STATUS_CODE_SUCCESS;
    }

    std::cout << "this is a nue CC event!" << std::endl;

    MCParticleToPfoMap mcElectronToPfoMap;
    this->FindElectronToPfoMatches(pShowerPfoList, mcElectronToPfoMap);

    if (mcElectronToPfoMap.empty())
        return STATUS_CODE_SUCCESS;

    ////////////////////////
    /*
    PfoList pfolist;
    for (const auto &entry : mcElectronToPfoMap)
        pfolist.insert(pfolist.begin(), entry.second.begin(), entry.second.end());

    PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &pfolist, "matched electron pfos", RED, true, true);
    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    */
    ////////////////////////

    LArMCParticleHelper::MCContributionMap mcElectronToHitMap;
    this->FindMissingElectronHits(pCaloHitList, mcElectronToPfoMap, mcElectronToHitMap);

    /////////////////////
    /*
    CaloHitList caloHitList;
    for (const auto &entry : mcElectronToHitMap)
        caloHitList.insert(caloHitList.end(), entry.second.begin(), entry.second.end());

    PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &caloHitList, "true electron hits", BLUE);
    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    */
    ////////////////////

    this->FindBestElectronToPfoMatch(mcElectronToPfoMap, mcElectronToHitMap);

    if (mcElectronToPfoMap.empty())
        return STATUS_CODE_SUCCESS;

    ////////////////////////
    /*
    PfoList bestpfolist;
    for (const auto &entry : mcElectronToPfoMap)
        bestpfolist.insert(bestpfolist.begin(), entry.second.begin(), entry.second.end());

    PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &bestpfolist, "BEST matched electron pfo", RED, true, true);
    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    
    ////////////////////////
    ////////////////////////
    for (const auto &entry : mcElectronToPfoMap)
    {
        if (entry.second.empty())
            continue;

        const CartesianVector &mcVertex(entry.first->GetVertex());
        const CartesianVector mcVertexU(LArGeometryHelper::ProjectPosition(this->GetPandora(), mcVertex, TPC_VIEW_U));
        const CartesianVector mcVertexV(LArGeometryHelper::ProjectPosition(this->GetPandora(), mcVertex, TPC_VIEW_V));
        const CartesianVector mcVertexW(LArGeometryHelper::ProjectPosition(this->GetPandora(), mcVertex, TPC_VIEW_W));

        const ParticleFlowObject *const pPfo(entry.second.front());
        CaloHitList caloHitListU, caloHitListV, caloHitListW;
        LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_U, caloHitListU);
        LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_V, caloHitListV);
        LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, caloHitListW);

        float closestDistanceU(std::numeric_limits<float>::max());
        float closestDistanceV(std::numeric_limits<float>::max());
        float closestDistanceW(std::numeric_limits<float>::max());

        for (const CaloHit *const pCaloHit : caloHitListU)
        {
            const CartesianVector &position(pCaloHit->GetPositionVector());
            const float separation((position - mcVertexU).GetMagnitude());

            if (separation < closestDistanceU)
                closestDistanceU = separation;
        }

        for (const CaloHit *const pCaloHit : caloHitListV)
        {
            const CartesianVector &position(pCaloHit->GetPositionVector());
            const float separation((position - mcVertexV).GetMagnitude());

            if (separation < closestDistanceV)
                closestDistanceV = separation;
        }

        for (const CaloHit *const pCaloHit : caloHitListW)
        {
            const CartesianVector &position(pCaloHit->GetPositionVector());
            const float separation((position - mcVertexW).GetMagnitude());

            if (separation < closestDistanceW)
                closestDistanceW = separation;
        }


        std::cout << "closestDistanceU: " << closestDistanceU << std::endl;
        std::cout << "closestDistanceV: " << closestDistanceV << std::endl;
        std::cout << "closestDistanceW: " << closestDistanceW << std::endl;
    }
    */
    ///////////////////////////////////////////

    this->FillOwnershipMaps(allPfoList, clusterList);

    this->RefineElectronPfos(mcElectronToHitMap, mcElectronToPfoMap);

    return STATUS_CODE_SUCCESS;
}


//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingCCElectronRefinementAlgorithm::FindElectronToPfoMatches(const PfoList *const pShowerPfoList, MCParticleToPfoMap &mcElectronToPfoMap) const
{
    for (const ParticleFlowObject *const pPfo : *pShowerPfoList)
    {
        const MCParticle *pMCParticle(LArMCParticleHelper::GetMainMCParticle(pPfo));
        const bool isLeadingElectron((std::abs(pMCParticle->GetParticleId()) == 11) && (LArMCParticleHelper::GetPrimaryMCParticle(pMCParticle) == pMCParticle));

        if (isLeadingElectron)
            mcElectronToPfoMap[pMCParticle].push_back(pPfo);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingCCElectronRefinementAlgorithm::FindMissingElectronHits(const CaloHitList *const pCaloHitList, const MCParticleToPfoMap &mcElectronToPfoMap, 
    LArMCParticleHelper::MCContributionMap &mcElectronToHitMap) const
{
    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        MCParticleVector contributingMCParticleVector;
        const MCParticleWeightMap &weightMap(pCaloHit->GetMCParticleWeightMap());

        for (const auto &mapEntry : weightMap)
            contributingMCParticleVector.push_back(mapEntry.first);

        std::sort(contributingMCParticleVector.begin(), contributingMCParticleVector.end(), PointerLessThan<MCParticle>());

        float highestWeight(0.f);
        const MCParticle *highestElectronContributor(nullptr);

        for (const MCParticle *const pMCParticle : contributingMCParticleVector)
        {
            const float weight(weightMap.at(pMCParticle));

            if (mcElectronToPfoMap.find(pMCParticle) != mcElectronToPfoMap.end())
            {
                if (weight > highestWeight)
                {
                    highestWeight = weight;
                    highestElectronContributor = pMCParticle;
                }
            }
        }

        if (highestElectronContributor)
            mcElectronToHitMap[highestElectronContributor].push_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingCCElectronRefinementAlgorithm::FindBestElectronToPfoMatch(MCParticleToPfoMap &mcElectronToPfoMap, LArMCParticleHelper::MCContributionMap &mcElectronToHitMap) const
{
    MCParticleToPfoMap tempMCElectronToPfoMap(mcElectronToPfoMap);

    mcElectronToPfoMap.clear();

    MCParticleVector mcElectronList;

    for (const auto &mapEntry : tempMCElectronToPfoMap)
        mcElectronList.push_back(mapEntry.first);

    std::sort(mcElectronList.begin(), mcElectronList.end(), PointerLessThan<MCParticle>());

    for (const MCParticle *const pMCElectron : mcElectronList)
    {
        const CaloHitList &mcElectronHits(mcElectronToHitMap.at(pMCElectron));
        const PfoList &matchedPfoList(tempMCElectronToPfoMap.at(pMCElectron));
        PfoVector matchedPfoVector(matchedPfoList.begin(), matchedPfoList.end());

        std::sort(matchedPfoVector.begin(), matchedPfoVector.end(), LArPfoHelper::SortByNHits);

        float bestCompleteness(-1.0);
        const ParticleFlowObject *pBestMatchedPfo(nullptr);

        for (const ParticleFlowObject *const pMatchedPfo : matchedPfoVector)
        {
            CaloHitList matchedPfoHitList;

            LArPfoHelper::GetCaloHits(pMatchedPfo, TPC_VIEW_U, matchedPfoHitList);
            LArPfoHelper::GetCaloHits(pMatchedPfo, TPC_VIEW_V, matchedPfoHitList);
            LArPfoHelper::GetCaloHits(pMatchedPfo, TPC_VIEW_W, matchedPfoHitList);

            const CaloHitList sharedHitList(LArMCParticleHelper::GetSharedHits(matchedPfoHitList, mcElectronHits));
            const float completeness(static_cast<float>(sharedHitList.size()) / static_cast<float>(mcElectronHits.size()));
            const float purity(static_cast<float>(sharedHitList.size()) / static_cast<float>(matchedPfoHitList.size()));

            if (completeness < 0.33)
                continue;

            if (purity < 0.5)
                continue;

            if (completeness > bestCompleteness)
            {
                bestCompleteness = completeness;
                pBestMatchedPfo = pMatchedPfo;
            }
        }

        if (pBestMatchedPfo)
            mcElectronToPfoMap[pMCElectron].push_back(pBestMatchedPfo);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingCCElectronRefinementAlgorithm::FillOwnershipMaps(const PfoList &allPfoList, const ClusterList &clusterList)
{
    for (const ParticleFlowObject *const pPfo : allPfoList)
    {
        ClusterList twoDClusterList;
        LArPfoHelper::GetClusters(pPfo, TPC_VIEW_U, twoDClusterList);
        LArPfoHelper::GetClusters(pPfo, TPC_VIEW_V, twoDClusterList);
        LArPfoHelper::GetClusters(pPfo, TPC_VIEW_W, twoDClusterList);

        for (const Cluster *const pCluster : twoDClusterList)
        {
            CaloHitList caloHitList;
            pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

            CaloHitList isolated(pCluster->GetIsolatedCaloHitList());

            for (const CaloHit *const pIsolated : isolated)
            {
                if (std::find(caloHitList.begin(), caloHitList.end(), pIsolated) == caloHitList.end())
                {
                    std::cout << "adding isolated hit..." << std::endl;
                    caloHitList.push_back(pIsolated);
                }
            }

            const HitType hitType(caloHitList.front()->GetHitType());
            HitToClusterMap &hitToClusterMap(hitType == TPC_VIEW_U ? m_hitToClusterMapU : hitType == TPC_VIEW_V ? m_hitToClusterMapV : m_hitToClusterMapW);
            ClusterToPfoMap &clusterToPfoMap(hitType == TPC_VIEW_U ? m_clusterToPfoMapU : hitType == TPC_VIEW_V ? m_clusterToPfoMapV : m_clusterToPfoMapW);

            for (const CaloHit *const pCaloHit : caloHitList)
                hitToClusterMap[pCaloHit] = pCluster;

            clusterToPfoMap[pCluster] = pPfo;
        }
    }

    //ClusterList notAvailable;
    for (const Cluster *const pCluster : clusterList)
    {
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
        HitToClusterMap &hitToClusterMap(hitType == TPC_VIEW_U ? m_hitToClusterMapU : hitType == TPC_VIEW_V ? m_hitToClusterMapV : m_hitToClusterMapW);
        ClusterToPfoMap &clusterToPfoMap(hitType == TPC_VIEW_U ? m_clusterToPfoMapU : hitType == TPC_VIEW_V ? m_clusterToPfoMapV : m_clusterToPfoMapW);

        if (clusterToPfoMap.find(pCluster) != clusterToPfoMap.end())
            continue;

        if (!pCluster->IsAvailable())
        {
            std::cout << "CLUSTER IS NOT AVAILABLE ISOBE" << std::endl;
            //notAvailable.push_back(pCluster);
        }

        CaloHitList caloHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

        CaloHitList isolated(pCluster->GetIsolatedCaloHitList());

        for (const CaloHit *const pIsolated : isolated)
        {
            if (std::find(caloHitList.begin(), caloHitList.end(), pIsolated) == caloHitList.end())
            {
                std::cout << "adding isolated hit..." << std::endl;
                caloHitList.push_back(pIsolated);
            }
        }

        for (const CaloHit *const pCaloHit : caloHitList)
            hitToClusterMap[pCaloHit] = pCluster;
    }

    //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &notAvailable, "not available", BLACK);
    //PandoraMonitoringApi::ViewEvent(this->GetPandora());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingCCElectronRefinementAlgorithm::RefineElectronPfos(LArMCParticleHelper::MCContributionMap &mcElectronToHitMap, 
    MCParticleToPfoMap &mcElectronToPfoMap) const
{
    MCParticleVector mcElectronList;

    for (auto &entry : mcElectronToHitMap)
        mcElectronList.push_back(entry.first);

    std::sort(mcElectronList.begin(), mcElectronList.end(), LArMCParticleHelper::SortByMomentum);

    for (const MCParticle *const pMCElectron : mcElectronList)
    {
        const CaloHitList &caloHitList(mcElectronToHitMap.at(pMCElectron));
        const ParticleFlowObject *const pMatchedPfo(mcElectronToPfoMap.at(pMCElectron).front());

        // Parameterise extension cone
        const CartesianVector &mcVertex(pMCElectron->GetVertex());
        const CartesianVector mcVertexU(LArGeometryHelper::ProjectPosition(this->GetPandora(), mcVertex, TPC_VIEW_U));
        const CartesianVector mcVertexV(LArGeometryHelper::ProjectPosition(this->GetPandora(), mcVertex, TPC_VIEW_V));
        const CartesianVector mcVertexW(LArGeometryHelper::ProjectPosition(this->GetPandora(), mcVertex, TPC_VIEW_W));

        CartesianVector mcDirection(pMCElectron->GetMomentum());

        if (mcDirection.GetMagnitude() < std::numeric_limits<float>::epsilon())
            return;

        mcDirection = mcDirection * (1.f / mcDirection.GetMagnitude());
        const CartesianVector mcDirectionU(LArGeometryHelper::ProjectPosition(this->GetPandora(), mcDirection, TPC_VIEW_U));
        const CartesianVector mcDirectionV(LArGeometryHelper::ProjectPosition(this->GetPandora(), mcDirection, TPC_VIEW_V));
        const CartesianVector mcDirectionW(LArGeometryHelper::ProjectPosition(this->GetPandora(), mcDirection, TPC_VIEW_W));

        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &mcVertexU, "mcVertexU", BLACK, 1);
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &mcVertexV, "mcVertexV", BLACK, 1);
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &mcVertexW, "mcVertexW", BLACK, 1);

        const float coneLengthU(this->FindConeLength(pMatchedPfo, mcVertexU, mcDirectionU, TPC_VIEW_U));
        const float coneLengthV(this->FindConeLength(pMatchedPfo, mcVertexV, mcDirectionV, TPC_VIEW_V));
        const float coneLengthW(this->FindConeLength(pMatchedPfo, mcVertexW, mcDirectionW, TPC_VIEW_W));

        //PandoraMonitoringApi::ViewEvent(this->GetPandora());

        for (const CaloHit *const pCaloHit : caloHitList)
        {
            const HitType hitType(pCaloHit->GetHitType());
            const CartesianVector viewMCVertex(hitType == TPC_VIEW_U ? mcVertexU : hitType == TPC_VIEW_V ? mcVertexV : mcVertexW);
            const CartesianVector viewMCDirection(hitType == TPC_VIEW_U ? mcDirectionU : hitType == TPC_VIEW_V ? mcDirectionV : mcDirectionW);
            const float viewConeLength(hitType == TPC_VIEW_U ? coneLengthU : hitType == TPC_VIEW_V ? coneLengthV : coneLengthW);

            if (!this->DoesPassCut(pCaloHit, viewMCVertex, viewMCDirection, viewConeLength))
                continue;

            // Now remove...
            const HitToClusterMap &hitToClusterMap(hitType == TPC_VIEW_U ? m_hitToClusterMapU : hitType == TPC_VIEW_V ? m_hitToClusterMapV : m_hitToClusterMapW);
            const ClusterToPfoMap &clusterToPfoMap(hitType == TPC_VIEW_U ? m_clusterToPfoMapU : hitType == TPC_VIEW_V ? m_clusterToPfoMapV : m_clusterToPfoMapW);
            std::string clusterListName(hitType == TPC_VIEW_U ? "ClustersU" : hitType == TPC_VIEW_V ? "ClustersV" : "ClustersW");

            const Cluster *pParentCluster(nullptr);
            const ParticleFlowObject *pParentPfo(nullptr);

            if (hitToClusterMap.find(pCaloHit) != hitToClusterMap.end())
            {
                pParentCluster = hitToClusterMap.at(pCaloHit);

                if (clusterToPfoMap.find(pParentCluster) != clusterToPfoMap.end())
                    pParentPfo = clusterToPfoMap.at(pParentCluster);
            }

            if (pParentPfo == pMatchedPfo)
                continue;

            //CartesianVector hitPosition(pCaloHit->GetPositionVector());
            //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "hit to add", VIOLET, 2);
            //std::cout << "AAAAAAA" << std::endl;
            if (pParentCluster)
            {
                //std::cout << "BBBBBBBBBBBBBB" << std::endl;
                const StatusCode statusCodeCluster(PandoraContentApi::RemoveFromCluster(*this, pParentCluster, pCaloHit));

                if (statusCodeCluster != STATUS_CODE_SUCCESS)
                {
                    //std::cout << "CCCCCCCCCCCC" << std::endl;
                    if (statusCodeCluster != STATUS_CODE_NOT_ALLOWED)
                    {
                        std::cout << "jam" << std::endl;
                        throw StatusCodeException(statusCodeCluster);
                    }

                    if (pParentPfo)
                    {
                        const StatusCode statusCodePfo(PandoraContentApi::RemoveFromPfo(*this, pParentPfo, pParentCluster));

                        unsigned int nHits(LArPfoHelper::GetNumberOfTwoDHits(pParentPfo));

                        if (nHits == 0)
                            std::cout << "CheatingCCElectronRefinementAlgorithm: ISOBEL - PFO HAS ZERO HITS" << std::endl;

                        if (statusCodePfo != STATUS_CODE_SUCCESS)
                        {
                            std::cout << "pfo jam" << std::endl;
                            throw;
                            // work out how many clusters pfo has?? if none...
                            //{
                            //PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pParentPfo));
                            //}
                        }
                    }

                    //std::cout << "EEEEEEEEEE" << std::endl;
                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, clusterListName));
                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pParentCluster));
                }
            }

            if (!PandoraContentApi::IsAvailable(*this, pCaloHit))
                std::cout << "CALO HIT IS NOT AVAILABLE!!" << std::endl;

            ClusterList matchedPfoCluster;
            LArPfoHelper::GetClusters(pMatchedPfo, hitType, matchedPfoCluster);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, matchedPfoCluster.front(), pCaloHit));
        }

        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------  

float CheatingCCElectronRefinementAlgorithm::FindConeLength(const ParticleFlowObject *const pPfo, const CartesianVector &mcVertex, 
    const CartesianVector &mcDirection, const HitType hitType) const
{
    CaloHitList caloHitList;
    LArPfoHelper::GetCaloHits(pPfo, hitType, caloHitList);

    CartesianVector closestPosition(0.f,0.f,0.f);
    float closestImpactL(std::numeric_limits<float>::max());

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const CartesianVector &position(pCaloHit->GetPositionVector());
        const CartesianVector displacement(position - mcVertex);
        const float impactL(mcDirection.GetDotProduct(displacement));

        if (impactL < 0.f)
            continue;

        if (impactL < closestImpactL)
        {
            closestPosition = pCaloHit->GetPositionVector();
            closestImpactL = impactL;
        }
    }

    //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &closestPosition, "closest impact parameter hit", RED, 1);

    return closestImpactL;
}

//------------------------------------------------------------------------------------------------------------------------------------------  

bool CheatingCCElectronRefinementAlgorithm::DoesPassCut(const CaloHit *const pCaloHit, const CartesianVector &mcVertex, 
    const CartesianVector &mcDirection, const float coneLength) const
{
    const CartesianVector &position(pCaloHit->GetPositionVector());
    const CartesianVector displacement(position - mcVertex);
    const float impactL(mcDirection.GetDotProduct(displacement));

    if (impactL < 0.f)
        return false;

    if (impactL > coneLength)
        return false;

    if (displacement.GetOpeningAngle(mcDirection) > (m_maxOpeningAngle * 3.14 / 180))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingCCElectronRefinementAlgorithm::FillLists(const MCParticleList *&pMCParticleList, const PfoList *&pShowerPfoList, PfoList &allPfoList, 
    const CaloHitList *&pCaloHitList, ClusterList &clusterList) const
{
    // Get MC Particles
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    if (!pMCParticleList || pMCParticleList->empty())
    {
        std::cout << "CheatingCCElectronRefinementAlgorithm: No MC particle list found, returning..." << std::endl;
        return STATUS_CODE_FAILURE;
    }

    // Get all Pfos
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_showerPfoListName, pShowerPfoList));

    if (!pShowerPfoList || pShowerPfoList->empty())
    {
        std::cout << "CheatingCCElectronRefinementAlgorithm: No shower pfo list found, returning..." << std::endl;
        return STATUS_CODE_FAILURE;
    }

    const PfoList *pTrackPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "TrackParticles3D", pTrackPfoList));

    PfoList combinedPfoList(*pShowerPfoList);
    if (!(!pTrackPfoList || pTrackPfoList->empty()))
    {
        combinedPfoList.insert(combinedPfoList.end(), pTrackPfoList->begin(), pTrackPfoList->end());
    }

    LArPfoHelper::GetAllConnectedPfos(combinedPfoList, allPfoList);

    // Get CaloHits
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    if (!pCaloHitList || pCaloHitList->empty())
    {
        std::cout << "CheatingCCElectronRefinementAlgorithm: No calo hit list found, returning..." << std::endl;
        return STATUS_CODE_FAILURE;
    }

    for (const std::string &clusterListName : m_clusterListNames)
    {
        const ClusterList *pClusterList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));

        if (!pClusterList || pClusterList->empty())
        {
            std::cout << "CheatingCCElectronRefinementAlgorithm: No cluster list found, returning..." << std::endl;
            return STATUS_CODE_FAILURE;
        }

        clusterList.insert(clusterList.end(), pClusterList->begin(), pClusterList->end());
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingCCElectronRefinementAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_showerPfoListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadVectorOfValues(xmlHandle, "ClusterListNames", m_clusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "MaxOpeningAngle", m_maxOpeningAngle));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
