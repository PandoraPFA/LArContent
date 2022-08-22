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
    m_completenessMode(false),
    m_thresholdCompleteness(0.33f),
    m_thresholdPurity(0.5f),
    m_maxOpeningAngle(5.f),
    m_extensionMode(true),
    m_addInSubdominantHits(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingCCElectronRefinementAlgorithm::Run()
{
    //PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_DEFAULT, -1.f, 1.f, 1.f);

    // Get Shower Pfos
    const PfoList *pShowerPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_showerPfoListName, pShowerPfoList));

    if (!pShowerPfoList || pShowerPfoList->empty())
    {
        std::cout << "CheatingCCElectronRefinementAlgorithm: No shower pfo list found, returning..." << std::endl;
        return STATUS_CODE_FAILURE;
    }

    // Get MC Particles
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    if (!pMCParticleList || pMCParticleList->empty())
    {
        std::cout << "CheatingCCElectronRefinementAlgorithm: No MC particle list found, returning..." << std::endl;
        return STATUS_CODE_FAILURE;
    }

    if (!LArMCParticleHelper::IsCCNuEvent(pMCParticleList, 12))
        return STATUS_CODE_SUCCESS; 

    // Get CaloHits
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    if (!pCaloHitList || pCaloHitList->empty())
    {
        std::cout << "CheatingCCElectronRefinementAlgorithm: No calo hit list found, returning..." << std::endl;
        return STATUS_CODE_FAILURE;
    }

    // Fill electron -> hit map
    LArMCParticleHelper::MCContributionMap electronHitMap;
    this->FillElectronHitMap(pCaloHitList, electronHitMap);

    if (electronHitMap.empty())
        return STATUS_CODE_SUCCESS;

    MCParticleVector mcElectronVector;

    for (const auto &entry : electronHitMap)
        mcElectronVector.push_back(entry.first);

    std::sort(mcElectronVector.begin(), mcElectronVector.end(), LArMCParticleHelper::SortByMomentum);

    for (const MCParticle *const pMCElectron : mcElectronVector)
    {
        const ParticleFlowObject *const pElectronPfo(this->FindBestPfoMatch(pMCElectron, pShowerPfoList, electronHitMap));

        if (!pElectronPfo)
            continue;

        //PfoList displayPfo;
        //displayPfo.push_back(pElectronPfo);
        //PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &displayPfo, "before", RED, true, true);
        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
        
        this->FillOwnershipMaps(pShowerPfoList);

        this->RefineElectronPfos(pElectronPfo, pMCElectron, electronHitMap);

        //PfoList displayPfo2;
        //displayPfo2.push_back(pElectronPfo);
        //PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &displayPfo2, "after", BLUE, true, true);
        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingCCElectronRefinementAlgorithm::FillElectronHitMap(const CaloHitList *const pCaloHitList, LArMCParticleHelper::MCContributionMap &mcElectronToHitMap) const
{
    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        MCParticleVector contributingMCParticleVector;
        const MCParticleWeightMap &weightMap(pCaloHit->GetMCParticleWeightMap());

        for (const auto &mapEntry : weightMap)
        {
            const MCParticle *pMCParticle(mapEntry.first);
            const bool isLeadingElectron((std::abs(pMCParticle->GetParticleId()) == 11) && (LArMCParticleHelper::GetPrimaryMCParticle(pMCParticle) == pMCParticle));

            if (isLeadingElectron)
                contributingMCParticleVector.push_back(pMCParticle);
        }

        std::sort(contributingMCParticleVector.begin(), contributingMCParticleVector.end(), PointerLessThan<MCParticle>());

        float highestWeight(-1.f);
        const MCParticle *highestElectronContributor(nullptr);

        for (const MCParticle *const pMCParticle : contributingMCParticleVector)
        {
            if (weightMap.find(pMCParticle) == weightMap.end())
                continue;

            const float weight(weightMap.at(pMCParticle));

            if (weight > highestWeight)
            {
                highestWeight = weight;
                highestElectronContributor = pMCParticle;
            }
        }

        if (highestElectronContributor)
            mcElectronToHitMap[highestElectronContributor].push_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const ParticleFlowObject *CheatingCCElectronRefinementAlgorithm::FindBestPfoMatch(const MCParticle *const pMCElectron, const PfoList *const pShowerPfoList, 
    LArMCParticleHelper::MCContributionMap &mcElectronToHitMap) const
{
    float bestCompleteness(-1.0), furthestSeparation(-std::numeric_limits<float>::max());
    const ParticleFlowObject *pBestCompletenessPfo(nullptr), *pFurthestPfo(nullptr);

    for (const ParticleFlowObject *const pShowerPfo : *pShowerPfoList)
    {
        CaloHitList pfoHitList;
        ClusterList pfo3DClusterList;

        LArPfoHelper::GetCaloHits(pShowerPfo, TPC_VIEW_U, pfoHitList);
        LArPfoHelper::GetCaloHits(pShowerPfo, TPC_VIEW_V, pfoHitList);
        LArPfoHelper::GetCaloHits(pShowerPfo, TPC_VIEW_W, pfoHitList);
        LArPfoHelper::GetClusters(pShowerPfo, TPC_3D, pfo3DClusterList);

        if (pfoHitList.empty())
            continue;

        const CaloHitList &mcHitList(mcElectronToHitMap.at(pMCElectron));
        const CaloHitList sharedHitList(LArMCParticleHelper::GetSharedHits(pfoHitList, mcHitList));

        if (sharedHitList.empty())
            continue;

        const float completeness(static_cast<float>(sharedHitList.size()) / static_cast<float>(mcHitList.size()));
        const float purity(static_cast<float>(sharedHitList.size()) / static_cast<float>(pfoHitList.size()));

        if (completeness < m_thresholdCompleteness)
            continue;

        if (purity < m_thresholdPurity)
            continue;

        if (completeness > bestCompleteness)
        {
            bestCompleteness = completeness;
            pBestCompletenessPfo = pShowerPfo;
        }

        if (pfo3DClusterList.empty())
            continue;

        const float separation(LArClusterHelper::GetClosestDistance(pMCElectron->GetVertex(), pfo3DClusterList.front()));

        if (separation > furthestSeparation)
        {
            furthestSeparation = separation;
            pFurthestPfo = pShowerPfo;
        }
    }

    if (!pFurthestPfo && pBestCompletenessPfo)
        return pBestCompletenessPfo;

    if (!pBestCompletenessPfo && pFurthestPfo)
        return pFurthestPfo;

    if (m_completenessMode)
        return pBestCompletenessPfo;
    else
        return pFurthestPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingCCElectronRefinementAlgorithm::FillOwnershipMaps(const PfoList *const pShowerPfoList)
{
    m_hitToClusterMapU.clear(); m_hitToClusterMapV.clear(); m_hitToClusterMapW.clear();
    m_clusterToPfoMapU.clear(); m_clusterToPfoMapV.clear(); m_clusterToPfoMapW.clear();

    const PfoList *pTrackPfoList(nullptr);
    PandoraContentApi::GetList(*this, m_trackPfoListName, pTrackPfoList);

    PfoList combinedPfoList(*pShowerPfoList);
    if (!(!pTrackPfoList || pTrackPfoList->empty()))
    {
        combinedPfoList.insert(combinedPfoList.end(), pTrackPfoList->begin(), pTrackPfoList->end());
    }

    PfoList allPfoList;
    LArPfoHelper::GetAllConnectedPfos(combinedPfoList, allPfoList);

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
                    caloHitList.push_back(pIsolated);
            }

            const HitType hitType(caloHitList.front()->GetHitType());
            HitToClusterMap &hitToClusterMap(hitType == TPC_VIEW_U ? m_hitToClusterMapU : hitType == TPC_VIEW_V ? m_hitToClusterMapV : m_hitToClusterMapW);
            ClusterToPfoMap &clusterToPfoMap(hitType == TPC_VIEW_U ? m_clusterToPfoMapU : hitType == TPC_VIEW_V ? m_clusterToPfoMapV : m_clusterToPfoMapW);

            for (const CaloHit *const pCaloHit : caloHitList)
                hitToClusterMap[pCaloHit] = pCluster;

            clusterToPfoMap[pCluster] = pPfo;
        }
    }

    ClusterList clusterList;
    for (const std::string &clusterListName : m_clusterListNames)
    {
        const ClusterList *pClusterList(nullptr);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));

        if (!pClusterList || pClusterList->empty())
        {
            std::cout << "CheatingCCElectronRefinementAlgorithm: No cluster list found, returning..." << std::endl;
            throw;
        }

        clusterList.insert(clusterList.end(), pClusterList->begin(), pClusterList->end());
    }

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
            throw;
        }

        CaloHitList caloHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

        CaloHitList isolated(pCluster->GetIsolatedCaloHitList());

        for (const CaloHit *const pIsolated : isolated)
        {
            if (std::find(caloHitList.begin(), caloHitList.end(), pIsolated) == caloHitList.end())
                caloHitList.push_back(pIsolated);
        }

        for (const CaloHit *const pCaloHit : caloHitList)
            hitToClusterMap[pCaloHit] = pCluster;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingCCElectronRefinementAlgorithm::RefineElectronPfos(const ParticleFlowObject *const pElectronPfo, const MCParticle *const pMCElectron, 
    LArMCParticleHelper::MCContributionMap &mcElectronToHitMap)
{
    const CaloHitList &caloHitList(mcElectronToHitMap.at(pMCElectron));

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

    const float coneLengthU(this->FindConeLength(pElectronPfo, mcVertexU, mcDirectionU, TPC_VIEW_U));
    const float coneLengthV(this->FindConeLength(pElectronPfo, mcVertexV, mcDirectionV, TPC_VIEW_V));
    const float coneLengthW(this->FindConeLength(pElectronPfo, mcVertexW, mcDirectionW, TPC_VIEW_W));

    //PandoraMonitoringApi::ViewEvent(this->GetPandora());

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        // Only add in dominant hits (if this is true)
        if (!m_addInSubdominantHits)
        {
            try
            {
                const MCParticle *pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));

                if (pMCParticle != pMCElectron)
                    continue;
            }
            catch (...)
            {
                continue;
            }
        }

        const HitType hitType(pCaloHit->GetHitType());
        const CartesianVector viewMCVertex(hitType == TPC_VIEW_U ? mcVertexU : hitType == TPC_VIEW_V ? mcVertexV : mcVertexW);
        const CartesianVector viewMCDirection(hitType == TPC_VIEW_U ? mcDirectionU : hitType == TPC_VIEW_V ? mcDirectionV : mcDirectionW);
        const float viewConeLength(hitType == TPC_VIEW_U ? coneLengthU : hitType == TPC_VIEW_V ? coneLengthV : coneLengthW);

        if (!this->DoesPassCut(pCaloHit, viewMCVertex, viewMCDirection, viewConeLength))
            continue;

        if (m_extensionMode)
        {
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

            if (pParentPfo == pElectronPfo)
                continue;

            //CartesianVector hitPosition(pCaloHit->GetPositionVector());
            //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "hit to add", VIOLET, 2);

            if (pParentCluster)
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, clusterListName));

                CaloHitList clusterNormalHitList;
                pParentCluster->GetOrderedCaloHitList().FillCaloHitList(clusterNormalHitList);
                const CaloHitList clusterIsolatedHitList(pParentCluster->GetIsolatedCaloHitList());
                const bool isIsolated(std::find(clusterIsolatedHitList.begin(), clusterIsolatedHitList.end(), pCaloHit) != clusterIsolatedHitList.end());

                if (!isIsolated && (clusterNormalHitList.size() == 1) && !(clusterIsolatedHitList.empty()))
                {
                    const HitType isolatedHitType(LArClusterHelper::GetClusterHitType(pParentCluster));
                    HitToClusterMap &isolatedHitToClusterMap(isolatedHitType == TPC_VIEW_U ? m_hitToClusterMapU : hitType == TPC_VIEW_V ? m_hitToClusterMapV : m_hitToClusterMapW);

                    for (const CaloHit * const pIsolatedHit : clusterIsolatedHitList)
                    {
                        isolatedHitToClusterMap.erase(pIsolatedHit);
                        const StatusCode isolatedStatusCode(PandoraContentApi::RemoveIsolatedFromCluster(*this, pParentCluster, pIsolatedHit));

                        if (isolatedStatusCode != STATUS_CODE_SUCCESS)
                        {
                            std::cout << "ISOBEL CANNOT REMOVE ISOLATED HIT?" << std::endl;
                            throw;
                        }
                    }
                }

                const StatusCode statusCodeCluster(isIsolated ? PandoraContentApi::RemoveIsolatedFromCluster(*this, pParentCluster, pCaloHit) : 
                    PandoraContentApi::RemoveFromCluster(*this, pParentCluster, pCaloHit));

                if (statusCodeCluster != STATUS_CODE_SUCCESS)
                {
                    if (statusCodeCluster != STATUS_CODE_NOT_ALLOWED)
                    {
                        std::cout << "CheatingCCElectronRefinementAlgorithm: cluster jam" << std::endl;
                        throw StatusCodeException(statusCodeCluster);
                    }

                    if (pParentPfo)
                    {
                        const StatusCode statusCodePfo(PandoraContentApi::RemoveFromPfo(*this, pParentPfo, pParentCluster));
                        const unsigned int nHits(LArPfoHelper::GetNumberOfTwoDHits(pParentPfo));

                        if (nHits == 0)
                            std::cout << "CheatingCCElectronRefinementAlgorithm: ISOBEL - PFO HAS ZERO HITS" << std::endl;

                        if (statusCodePfo != STATUS_CODE_SUCCESS)
                        {
                            std::cout << "CheatingCCElectronRefinementAlgorithm: pfo jam" << std::endl;
                            throw;
                        }
                    }

                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pParentCluster));;
                }
            }

            if (!PandoraContentApi::IsAvailable(*this, pCaloHit))
            {
                std::cout << "CALO HIT IS NOT AVAILABLE!!" << std::endl;
                throw;
            }

            ClusterList electronClusters;
            LArPfoHelper::GetClusters(pElectronPfo, hitType, electronClusters);

            if (electronClusters.empty())
            {
                const Cluster *pNewCluster(this->CreateCluster(pCaloHit, hitType));

                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,  PandoraContentApi::AddToPfo(*this, pElectronPfo, pNewCluster));
            }
            else
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, electronClusters.front(), pCaloHit));
            }
        }
    }

    //PandoraMonitoringApi::ViewEvent(this->GetPandora());
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

const Cluster *CheatingCCElectronRefinementAlgorithm::CreateCluster(const CaloHit *const pCaloHit, const HitType hitType) const
{
    const ClusterList *pTemporaryList(nullptr);
    std::string temporaryListName, currentListName;

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*this, currentListName));

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        PandoraContentApi::CreateTemporaryListAndSetCurrent<ClusterList>(*this, pTemporaryList, temporaryListName));

    const Cluster *pCluster(nullptr);
    PandoraContentApi::Cluster::Parameters parameters;
    parameters.m_caloHitList.push_back(pCaloHit);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pCluster));

    PandoraContentApi::Cluster::Metadata metadata;
    metadata.m_particleId = 11;

    if (metadata.m_particleId.IsInitialized())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::AlterMetadata(*this, pCluster, metadata));

    const std::string &clusterListName(hitType == TPC_VIEW_U ? "ClustersU" : hitType == TPC_VIEW_V ? "ClustersV" : "ClustersW");

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, temporaryListName, clusterListName));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, currentListName));

    return pCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingCCElectronRefinementAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ShowerPfoListName", m_showerPfoListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrackPfoListName", m_trackPfoListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadVectorOfValues(xmlHandle, "ClusterListNames", m_clusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "CompletenessMode", m_completenessMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "ThresholdCompleteness", m_thresholdCompleteness));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "ThresholdPurity", m_thresholdPurity));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "MaxOpeningAngle", m_maxOpeningAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "ExtensionMode", m_extensionMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "AddInSubdominantHits", m_addInSubdominantHits));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
