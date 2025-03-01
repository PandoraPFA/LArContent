/**
 *  @file   larpandoracontent/LArReclustering/ThreeDMultiReclusteringAlgorithm.cc
 *
 *  @brief  Implementation file for the reclustering algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArReclustering/ThreeDMultiReclusteringAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ThreeDMultiReclusteringAlgorithm::ThreeDMultiReclusteringAlgorithm() :
    m_pFomAlgTool {nullptr},
    m_clusterUListName {"ClustersU"},
    m_clusterVListName {"ClustersV"},
    m_clusterWListName {"ClustersW"}
{ }

//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDMultiReclusteringAlgorithm::~ThreeDMultiReclusteringAlgorithm()
{ }

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDMultiReclusteringAlgorithm::Run()
{
    const PfoList *pPfos {nullptr};
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=,
                                    PandoraContentApi::GetList(*this, m_pfoListName, pPfos));
    if (!pPfos || pPfos->empty())
        return STATUS_CODE_SUCCESS;

    // Ask the FOM alg tool which pfos to recluster
    PfoList pfosToRecluster;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, m_pFomAlgTool->GetPfosToRecluster(pPfos, pfosToRecluster));
    if (pfosToRecluster.empty())
        return STATUS_CODE_SUCCESS;

    // Remove clusters from the pfos, keeping track of their original pfo
    std::map<HitType, ClusterList> viewToFreeClusters = { {TPC_3D, {}}, {TPC_VIEW_U, {}}, {TPC_VIEW_V, {}}, {TPC_VIEW_W, {}} };
    std::map<const Pfo *const, ClusterList> pfoToFreeClusters;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, FreeClustersFromPfos(pfosToRecluster, viewToFreeClusters, pfoToFreeClusters));
    if (viewToFreeClusters.at(TPC_3D).empty())
        return STATUS_CODE_SUCCESS;

    // Run reclustering algs over the free 3D clusters to see if there is a superior 3D clustering
    // NOTE Temporarily replacing the current list doesn't seem to work
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_cluster3DListName));
    std::string originalClusterListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                             PandoraContentApi::InitializeReclustering(*this, TrackList(), viewToFreeClusters.at(TPC_3D), originalClusterListName));

    float bestReclusterFom;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, m_pFomAlgTool->CalcClusteringFom(viewToFreeClusters.at(TPC_3D), bestReclusterFom));
    std::string bestReclusterListName {originalClusterListName};
    ClusterList bestReclusterList;

    for (const std::string &clusteringAlg : m_clusteringAlgs)
    {
        std::string reclusterListName;
        const ClusterList *pReclusterList {nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                                 PandoraContentApi::RunClusteringAlgorithm(*this, clusteringAlg, pReclusterList, reclusterListName));
        if (!pReclusterList || pReclusterList->empty())
            continue;

        float reclusterFom;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, m_pFomAlgTool->CalcClusteringFom(*pReclusterList, reclusterFom));
        if (reclusterFom > bestReclusterFom)
        {
            bestReclusterFom = reclusterFom;
            bestReclusterListName = reclusterListName;
            bestReclusterList = *pReclusterList;
        }
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::TemporarilyReplaceCurrentList<Cluster>(*this, bestReclusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndReclustering(*this, bestReclusterListName));

    if (bestReclusterListName == originalClusterListName) // Return the clusters to the original pfos as if nothing happened...
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, AddClustersToPfos(pfoToFreeClusters));
        return STATUS_CODE_SUCCESS;
    }
    else
    {
        // Delete the original pfos in preparation for making new ones
        // ATTN This will break any existing particle hierarchy
        for (const Pfo *const pPfo : pfosToRecluster)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pPfo, m_pfoListName));
    }

    // The 3D clusters are reclustered, now need to manually recluster the 2D clusters to be consistent with the 3D reclustering
    // First, remove the original 2D clusters
    std::map<HitType, CaloHitList> viewToFreeCaloHits2D = { {TPC_VIEW_U, {}}, {TPC_VIEW_V, {}}, {TPC_VIEW_W, {}} };
    for (const auto &[view, freeClusters] : viewToFreeClusters)
    {
        if (view == TPC_3D)
            continue;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, FreeCaloHitsFromClusters(freeClusters, view, viewToFreeCaloHits2D.at(view)));
    }

    // Now make new 2D clusters following the new 3D clusters, and make new corresponding pfos
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, BuildNewPfos(bestReclusterList, viewToFreeCaloHits2D));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDMultiReclusteringAlgorithm::FreeClustersFromPfos(const PfoList &pfos,
                                                                  std::map<HitType, ClusterList> &viewToFreeClusters,
                                                                  std::map<const Pfo *const, ClusterList> &pfoToFreeClusters)
{
    for (const Pfo *const pPfo : pfos)
    {
        ClusterList thisPfoClusters;
        for (auto &[view, freeClusters] : viewToFreeClusters)
        {
            ClusterList clusters;
            LArPfoHelper::GetClusters(pPfo, view, clusters);
            for (const Cluster *const pCluster : clusters)
            {
                freeClusters.push_back(pCluster);
                thisPfoClusters.push_back(pCluster);
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromPfo(*this, pPfo, pCluster));
            }
        }
        pfoToFreeClusters.insert(std::make_pair(pPfo, thisPfoClusters));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDMultiReclusteringAlgorithm::AddClustersToPfos(std::map<const Pfo *const, ClusterList> &pfoToClusters)
{
    for (const auto &[pPfo, clusters] : pfoToClusters)
        for (const Cluster *const pCluster : clusters)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pPfo, pCluster));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDMultiReclusteringAlgorithm::FreeCaloHitsFromClusters(const ClusterList &clusters,
                                                                      const HitType &view,
                                                                      CaloHitList &freeCaloHits)
{
    std::string clusterListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, GetClusterListName(view, clusterListName));
    for (const Cluster *const pCluster: clusters)
    {
        CaloHitList caloHits;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHits);
        freeCaloHits.insert(freeCaloHits.end(), caloHits.begin(), caloHits.end());
        freeCaloHits.insert(freeCaloHits.end(), pCluster->GetIsolatedCaloHitList().begin(), pCluster->GetIsolatedCaloHitList().end());
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pCluster, clusterListName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDMultiReclusteringAlgorithm::BuildNewPfos(const ClusterList &clusters, std::map<HitType, CaloHitList> &viewToFreeCaloHits2D)
{
    // Cluster any 2D hits directly associated with 3D hits in the same way as the 3D clusters
    std::map<HitType, ClusterList> viewToNewClusters2D = { {TPC_VIEW_U, {}}, {TPC_VIEW_V, {}}, {TPC_VIEW_W, {}} };
    std::vector<ClusterList> newPfoClusters;
    for (const Cluster *const pCluster3D : clusters)
    {
        ClusterList newClusters;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, Build2DClustersFrom3D(pCluster3D, viewToFreeCaloHits2D, newClusters));
        for (const Cluster *const pCluster : newClusters)
            viewToNewClusters2D.at(LArClusterHelper::GetClusterHitType(pCluster)).push_back(pCluster);
        newClusters.push_back(pCluster3D);
        newPfoClusters.push_back(newClusters);
    }

    // Put any remaining 2D hits into the nearest new 2D cluster made in the last step
    // If no 2D clusters for a view could be made in the previous step (happens rarely when no 3D hits get made from one of the cluster views),
    // just mop up the remaining 2D hits into the old 2D clusters
    // NOTE Hits added to 2D clusters this way are added as isolated
    for (const auto &[view, caloHits2D] : viewToFreeCaloHits2D)
    {
        if (viewToNewClusters2D.at(view).empty())
        {
            std::string clusterListName;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, GetClusterListName(view, clusterListName));
            const ClusterList *pOldClusters2D {nullptr};
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, clusterListName, pOldClusters2D));
            // This means no hits in this view made any 3D hits, this happens very rarely for events with few hits.
            // In this case, just leave these hits unclustered. Could also consider putting them in their own isolated hit-only cluster?
            if (pOldClusters2D->empty())
                continue;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, MopUpCaloHits(caloHits2D, *pOldClusters2D, true));
            continue;
        }
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, MopUpCaloHits(caloHits2D, viewToNewClusters2D.at(view), true));
    }

    // Create the new Pfos
    // ATTN New pfos will: have no vertex, not be characterised as track/shower, not be part of any particle hierarchy.
    //      If desired, these these things will need to be added in succeeding algs.
    const PfoList *pNewPfoList {nullptr};
    std::string tempPfoListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                             PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pNewPfoList, tempPfoListName));
    for (const ClusterList &pfoClusters : newPfoClusters)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, CreatePfoFromClusters(pfoClusters));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, tempPfoListName, m_pfoListName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDMultiReclusteringAlgorithm::Build2DClustersFrom3D(const Cluster *const pCluster3D,
                                                                   std::map<HitType, CaloHitList> &viewToFreeCaloHits2D,
                                                                   ClusterList &newClusters2D)
{
    std::map<HitType, PandoraContentApi::Cluster::Parameters> viewToParamsCluster2D = {
        {TPC_VIEW_U, {}}, {TPC_VIEW_V, {}}, {TPC_VIEW_W, {}} };

    CaloHitList cluster3DCaloHits3D;
    pCluster3D->GetOrderedCaloHitList().FillCaloHitList(cluster3DCaloHits3D);
    for (const CaloHit *const pCaloHit3D : cluster3DCaloHits3D)
    {
        const CaloHit *const pCaloHit3DParent2D {static_cast<const CaloHit *>(pCaloHit3D->GetParentAddress())};
        const HitType view {pCaloHit3DParent2D->GetHitType()};
        if (std::find(viewToFreeCaloHits2D.at(view).begin(), viewToFreeCaloHits2D.at(view).end(), pCaloHit3DParent2D) ==
            viewToFreeCaloHits2D.at(view).end())
            return STATUS_CODE_FAILURE;
        viewToParamsCluster2D.at(view).m_caloHitList.push_back(pCaloHit3DParent2D);
        viewToFreeCaloHits2D.at(view).remove(pCaloHit3DParent2D);
    }

    CaloHitList cluster3DIsoCaloHits3D {pCluster3D->GetIsolatedCaloHitList()};
    for (const CaloHit *const pCaloHit3D : cluster3DIsoCaloHits3D)
    {
        const CaloHit *const pCaloHit3DParent2D {static_cast<const CaloHit *>(pCaloHit3D->GetParentAddress())};
        const HitType view {pCaloHit3DParent2D->GetHitType()};
        if (std::find(viewToFreeCaloHits2D.at(view).begin(), viewToFreeCaloHits2D.at(view).end(), pCaloHit3DParent2D) ==
            viewToFreeCaloHits2D.at(view).end())
            return STATUS_CODE_FAILURE;
        viewToParamsCluster2D.at(view).m_isolatedCaloHitList.push_back(pCaloHit3DParent2D);
        viewToFreeCaloHits2D.at(view).remove(pCaloHit3DParent2D);
    }

    for (const auto &[view, params] : viewToParamsCluster2D)
    {
        if (params.m_caloHitList.empty() && params.m_isolatedCaloHitList.empty())
            continue;

        std::string tempListName;
        const ClusterList *pTempClusterList {nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                                 PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pTempClusterList, tempListName));
        const Cluster *pNewCluster2D {nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, params, pNewCluster2D));
        std::string cluster2DListName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, GetClusterListName(view, cluster2DListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, tempListName, cluster2DListName));
        newClusters2D.push_back(pNewCluster2D);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDMultiReclusteringAlgorithm::MopUpCaloHits(const CaloHitList &caloHits, const ClusterList &clusters, bool addAsIso)
{
    for (const CaloHit *const pCaloHit : caloHits)
    {
        const Cluster *pNearestCluster {nullptr};
        double minimumDistance {std::numeric_limits<float>::max()};
        for (const Cluster *const pCluster : clusters)
        {
            double dist = LArClusterHelper::GetClosestDistance(pCaloHit->GetPositionVector(), pCluster);
            if (dist < minimumDistance)
            {
                minimumDistance = dist;
                pNearestCluster = pCluster;
            }
        }
        if (!pNearestCluster)
            return STATUS_CODE_FAILURE;

        if (addAsIso)
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddIsolatedToCluster(*this, pNearestCluster, pCaloHit));
        }
        else
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pNearestCluster, pCaloHit));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDMultiReclusteringAlgorithm::CreatePfoFromClusters(const ClusterList &clusters)
{
    PandoraContentApi::ParticleFlowObject::Parameters params;
    params.m_clusterList = clusters;
    params.m_particleId = MU_MINUS; // Arbitrary choice to mark all new pfos as track-like
    params.m_charge = PdgTable::GetParticleCharge(params.m_particleId.Get());
    params.m_mass = PdgTable::GetParticleMass(params.m_particleId.Get());
    params.m_energy = 0.f;
    params.m_momentum = CartesianVector(0.f, 0.f, 0.f);
    const Pfo *pPfo {nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, params, pPfo));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDMultiReclusteringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "PfoListName", m_pfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "Cluster3DListName", m_cluster3DListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterUListName", m_clusterUListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterVListName", m_clusterVListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterWListName", m_clusterWListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle,
        "ClusteringAlgorithms", m_clusteringAlgs));
    AlgorithmTool *pAlgTool {nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle,
        "FomAlgorithmTool", pAlgTool));
    m_pFomAlgTool = dynamic_cast<ThreeDReclusteringFigureOfMeritBaseTool *>(pAlgTool);

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
