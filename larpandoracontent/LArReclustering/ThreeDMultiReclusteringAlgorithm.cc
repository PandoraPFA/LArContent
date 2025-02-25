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
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDMultiReclusteringAlgorithm::~ThreeDMultiReclusteringAlgorithm()
{
}

StatusCode ThreeDMultiReclusteringAlgorithm::DebugCurrentLists()
{
    const PfoList *pPfos {nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_pfoListName, pPfos));
    const ClusterList *pClusters3D {nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_cluster3DListName, pClusters3D));
    const ClusterList *pClustersU {nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_clusterUListName, pClustersU));
    const ClusterList *pClustersV {nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_clusterVListName, pClustersV));
    const ClusterList *pClustersW {nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_clusterWListName, pClustersW));

    std::cout << "====\n";
    std::cout << "--- " << pPfos->size() << " pfos ---\n";
    int totalCaloHits3D {0}, totalCaloHitsU {0}, totalCaloHitsV {0}, totalCaloHitsW {0};
    int totalIsoCaloHits3D {0}, totalIsoCaloHitsU {0}, totalIsoCaloHitsV {0}, totalIsoCaloHitsW {0};
    for (const Pfo *const pPfo : *pPfos)
    {
        ClusterList clusters3D;
        LArPfoHelper::GetClusters(pPfo, TPC_3D, clusters3D);
        ClusterList clustersU;
        LArPfoHelper::GetClusters(pPfo, TPC_VIEW_U, clustersU);
        ClusterList clustersV;
        LArPfoHelper::GetClusters(pPfo, TPC_VIEW_V, clustersV);
        ClusterList clustersW;
        LArPfoHelper::GetClusters(pPfo, TPC_VIEW_W, clustersW);
        std::cout << pPfo << ": " << clusters3D.size() << " - ";
        if (clusters3D.size() != 0)
        {
            std::cout << clusters3D.front()->GetNCaloHits() << " / " << clusters3D.front()->GetNIsolatedCaloHits() << " (";
            totalCaloHits3D += clusters3D.front()->GetNCaloHits();
            totalIsoCaloHits3D += clusters3D.front()->GetNIsolatedCaloHits();
        }
        else
            std::cout << "0 / 0 (";
        if (clustersU.size() != 0)
        {
            std::cout << clustersU.size() << " - " << clustersU.front()->GetNCaloHits() << " / " << clustersU.front()->GetNIsolatedCaloHits() << ", ";
            totalCaloHitsU += clustersU.front()->GetNCaloHits();
            totalIsoCaloHitsU += clustersU.front()->GetNIsolatedCaloHits();
        }
        else
            std::cout << "0 / 0, ";
        if (clustersV.size() != 0)
        {
            std::cout << clustersV.size() << " - " << clustersV.front()->GetNCaloHits() << " / " << clustersV.front()->GetNIsolatedCaloHits() << ", ";
            totalCaloHitsV += clustersV.front()->GetNCaloHits();
            totalIsoCaloHitsV += clustersV.front()->GetNIsolatedCaloHits();
        }
        else
            std::cout << "0 / 0, ";
        if (clustersW.size() != 0)
        {
            std::cout << clustersW.size() << " - " << clustersW.front()->GetNCaloHits() << " / " << clustersW.front()->GetNIsolatedCaloHits() << ")\n";
            totalCaloHitsW += clustersW.front()->GetNCaloHits();
            totalIsoCaloHitsW += clustersW.front()->GetNIsolatedCaloHits();
        }
        else
            std::cout << "0 / 0)\n";
    }
    std::cout << " ---> total: " << totalCaloHits3D << " / " << totalIsoCaloHits3D << " ("
                                 << totalCaloHitsU << " / " << totalIsoCaloHitsU << ", "
                                 << totalCaloHitsV << " / " << totalIsoCaloHitsV << ", "
                                 << totalCaloHitsW << " / " << totalIsoCaloHitsW << ")\n";
    std::cout << "--- " << pClusters3D->size() << " 3D clusters ---\n";
    int totalCaloHits {0}, totalIsoCaloHits {0};
    for (const Cluster *const pCluster : *pClusters3D)
    {
        std::cout << pCluster << ": " << pCluster->GetNCaloHits() << " / " << pCluster->GetNIsolatedCaloHits() << "\n";
        totalCaloHits += pCluster->GetNCaloHits();
        totalIsoCaloHits += pCluster->GetNIsolatedCaloHits();
    }
    std::cout << " ---> total: " << totalCaloHits << " / " << totalIsoCaloHits << "\n";
    totalCaloHits = 0; totalIsoCaloHits = 0;
    std::cout << "--- " << pClustersU->size() << " 2D U clusters ---\n";
    for (const Cluster *const pCluster : *pClustersU)
    {
        std::cout << pCluster << ": " << pCluster->GetNCaloHits() << " / " << pCluster->GetNIsolatedCaloHits() << "\n";
        totalCaloHits += pCluster->GetNCaloHits();
        totalIsoCaloHits += pCluster->GetNIsolatedCaloHits();
    }
    std::cout << " ---> total: " << totalCaloHits << " / " << totalIsoCaloHits << "\n";
    totalCaloHits = 0; totalIsoCaloHits = 0;
    std::cout << "--- " << pClustersV->size() << " 2D V clusters ---\n";
    for (const Cluster *const pCluster : *pClustersV)
    {
        std::cout << pCluster << ": " << pCluster->GetNCaloHits() << " / " << pCluster->GetNIsolatedCaloHits() << "\n";
        totalCaloHits += pCluster->GetNCaloHits();
        totalIsoCaloHits += pCluster->GetNIsolatedCaloHits();
    }
    std::cout << " ---> total: " << totalCaloHits << " / " << totalIsoCaloHits << "\n";
    totalCaloHits = 0; totalIsoCaloHits = 0;
    std::cout << "--- " << pClustersW->size() << " 2D W clusters ---\n";
    for (const Cluster *const pCluster : *pClustersW)
    {
        std::cout << pCluster << ": " << pCluster->GetNCaloHits() << " / " << pCluster->GetNIsolatedCaloHits() << "\n";
        totalCaloHits += pCluster->GetNCaloHits();
        totalIsoCaloHits += pCluster->GetNIsolatedCaloHits();
    }
    std::cout << " ---> total: " << totalCaloHits << " / " << totalIsoCaloHits << "\n";
    std::cout << "====\n";

    return STATUS_CODE_SUCCESS;
}

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

    // DebugCurrentLists();

    // Remove clusters from the pfos, keeping track of their original pfo
    std::map<HitType, ClusterList> viewToFreeClusters = { {TPC_3D, {}}, {TPC_VIEW_U, {}}, {TPC_VIEW_V, {}}, {TPC_VIEW_W, {}} };
    std::map<const Pfo *const, ClusterList> pfoToFreeClusters;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, FreeClustersFromPfos(pfosToRecluster, viewToFreeClusters, pfoToFreeClusters));

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

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        PandoraContentApi::TemporarilyReplaceCurrentList<Cluster>(*this, bestReclusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndReclustering(*this, bestReclusterListName));

    if (bestReclusterListName == originalClusterListName) // Return the clusters to the original pfos as if nothing happened...
    {
        // std::cout << "No better reclustering found\n";
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, AddClustersToPfos(pfoToFreeClusters));
        // DebugCurrentLists();
        return STATUS_CODE_SUCCESS;
    }
    else // Delete the original pfos in preparation for making new ones
    {
        for (const Pfo *const pPfo : pfosToRecluster)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pPfo, m_pfoListName));
    }
    // std::cout << "Reclustering found: " << pfosToRecluster.size() << " -> " << bestReclusterList.size() << "\n";

    // The 3D clusters are reclustered, now need to manually recluster the 2D clusters to be consistent with the 3D reclustering
    // First, remove the original 2D clusters
    std::map<HitType, CaloHitList> viewToFreeCaloHits2D = { {TPC_VIEW_U, {}}, {TPC_VIEW_V, {}}, {TPC_VIEW_W, {}} };
    std::map<HitType, CaloHitList> viewToFreeIsoCaloHits2D = { {TPC_VIEW_U, {}}, {TPC_VIEW_V, {}}, {TPC_VIEW_W, {}} };
    for (const auto &[view, freeClusters] : viewToFreeClusters)
    {
        if (view == TPC_3D)
            continue;

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
            FreeCaloHitsFromClusters(freeClusters, view, viewToFreeCaloHits2D.at(view), viewToFreeIsoCaloHits2D.at(view)));
    }

    // Now make new 2D clusters following the new 3D clusters, and make new corresponding pfos
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, BuildNewPfos(bestReclusterList, viewToFreeCaloHits2D, viewToFreeIsoCaloHits2D));

    // DebugCurrentLists();

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
                                                                      CaloHitList &freeCaloHits,
                                                                      CaloHitList &freeIsoCaloHits)
{
    std::string clusterListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, GetClusterListName(view, clusterListName));
    for (const Cluster *const pCluster: clusters)
    {
        CaloHitList caloHits;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHits);
        freeCaloHits.insert(freeCaloHits.end(), caloHits.begin(), caloHits.end());
        freeIsoCaloHits.insert(freeIsoCaloHits.end(), pCluster->GetIsolatedCaloHitList().begin(), pCluster->GetIsolatedCaloHitList().end());
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pCluster, clusterListName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDMultiReclusteringAlgorithm::BuildNewPfos(const ClusterList &clusters,
                                                          std::map<HitType, CaloHitList> &viewToFreeCaloHits2D,
                                                          std::map<HitType, CaloHitList> &viewToFreeIsoCaloHits2D)
{
    // Cluster any 2D hits directly associated with 3D hits in the same way as the 3D clusters
    std::map<HitType, ClusterList> viewToNewClusters2D = { {TPC_VIEW_U, {}}, {TPC_VIEW_V, {}}, {TPC_VIEW_W, {}} };
    std::vector<ClusterList> newPfoClusters;
    for (const Cluster *const pCluster3D : clusters)
    {
        ClusterList newClusters;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
            Build2DClustersFrom3D(pCluster3D, viewToFreeCaloHits2D, viewToFreeIsoCaloHits2D, newClusters));
        for (const Cluster *const pCluster : newClusters)
            viewToNewClusters2D.at(LArClusterHelper::GetClusterHitType(pCluster)).push_back(pCluster);
        newClusters.push_back(pCluster3D);
        newPfoClusters.push_back(newClusters);
    }

    // Put any remaining 2D hits into the nearest new 2D cluster made in the last step
    for (const auto &[view, caloHits2D] : viewToFreeCaloHits2D)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, MopUpCaloHits(caloHits2D, viewToNewClusters2D.at(view), false));
    for (const auto &[view, caloHits2D] : viewToFreeIsoCaloHits2D)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, MopUpCaloHits(caloHits2D, viewToNewClusters2D.at(view), true));

    // Create the new Pfos
    // NOTE New pfo characterisation will be required
    const PfoList *pNewPfoList {nullptr};
    std::string tempPfoListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pNewPfoList, tempPfoListName));
    for (const ClusterList &pfoClusters : newPfoClusters)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, CreatePfoFromClusters(pfoClusters));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        PandoraContentApi::SaveList<Pfo>(*this, tempPfoListName, m_pfoListName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDMultiReclusteringAlgorithm::Build2DClustersFrom3D(const Cluster *const pCluster3D,
                                                                   std::map<HitType, CaloHitList> &viewToFreeCaloHits2D,
                                                                   std::map<HitType, CaloHitList> &viewToFreeIsoCaloHits2D,
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
        if (std::find(viewToFreeCaloHits2D.at(view).begin(), viewToFreeCaloHits2D.at(view).end(), pCaloHit3DParent2D) !=
            viewToFreeCaloHits2D.at(view).end())
            viewToFreeCaloHits2D.at(view).remove(pCaloHit3DParent2D);
        else if (std::find(viewToFreeIsoCaloHits2D.at(view).begin(), viewToFreeIsoCaloHits2D.at(view).end(), pCaloHit3DParent2D) !=
            viewToFreeIsoCaloHits2D.at(view).end())
            viewToFreeIsoCaloHits2D.at(view).remove(pCaloHit3DParent2D);
        else
            return STATUS_CODE_FAILURE;
        viewToParamsCluster2D.at(view).m_caloHitList.push_back(pCaloHit3DParent2D);
    }

    CaloHitList cluster3DIsoCaloHits3D {pCluster3D->GetIsolatedCaloHitList()};
    for (const CaloHit *const pCaloHit3D : cluster3DIsoCaloHits3D)
    {
        const CaloHit *const pCaloHit3DParent2D {static_cast<const CaloHit *>(pCaloHit3D->GetParentAddress())};
        const HitType view {pCaloHit3DParent2D->GetHitType()};
        if (std::find(viewToFreeCaloHits2D.at(view).begin(), viewToFreeCaloHits2D.at(view).end(), pCaloHit3DParent2D) !=
            viewToFreeCaloHits2D.at(view).end())
            viewToFreeCaloHits2D.at(view).remove(pCaloHit3DParent2D);
        else if (std::find(viewToFreeIsoCaloHits2D.at(view).begin(), viewToFreeIsoCaloHits2D.at(view).end(), pCaloHit3DParent2D) !=
            viewToFreeIsoCaloHits2D.at(view).end())
            viewToFreeIsoCaloHits2D.at(view).remove(pCaloHit3DParent2D);
        else
            return STATUS_CODE_FAILURE;
        viewToParamsCluster2D.at(view).m_isolatedCaloHitList.push_back(pCaloHit3DParent2D);
    }

    for (const auto &[view, params] : viewToParamsCluster2D)
    {
        if (params.m_caloHitList.empty())
            continue;

        std::string cluster2DListName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, GetClusterListName(view, cluster2DListName));

        std::string tempListName;
        const ClusterList *pTempClusterList {nullptr};
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
            PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pTempClusterList, tempListName));
        const Cluster *pNewCluster2D {nullptr};
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, params, pNewCluster2D));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
            PandoraContentApi::SaveList<Cluster>(*this, tempListName, cluster2DListName));
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
        if (addAsIso)
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddIsolatedToCluster(*this, pNearestCluster, pCaloHit));
        }
        else
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pNearestCluster, pCaloHit));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDMultiReclusteringAlgorithm::CreatePfoFromClusters(const ClusterList &clusters)
{
    PandoraContentApi::ParticleFlowObject::Parameters params;
    params.m_clusterList = clusters;
    params.m_particleId = MU_MINUS;
    params.m_charge = PdgTable::GetParticleCharge(params.m_particleId.Get());
    params.m_mass = PdgTable::GetParticleMass(params.m_particleId.Get());
    params.m_energy = 0.f;
    params.m_momentum = CartesianVector(0.f, 0.f, 0.f);
    const Pfo *pPfo {nullptr};
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, params, pPfo));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDMultiReclusteringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "Cluster3DListName", m_cluster3DListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ClusterUListName", m_clusterUListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ClusterVListName", m_clusterVListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ClusterWListName", m_clusterWListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        XmlHelper::ProcessAlgorithmList(*this, xmlHandle, "ClusteringAlgorithms", m_clusteringAlgs));
    AlgorithmTool *pAlgTool {nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "FomAlgorithmTool", pAlgTool));
    m_pFomAlgTool = dynamic_cast<ThreeDReclusteringFigureOfMeritBaseTool *>(pAlgTool);

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
