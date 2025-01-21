/**
 *  @file   larpandoracontent/LArMonitoring/EventClusterValidationAlgorithm.cc
 *
 *  @brief  Implementation of the event-level cluster validation.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/EventClusterValidationAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

using namespace pandora;

namespace lar_content
{

EventClusterValidationAlgorithm::ClusterMetrics::ClusterMetrics() :
    m_pMainMCParticles{std::vector<const pandora::MCParticle *> {}},
    m_nContributions{std::vector<int>{}},
    m_nRecoHits{std::vector<int>{}},
    m_purities{std::vector<float>{}},
    m_completenesses{std::vector<float>{}},
    m_fragmentationFractions{std::vector<float>{}}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

EventClusterValidationAlgorithm::EventClusterValidationAlgorithm() :
    m_eventNumber{0},
    m_caloHitListName{"CaloHitList2D"},
    m_minMCHitsPerView{0.f}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

EventClusterValidationAlgorithm::~EventClusterValidationAlgorithm()
{
    // PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName, m_fileName, "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventClusterValidationAlgorithm::Run()
{
    m_eventNumber++;

    const CaloHitList *pFullCaloHitList {};
    PANDORA_THROW_RESULT_IF(
        STATUS_CODE_SUCCESS, !=,
        PandoraContentApi::GetList(*this, m_caloHitListName, pFullCaloHitList)
    );

    std::map<HitType, CaloHitList> viewCaloHits;
    for (const CaloHit *const pCaloHit : *pFullCaloHitList)
        viewCaloHits[pCaloHit->GetHitType()].emplace_back(pCaloHit);

    // Collect the sum of hit weights for each MCParticle on each view
    std::map<HitType, MCSumHitWeightMap> mcIDSumHitWeight;
    GetMCParticleViewSumHits(pFullCaloHitList, mcIDSumHitWeight);

    // Gather clusters from each view list into one container
    std::map<HitType, ClusterList> viewClustersMap;
    for (std::string listName : m_clusterListNames)
    {
        const ClusterList *pClusterList {};
        PANDORA_RETURN_RESULT_IF(
            STATUS_CODE_SUCCESS, !=,
            PandoraContentApi::GetList(*this, listName, pClusterList)
        );
        for (const Cluster *const pCluster : *pClusterList)
        {
            const HitType view {LArClusterHelper::GetClusterHitType(pCluster)};
            viewClustersMap[view].emplace_back(pCluster);
        }
    }

    for (const auto &[view, clusterList] : viewClustersMap) {
        const CaloHitList caloHitList {viewCaloHits[view]};
        ClusterMetrics clusterMetrics;
        GetMetrics(clusterList, caloHitList, clusterMetrics);
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "event", m_eventNumber - 1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "view", static_cast<int>{view}));

        // TODO
        // Write everything out to tree, think I will need to go back an make explicit vectors for each MC particle propery
        // Write out mean purity and completeness too
        // Test this all works and returns purity and completeness in line with the ClusterValidations
        // Apply minMCHit threshold in GetMetrics
        // Calculate Rand index

        // PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "pdg", metrics.m_pMC->GetParticleId()));
        // PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "energy", metrics.m_pMC->GetEnergy()));
        // PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "tier", LArMCParticleHelper::GetHierarchyTier(metrics.m_pMC)));
        // PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "mom_x", mom.GetX()));
        // PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "mom_y", mom.GetY()));
        // PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "mom_z", mom.GetZ()));
        // PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "n_contribs", metrics.m_nContributions));
        // PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "n_mc_hits", metrics.m_nTotalMCHits));
        // PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "weighted_mc_hits", metrics.m_wTotalMCHits));
        // PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "n_reco_hits", metrics.m_nRecoHits));
        // PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "purity", metrics.m_purity));
        // PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "completeness", metrics.m_completeness));
        // PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "fragmentation", metrics.m_fragmentationFraction));
        // PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName));
    }


    // ClusterMetricsMap metricsMap;
    // this->GetMetrics(viewClustersMap, metricsMap);

    // if (m_writeFile)
    // {
    //     for (const auto &[pCluster, metrics] : metricsMap)
    //     {
    //         int view{LArClusterHelper::GetClusterHitType(pCluster)};
    //         const CartesianVector &mom{metrics.m_pMC->GetMomentum()};
    //         PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "event", m_eventNumber - 1));
    //         PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "view", view));
    //         PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "pdg", metrics.m_pMC->GetParticleId()));
    //         PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "energy", metrics.m_pMC->GetEnergy()));
    //         PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "tier", LArMCParticleHelper::GetHierarchyTier(metrics.m_pMC)));
    //         PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "mom_x", mom.GetX()));
    //         PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "mom_y", mom.GetY()));
    //         PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "mom_z", mom.GetZ()));
    //         PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "n_contribs", metrics.m_nContributions));
    //         PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "n_mc_hits", metrics.m_nTotalMCHits));
    //         PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "weighted_mc_hits", metrics.m_wTotalMCHits));
    //         PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "n_reco_hits", metrics.m_nRecoHits));
    //         PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "purity", metrics.m_purity));
    //         PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "completeness", metrics.m_completeness));
    //         PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "fragmentation", metrics.m_fragmentationFraction));
    //         PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName));
    //     }
    // }

    // if (m_visualize)
    // {
    //     PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    // }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventClusterValidationAlgorithm::GetMCParticleViewSumHits(
    const pandora::CaloHitList *const pCaloHitList, std::map<HitType, MCSumHitWeightMap> &mcIDSumHitWeight
) const
{
    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        const HitType view {pCaloHit->GetHitType()};
        if (mcIDSumHitWeight.find(view) == mcIDSumHitWeight.end())
            mcIDSumHitWeight[view] = MCSumHitWeightMap{};

        const MCParticleWeightMap &weightMap {pCaloHit->GetMCParticleWeightMap()};
        for (const auto &[pMC, weight] : weightMap)
        {
            if (mcIDSumHitWeight[view].find(pMC->GetUid()) == mcIDSumHitWeight[view].end())
                mcIDSumHitWeight[view][pMC->GetUid()] = 0.0f;
            mcIDSumHitWeight[view][pMC->GetUid()] += weight;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventClusterValidationAlgorithm::GetMetrics(
    const ClusterList &viewClusters, const CaloHitList &viewCaloHits, ClusterMetrics &metrics
) const
{
    // Adapted from Andy's code for calculating cluster purity and completeness (ClusterValidationAlgorithm)
    for (const Cluster *const pCluster : viewClusters)
    {
        const CaloHitList &isolatedHits {pCluster->GetIsolatedCaloHitList()};
        CaloHitList caloHits;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHits);
        caloHits.insert(caloHits.end(), isolatedHits.begin(), isolatedHits.end());

        // Aggregate MC particle weights from each calo hit in the cluster
        MCParticleWeightMap clusterMCWeightMap;
        float totalClusterWeight {0.f};
        for (const CaloHit *const pCaloHit : caloHits)
        {
            const MCParticleWeightMap &weightMap {pCaloHit->GetMCParticleWeightMap()};
            for (const auto &[pMC, weight] : weightMap)
            {
                if (clusterMCWeightMap.find(pMC) == clusterMCWeightMap.end())
                    clusterMCWeightMap[pMC] = 0.f;
                clusterMCWeightMap[pMC] += weight;
                totalClusterWeight += weight;
            }
        }
        if (totalClusterWeight > 0)
        {
            for (const auto &[pMC, weight] : clusterMCWeightMap)
                clusterMCWeightMap[pMC] /= totalClusterWeight;
        }

        // Determine which MC particle contributes the most weight across the cluster
        const MCParticle *pMainMC {nullptr};
        float maxWeight {std::numeric_limits<float>::lowest()};
        for (const auto &[pMC, weight] : clusterMCWeightMap)
        {
            if (weight > maxWeight)
            {
                pMainMC = pMC;
                maxWeight = weight;
            }
        }

        // Calculate the cluster purity
        float wMatchedHits {0};
        int nTotalRecoHits {0};
        for (const CaloHit *const pCaloHit : caloHits)
        {
            const MCParticleWeightMap &weightMap {pCaloHit->GetMCParticleWeightMap()};
            if (weightMap.find(pMainMC) != weightMap.end())
                wMatchedHits += weightMap.at(pMainMC);
            if (!weightMap.empty())
                nTotalRecoHits++;
        }
        if (nTotalRecoHits > 0)
        {
            metrics.m_purities.push_back(wMatchedHits / nTotalRecoHits);
            metrics.m_pMainMCParticles.push_back(pMainMC);
            metrics.m_nRecoHits.push_back(nTotalRecoHits);
        }

        // Calculate Cluster completeness
        float wTotalMCHits {0};
        for (const CaloHit *const pCaloHit : viewCaloHits)
        {
            const MCParticleWeightMap &weightMap {pCaloHit->GetMCParticleWeightMap()};
            if (weightMap.find(pMainMC) != weightMap.end())
                wTotalMCHits += weightMap.at(pMainMC);
        }

        if (wTotalMCHits > 0)
            metrics.m_completenesses.push_back(wMatchedHits / wTotalMCHits);
        float accountedWeight {0.f};
        int nContributions {0};
        for (const auto &[pMC, weight] : clusterMCWeightMap)
        {
            if (weight > 0.2f)
            {
                ++nContributions;
                accountedWeight += weight;
            }
        }
        metrics.m_nContributions.push_back(nContributions);
        metrics.m_fragmentationFractions.push_back(1.f - accountedWeight);
    }
}
// {
//     const CaloHitList *pFullCaloHitList{nullptr};
//     PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pFullCaloHitList));
//     for (const auto &[view, clusterList] : viewClustersMap)
//     {
//         CaloHitList viewCaloHitList;
//         for (const CaloHit *const pCaloHit : *pFullCaloHitList)
//             if (pCaloHit->GetHitType() == view)
//                 viewCaloHitList.emplace_back(pCaloHit);
//         for (const Cluster *const pCluster : clusterList)
//         {
//             MCParticleWeightMap clusterMCWeightMap;
//             const CaloHitList &isolatedHits{pCluster->GetIsolatedCaloHitList()};
//             CaloHitList caloHits;
//             pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHits);
//             caloHits.insert(caloHits.end(), isolatedHits.begin(), isolatedHits.end());
//             float totalClusterWeight{0.f};
//             for (const CaloHit *const pCaloHit : caloHits)
//             {
//                 const MCParticleWeightMap &weightMap{pCaloHit->GetMCParticleWeightMap()};
//                 for (const auto &[pMC, weight] : weightMap)
//                 {
//                     if (clusterMCWeightMap.find(pMC) == clusterMCWeightMap.end())
//                         clusterMCWeightMap[pMC] = 0.f;
//                     clusterMCWeightMap[pMC] += weight;
//                     totalClusterWeight += weight;
//                 }
//             }
//             if (totalClusterWeight > 0)
//                 for (const auto &[pMC, weight] : clusterMCWeightMap)
//                     clusterMCWeightMap[pMC] /= totalClusterWeight;
//             // Determine which MC particle contributes the most weight across the cluster
//             const MCParticle *pMainMC{metricsMap.find(pCluster) != metricsMap.end() ? metricsMap[pCluster].m_pMC : nullptr};
//             if (!pMainMC)
//             {
//                 float maxWeight{std::numeric_limits<float>::lowest()};
//                 for (const auto &[pMC, weight] : clusterMCWeightMap)
//                 {
//                     if (weight > maxWeight)
//                     {
//                         pMainMC = pMC;
//                         maxWeight = weight;
//                     }
//                 }
//             }
//             if (pMainMC)
//             {
//                 float wMatchedHits{0};
//                 int nTotalRecoHits{0};
//                 for (const CaloHit *const pCaloHit : caloHits)
//                 {
//                     const MCParticleWeightMap &weightMap{pCaloHit->GetMCParticleWeightMap()};
//                     if (weightMap.find(pMainMC) != weightMap.end())
//                     {
//                         wMatchedHits += weightMap.at(pMainMC);
//                     }
//                     if (!weightMap.empty())
//                         nTotalRecoHits++;
//                 }
//                 if (nTotalRecoHits > 0)
//                 {
//                     metricsMap[pCluster].m_purity = wMatchedHits / nTotalRecoHits;
//                     metricsMap[pCluster].m_pMC = pMainMC;
//                     metricsMap[pCluster].m_nRecoHits = nTotalRecoHits;
//                 }
//                 float wTotalMCHits{0};
//                 int nTotalMCHits{0};
//                 for (const CaloHit *const pCaloHit : viewCaloHitList)
//                 {
//                     const MCParticleWeightMap &weightMap{pCaloHit->GetMCParticleWeightMap()};
//                     if (weightMap.find(pMainMC) != weightMap.end())
//                     {
//                         wTotalMCHits += weightMap.at(pMainMC);
//                         nTotalMCHits++;
//                     }
//                 }

//                 if (wTotalMCHits > 0)
//                 {
//                     metricsMap[pCluster].m_completeness = wMatchedHits / wTotalMCHits;
//                     metricsMap[pCluster].m_nTotalMCHits = nTotalMCHits;
//                     metricsMap[pCluster].m_wTotalMCHits = wTotalMCHits;
//                 }
//                 float accountedWeight{0.f};
//                 for (const auto &[pMC, weight] : clusterMCWeightMap)
//                 {
//                     if (weight > 0.2f)
//                     {
//                         ++metricsMap[pCluster].m_nContributions;
//                         accountedWeight += weight;
//                     }
//                 }
//                 metricsMap[pCluster].m_fragmentationFraction = 1.f - accountedWeight;
//             }
//         }
//     }
// }

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventClusterValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(
        STATUS_CODE_SUCCESS, !=,
        XmlHelper::ReadValue(xmlHandle, "FileName", m_fileName)
    );
    PANDORA_RETURN_RESULT_IF(
        STATUS_CODE_SUCCESS, !=,
        XmlHelper::ReadValue(xmlHandle, "TreeName", m_treeName)
    );
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName)
    );
    PANDORA_RETURN_RESULT_IF(
        STATUS_CODE_SUCCESS, !=,
        XmlHelper::ReadVectorOfValues(xmlHandle, "ClusterListNames", m_clusterListNames)
    );
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinMCHitsPerView", m_minMCHitsPerView)
    );

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

