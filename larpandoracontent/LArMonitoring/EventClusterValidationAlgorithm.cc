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
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"

#include <numeric>

using namespace pandora;

namespace lar_content
{

EventClusterValidationAlgorithm::CaloHitParents::CaloHitParents() :
    m_pMainMC{nullptr},
    m_pCluster{nullptr}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

EventClusterValidationAlgorithm::ClusterMetrics::ClusterMetrics() :
    m_nRecoHits{std::vector<int>{}},
    m_purities{std::vector<float>{}},
    m_completenesses{std::vector<float>{}},
    m_nHits{0},
    m_nClusters{0},
    m_nMainMCs{0}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

EventClusterValidationAlgorithm::EventClusterValidationAlgorithm() :
    m_eventNumber{0},
    m_caloHitListName{"CaloHitList2D"},
    m_minMCHitsPerView{0}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

EventClusterValidationAlgorithm::~EventClusterValidationAlgorithm()
{
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName, m_fileName, "UPDATE"));

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName + "_meta", "min_mc_hits_per_view", m_minMCHitsPerView));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName + "_meta"));
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName + "_meta", m_fileName, "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventClusterValidationAlgorithm::Run()
{
    m_eventNumber++;

    // Gather hits by view
    const CaloHitList *pFullCaloHitList{nullptr};
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pFullCaloHitList));
    std::map<HitType, CaloHitList> viewToCaloHits;
    for (const CaloHit *const pCaloHit : *pFullCaloHitList)
        viewToCaloHits[pCaloHit->GetHitType()].emplace_back(pCaloHit);

    // Gather clusters by view
    std::map<HitType, ClusterList> viewToClusters;
    for (std::string listName : m_clusterListNames)
    {
        const ClusterList *pClusterList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pClusterList));
        for (const Cluster *const pCluster : *pClusterList)
        {
            const HitType view{LArClusterHelper::GetClusterHitType(pCluster)};
            viewToClusters[view].emplace_back(pCluster);
        }
    }

    const std::vector<ValidationType> valTypes{ValidationType::ALL, ValidationType::SHOWER, ValidationType::TRACK};
    for (const auto &[view, caloHits] : viewToCaloHits)
    {
        const ClusterList &clusters{viewToClusters[view]};

        // For each hit, get the main MCParticle and the cluster it belongs to (if the truth matching is not missing)
        // NOTE Using hitParents map to track which hits are being considered in metric calculations
        std::map<const CaloHit *const, CaloHitParents> hitParents;
        GetHitParents(caloHits, clusters, hitParents);

        // Drop any hits with a main MC particle with insufficient activity
        if (m_minMCHitsPerView > 0)
            ApplyMCParticleMinSumHits(hitParents);

        for (const ValidationType valType : valTypes)
        {
            // Apply the track/shower/all criteria
            std::map<const CaloHit *const, CaloHitParents> hitParentsValid{ApplyPDGCut(hitParents, valType)};

            ClusterMetrics metrics;
            GetMetrics(hitParentsValid, metrics);

            float adjustedRandI{CalcRandIndex(hitParentsValid)};

            std::string branchPrefix;
            if (valType == ValidationType::ALL)
                branchPrefix = "all_";
            else if (valType == ValidationType::SHOWER)
                branchPrefix = "shower_";
            else if (valType == ValidationType::TRACK)
                branchPrefix = "track_";
            SetBranches(metrics, adjustedRandI, branchPrefix);
        }

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "event", m_eventNumber - 1));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "view", static_cast<int>(view)));
        PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventClusterValidationAlgorithm::GetHitParents(
    const CaloHitList &caloHits, const ClusterList &clusters, std::map<const CaloHit *const, CaloHitParents> &hitParents) const
{
    for (const CaloHit *const pCaloHit : caloHits)
    {
        const MCParticle *pMainMC{nullptr};
        const MCParticleWeightMap &weightMap{pCaloHit->GetMCParticleWeightMap()};
        float maxWeight{std::numeric_limits<float>::lowest()};
        for (const auto &[pMC, weight] : weightMap)
        {
            if (weight > maxWeight)
                pMainMC = pMC;
        }
        if (pMainMC)
        {
            hitParents[pCaloHit] = CaloHitParents();
            hitParents[pCaloHit].m_pMainMC = pMainMC;
        }
    }

    for (const Cluster *const pCluster : clusters)
    {
        const CaloHitList &isolatedHits{pCluster->GetIsolatedCaloHitList()};
        CaloHitList clusterCaloHits;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(clusterCaloHits);
        clusterCaloHits.insert(clusterCaloHits.end(), isolatedHits.begin(), isolatedHits.end());
        for (const CaloHit *const pCaloHit : clusterCaloHits)
        {
            // Ignoring the calo hit if truth matching is missing
            if (hitParents.find(pCaloHit) == hitParents.end())
                continue;
            hitParents[pCaloHit].m_pCluster = pCluster;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventClusterValidationAlgorithm::ApplyMCParticleMinSumHits(std::map<const CaloHit *const, CaloHitParents> &hitParents) const
{
    std::map<const MCParticle *const, int> mcSumHits;
    for (const auto &[pCaloHit, parents] : hitParents)
    {
        const MCParticle *const pMainMC{parents.m_pMainMC};
        if (mcSumHits.find(pMainMC) == mcSumHits.end())
            mcSumHits[pMainMC] = 0;
        mcSumHits.at(pMainMC)++;
    }

    for (auto it = hitParents.begin(); it != hitParents.end();)
    {
        if (mcSumHits.at(it->second.m_pMainMC) < m_minMCHitsPerView)
            it = hitParents.erase(it);
        else
            it++;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::map<const CaloHit *const, EventClusterValidationAlgorithm::CaloHitParents> EventClusterValidationAlgorithm::ApplyPDGCut(
    std::map<const CaloHit *const, CaloHitParents> &hitParents, const ValidationType &valType) const
{
    std::map<const CaloHit *const, CaloHitParents> hitParentsValid;
    for (const auto &[pCaloHit, parents] : hitParents)
    {
        const MCParticle *const pMainMC{parents.m_pMainMC};
        if (valType != ValidationType::ALL)
        {
            const int pdg{std::abs(pMainMC->GetParticleId())};
            if ((valType == ValidationType::SHOWER && (pdg != PHOTON && pdg != E_MINUS)) ||
                (valType == ValidationType::TRACK && (pdg == PHOTON || pdg == E_MINUS)))
                continue;
        }
        hitParentsValid[pCaloHit] = parents;
    }

    return hitParentsValid;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventClusterValidationAlgorithm::GetMetrics(const std::map<const CaloHit *const, CaloHitParents> &hitParents, ClusterMetrics &metrics) const
{
    // Adapted from Andy's code for calculating cluster purity and completeness (ClusterValidationAlgorithm)
    std::set<const Cluster *> validClusters;
    std::set<const MCParticle *> validMCParticles;
    for (const auto &[pCaloHit, parents] : hitParents)
    {
        validClusters.insert(parents.m_pCluster);
        validMCParticles.insert(parents.m_pMainMC);
    }
    metrics.m_nHits = hitParents.size();
    metrics.m_nClusters = validClusters.size();
    metrics.m_nMainMCs = validMCParticles.size();

    for (const Cluster *const pCluster : validClusters)
    {
        const CaloHitList &isolatedHits{pCluster->GetIsolatedCaloHitList()};
        CaloHitList clusterCaloHits;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(clusterCaloHits);
        clusterCaloHits.insert(clusterCaloHits.end(), isolatedHits.begin(), isolatedHits.end());

        std::map<const MCParticle *const, int> mainMCParticleHits;
        float totalHits{0};
        for (const CaloHit *const pCaloHit : clusterCaloHits)
        {
            if (hitParents.find(pCaloHit) == hitParents.end())
                continue;

            const MCParticle *const pMC = hitParents.at(pCaloHit).m_pMainMC;
            if (mainMCParticleHits.find(pMC) == mainMCParticleHits.end())
                mainMCParticleHits[pMC] = 0;
            mainMCParticleHits.at(pMC)++;
            totalHits++;
        }
        if (totalHits == 0) // Something went wrong and won't be able to get metrics from this cluster
            continue;

        // Determine which MC particle contributes the most weight across the cluster
        const MCParticle *pMainMC{nullptr};
        int maxHits{0};
        for (const auto &[pMC, nHits] : mainMCParticleHits)
        {
            if (nHits > maxHits)
            {
                pMainMC = pMC;
                maxHits = nHits;
            }
        }
        metrics.m_purities.emplace_back(maxHits / totalHits);
        metrics.m_nRecoHits.emplace_back(totalHits);

        // Calculate cluster completeness
        int nTotalMainMCHits{0};
        for (const auto &[pCaloHit, parents] : hitParents)
        {
            if (pMainMC == parents.m_pMainMC)
                nTotalMainMCHits++;
        }
        metrics.m_completenesses.emplace_back(maxHits / nTotalMainMCHits);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float EventClusterValidationAlgorithm::CalcRandIndex(std::map<const CaloHit *const, CaloHitParents> &hitParents) const
{
    // Fill contingency table
    ContingencyTable cTable;
    for (const auto &[pCaloHit, parents] : hitParents)
    {
        const MCParticle *const pMC = parents.m_pMainMC;
        const Cluster *const pCluster = parents.m_pCluster;

        if (cTable.find(pCluster) == cTable.end() || cTable.at(pCluster).find(pMC) == cTable.at(pCluster).end())
            cTable[pCluster][pMC] = 0;
        cTable.at(pCluster).at(pMC)++;
    }

    return LArMonitoringHelper::CalcRandIndex(cTable);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventClusterValidationAlgorithm::SetBranches(
    [[maybe_unused]] ClusterMetrics &metrics, [[maybe_unused]] float randIndex, [[maybe_unused]] std::string branchPrefix) const
{
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, branchPrefix + "n_hits", metrics.m_nHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, branchPrefix + "n_clusters", metrics.m_nClusters));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, branchPrefix + "n_mainMCs", metrics.m_nMainMCs));
#ifdef MONITORING
    double meanPurity{-1.0}, meanCompleteness{-1.0}, meanNRecoHits{-1.0};
    if (!metrics.m_purities.empty() && !metrics.m_completenesses.empty() && !metrics.m_nRecoHits.empty())
    {
        meanPurity = std::accumulate(metrics.m_purities.begin(), metrics.m_purities.end(), 0.0) / metrics.m_purities.size();
        meanCompleteness = std::accumulate(metrics.m_completenesses.begin(), metrics.m_completenesses.end(), 0.0) / metrics.m_completenesses.size();
        meanNRecoHits = std::accumulate(metrics.m_nRecoHits.begin(), metrics.m_nRecoHits.end(), 0.0) / metrics.m_nRecoHits.size();
    }
#endif
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, branchPrefix + "mean_purity", meanPurity));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, branchPrefix + "mean_completeness", meanCompleteness));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, branchPrefix + "mean_n_reco_hits", meanNRecoHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, branchPrefix + "adjusted_rand_idx", randIndex));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventClusterValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "FileName", m_fileName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "treeName", m_treeName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "ClusterListNames", m_clusterListNames));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinMCHitsPerView", m_minMCHitsPerView));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
