/**
 *  @file   larpandoracontent/LArMonitoring/EventClusterValidationAlgorithm.cc
 *
 *  @brief  Implementation of the event-level cluster validation.
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/EventClusterValidationAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include <numeric>
#include <algorithm>

using namespace pandora;

namespace lar_content
{

EventClusterValidationAlgorithm::CaloHitParents::CaloHitParents() :
    m_pMainMC{nullptr},
    m_pCluster{nullptr},
    m_pClusterMainMC{nullptr}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

EventClusterValidationAlgorithm::ClusterMetrics::ClusterMetrics() :
    m_purity{0.},
    m_showerPurity{0.},
    m_trackPurity{0.},
    m_completeness{0.},
    m_showerCompleteness{0.},
    m_trackCompleteness{0.},
    m_ari{0.},
    m_showerAri{0.},
    m_trackAri{0.},
    m_nHits{0},
    m_nHitsNullCluster{0},
    m_nShowerTrueHits{0},
    m_nTrackTrueHits{0},
    m_nRecoClusters{0},
    m_nShowerRecoClusters{0},
    m_nTrackRecoClusters{0},
    m_nTrueClusters{0},
    m_nShowerTrueClusters{0},
    m_nTrackTrueClusters{0},
    m_nAriRecoClusters{0},
    m_nShowerAriRecoClusters{0},
    m_nTrackAriRecoClusters{0}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

EventClusterValidationAlgorithm::MatchedParticleMetrics::MatchedParticleMetrics() :
    m_pdg{std::vector<int>{}},
    m_causesShower{std::vector<int>{}},
    m_isPrimary{std::vector<int>{}},
    m_trueEnergy{std::vector<float>{}},
    m_nTrueHits{std::vector<int>{}},
    m_nMatchedCorrectHits{std::vector<int>{}},
    m_nMatchedTotalHits{std::vector<int>{}}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

EventClusterValidationAlgorithm::EventClusterValidationAlgorithm() :
    m_eventNumber{0},
    m_deltaRayLengthThresholdSquared{std::map<HitType, float>{}},
    m_deltaRayParentWeightThreshold{0.f},
    m_caloHitListNames{ { "CaloHitList2D" } },
    m_minMCHitsPerView{0},
    m_onlyRandIndex{false},
    m_foldShowers{false},
    m_handleDeltaRays{false},
    m_mergeShowerClustersForRandIndex{false},
    m_visualize{false},
    m_matchedParticleMetrics{false},
    m_dropNullClusterHits{false},
    m_hitWeightedPurityCompleteness{false},
    m_maximalMatching{false}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

EventClusterValidationAlgorithm::~EventClusterValidationAlgorithm()
{
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName, m_fileName, "UPDATE"));
    if (m_matchedParticleMetrics)
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName + "_matching", m_fileName, "UPDATE"));
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName + "_meta", "min_mc_hits_per_view", m_minMCHitsPerView));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName + "_meta", "fold_showers", m_foldShowers ? 1 : 0));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName + "_meta", "handle_delta_rays", m_handleDeltaRays ? 1 : 0));
    PANDORA_MONITORING_API(
        SetTreeVariable(
            this->GetPandora(), m_treeName + "_meta", "merge_shower_clusters_for_rand_index", m_mergeShowerClustersForRandIndex ? 1 : 0));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName + "_meta", "maximal_matching", m_maximalMatching ? 1 : 0));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName + "_meta"));
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName + "_meta", m_fileName, "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventClusterValidationAlgorithm::Run()
{
    m_eventNumber++;

    // Gather hits by view
    std::map<HitType, CaloHitList> viewToCaloHits;
    for (std::string listName : m_caloHitListNames)
    {
      const CaloHitList *pCaloHitList{nullptr};
      PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pCaloHitList));
      for (const CaloHit *const pCaloHit : *pCaloHitList)
      {
          viewToCaloHits[pCaloHit->GetHitType()].emplace_back(pCaloHit);
      }
    }

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

    for (const auto &[view, caloHits] : viewToCaloHits)
    {
        if (viewToClusters.find(view) == viewToClusters.end()) // Seen it happen for a 1 hit event
        {
            continue;
        }
        const ClusterList &clusters{viewToClusters.at(view)};

        ClusterMetrics clusterMetrics;
        MatchedParticleMetrics matchedParticleMetrics;

        // For each hit, get the main MCParticle and the cluster it belongs to (if the truth matching is not missing)
        // NOTE Using hitParents map to track which hits are being considered in metric calculations
        std::map<const CaloHit *const, CaloHitParents> hitParents;
        this->GetHitParents(caloHits, clusters, hitParents);
        clusterMetrics.m_nHits = hitParents.size();

        // Drop any hits with a main MC particle with insufficient activity
        this->ApplyMCParticleMinSumHits(hitParents);

        clusterMetrics.m_nHitsNullCluster = this->HandleNullClusterHits(hitParents);

        if (hitParents.empty())
        {
            continue;
        }

        this->GetClusterMainMC(hitParents);

        if (m_visualize)
        {
            std::cout << "-- view " << (view == TPC_VIEW_U ? "U" : (view == TPC_VIEW_V ? "V" : "W")) << "\n";
            std::cout << "--- Drawing target true clustering...\n";
            this->VisualizeTargetClusters(hitParents);
        }

        this->CalcRandIndex(hitParents, clusterMetrics);

        if (m_onlyRandIndex)
        {
            this->SetBranches(clusterMetrics, matchedParticleMetrics, view);
        }
        else
        {
            this->GetClusterMetrics(hitParents, clusterMetrics);
            if (m_matchedParticleMetrics)
            {
                this->GetMatchedParticleMetrics(hitParents, matchedParticleMetrics);
            }
            this->SetBranches(clusterMetrics, matchedParticleMetrics, view);
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventClusterValidationAlgorithm::GetHitParents(
        const CaloHitList &caloHits, const ClusterList &clusters, std::map<const CaloHit *const, CaloHitParents> &hitParents) const
{
    std::map<const MCParticle *const, const MCParticle *const> mcFoldTo;
    for (const CaloHit *const pCaloHit : caloHits)
    {
        MCParticleWeightMap weightMap{pCaloHit->GetMCParticleWeightMap()};

        if(m_foldShowers)
        {
            MCParticleWeightMap foldedWeightMap;
            for (const auto &[pMC, weight] : weightMap)
            {
                const MCParticle *pFoldedMC{nullptr};
                if (mcFoldTo.find(pMC) != mcFoldTo.end())
                {
                    pFoldedMC = mcFoldTo.at(pMC);
                }
                else
                {
                    pFoldedMC = this->FoldMCTo(pMC);
                    mcFoldTo.insert({pMC, pFoldedMC});
                }
                foldedWeightMap[pFoldedMC] += weight;
            }
            weightMap = foldedWeightMap;
        }

        const MCParticle *pMainMC{nullptr};
        float maxWeight{0.f};
        for (const auto &[pMC, weight] : weightMap)
        {
            if (weight > maxWeight)
            {
                pMainMC = pMC;
                maxWeight = weight;
            }
            if (weight == maxWeight) // tie-breaker (very unlikely)
            {
                if (LArMCParticleHelper::SortByMomentum(pMC, pMainMC))
                {
                    pMainMC = pMC;
                }
            }
        }
        if (pMainMC)
        {
            hitParents[pCaloHit] = CaloHitParents();
            hitParents[pCaloHit].m_pMainMC = pMainMC;
        }
    }

    if (m_handleDeltaRays)
    {
        for (auto &[pCaloHit, parents] : hitParents)
        {
            parents.m_pMainMC = this->FoldPotentialDeltaRayTo(pCaloHit, parents.m_pMainMC);
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
            if (hitParents.find(pCaloHit) == hitParents.end())
            {
                continue;
            }
            hitParents[pCaloHit].m_pCluster = pCluster;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const MCParticle* EventClusterValidationAlgorithm::FoldMCTo(const MCParticle *const pMC) const
{
    if (!this->IsEM(pMC))
    {
        return pMC;
    }

    bool hasAncestorElectron{false};
    const MCParticle *pCurrentMC{pMC};
    const MCParticle *pLeadingMC{pMC};
    while (!pCurrentMC->IsRootParticle())
    {
        const MCParticle *const pParentMC{*(pCurrentMC->GetParentList().begin())};
        const int parentPdg{std::abs(pParentMC->GetParticleId())};
        if (parentPdg == PHOTON || parentPdg == E_MINUS)
        {
            if (parentPdg == E_MINUS)
            {
                hasAncestorElectron = true;
            }
            pCurrentMC = pParentMC;
            continue;
        }
        pLeadingMC = pCurrentMC;
        break;
    }

    // Don't fold "showers" that consist of only compton scatters
    // Trying to prevents distant diffuse hits disconnected from a "real" bremm + pair production shower being clustered together
    if (hasAncestorElectron || std::abs(pLeadingMC->GetParticleId()) == E_MINUS)
    {
        return pLeadingMC;
    }
    else if (this->CausesShower(pLeadingMC, 0))
    {
        return pLeadingMC;
    }
    else
    {
        return pMC;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool EventClusterValidationAlgorithm::CausesShower(const MCParticle *const pMC, int nDescendentElectrons) const
{
    if (nDescendentElectrons > 1)
    {
        return true;
    }

    if (std::abs(pMC->GetParticleId()) == E_MINUS)
    {
        nDescendentElectrons++; // Including the parent particle, ie. the first in the recursion, as a descendent
    }
    for (const MCParticle *pChildMC : pMC->GetDaughterList())
    {
        if (this->CausesShower(pChildMC, nDescendentElectrons))
        {
            return true;
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool EventClusterValidationAlgorithm::IsEM(const pandora::MCParticle *const pMC) const
{
    const int pdg{std::abs(pMC->GetParticleId())};
    return (pdg == E_MINUS || pdg == PHOTON);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const MCParticle* EventClusterValidationAlgorithm::FoldPotentialDeltaRayTo(const CaloHit *const pCaloHit, const MCParticle *const pMC) const
{
    // Not an electron -> not a delta ray -> do nothing
    if (pMC->IsRootParticle() || pMC->GetParticleId() != E_MINUS)
    {
        return pMC;
    }

    // Did not come from a track-like particle -> not a delta ray -> do nothing
    const MCParticle *const pParentMC{*(pMC->GetParentList().begin())};
    const int parentPdg{std::abs(pParentMC->GetParticleId())};
    if (parentPdg == PHOTON || parentPdg == E_MINUS || PdgTable::GetParticleCharge(parentPdg) == 0)
    {
        return pMC;
    }

    // Delta ray that does not start a shower and is short -> fold into parent particle
    if (!this->CausesShower(pMC, 0) &&
        (pMC->GetVertex() - pMC->GetEndpoint()).GetMagnitudeSquared() < m_deltaRayLengthThresholdSquared.at(pCaloHit->GetHitType()))
    {
        return pParentMC;
    }

    // Now have a delta ray that we would like to cluster but only the hits that are not overlapping with the parent particle
    float parentWeight{std::numeric_limits<float>::lowest()};
    const MCParticleWeightMap &weightMap{pCaloHit->GetMCParticleWeightMap()};
    for (const auto &[pContributingMC, weight] : weightMap)
    {
        if (pContributingMC == pParentMC)
        {
            parentWeight = weight;
            break;
        }
    }
    if (parentWeight > m_deltaRayParentWeightThreshold)
    {
        return pParentMC;
    }
    return pMC;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventClusterValidationAlgorithm::ApplyMCParticleMinSumHits(std::map<const CaloHit *const, CaloHitParents> &hitParents) const
{
    if (m_minMCHitsPerView == 0)
    {
        return;
    }

    std::map<const MCParticle *const, int> mcSumHits;
    for (const auto &[pCaloHit, parents] : hitParents)
    {
        mcSumHits[parents.m_pMainMC]++;
    }

    for (auto it = hitParents.begin(); it != hitParents.end();)
    {
        if (mcSumHits.at(it->second.m_pMainMC) < m_minMCHitsPerView)
        {
            it = hitParents.erase(it);
        }
        else
        {
            it++;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventClusterValidationAlgorithm::GetClusterMainMC(std::map<const CaloHit *const, CaloHitParents> &hitParents) const
{
    std::map<const Cluster *const, std::map<const MCParticle *const, int>> clusterMainMCNHits;
    for (const auto &[pCaloHit, parents] : hitParents)
    {
        const Cluster *const pCluster{parents.m_pCluster};
        const MCParticle *const pMC{parents.m_pMainMC};
        clusterMainMCNHits[pCluster][pMC]++;
    }

    std::map<const Cluster *const, const MCParticle *const> clusterMainMCs;
    for (const auto &[pCluster, mainMCNHits] : clusterMainMCNHits)
    {
        const MCParticle *pClusterMainMC{nullptr};
        int maxHits{0};
        for (const auto &[pMC, nHits] : mainMCNHits)
        {
            if (nHits > maxHits)
            {
                pClusterMainMC = pMC;
                maxHits = nHits;
            }
            else if (nHits == maxHits) // tie-breaker
            {
                if (LArMCParticleHelper::SortByMomentum(pMC, pClusterMainMC))
                {
                    pClusterMainMC = pMC;
                }
            }
        }
        clusterMainMCs.insert({pCluster, pClusterMainMC});
    }

    for (auto &[pCaloHit, parents] : hitParents)
    {
        parents.m_pClusterMainMC = clusterMainMCs.at(parents.m_pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

int EventClusterValidationAlgorithm::HandleNullClusterHits(std::map<const CaloHit *const, CaloHitParents> &hitParents) const
{
    int nNullClusterHits{0};
    for (auto it = hitParents.begin(); it != hitParents.end();)
    {
        if (!it->second.m_pCluster)
        {
            if (m_dropNullClusterHits)
            {
                it = hitParents.erase(it);
            }
            else
            {
                it++;
            }
            nNullClusterHits++;
        }
        else
        {
            it++;
        }
    }

    return nNullClusterHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventClusterValidationAlgorithm::GetClusterMetrics(
    const std::map<const CaloHit *const, CaloHitParents> &hitParents, ClusterMetrics &metrics) const
{
    std::map<const MCParticle *const, std::set<const CaloHit *>> trueClusters;
    std::map<const Cluster *const, std::set<const CaloHit *>> recoClusters;
    std::map<const Cluster *const, const MCParticle *const> recoClusterToMainMC;
    for (const auto &[pCaloHit, parents] : hitParents)
    {
        if (trueClusters.find(parents.m_pMainMC) == trueClusters.end())
        {
            trueClusters.insert({parents.m_pMainMC, { pCaloHit }});
        }
        else
        {
            trueClusters.at(parents.m_pMainMC).insert(pCaloHit);
        }
        if (recoClusters.find(parents.m_pCluster) == recoClusters.end())
        {
            recoClusters.insert({parents.m_pCluster, { pCaloHit }});
            recoClusterToMainMC.insert({parents.m_pCluster, parents.m_pClusterMainMC});
        }
        else
        {
            recoClusters.at(parents.m_pCluster).insert(pCaloHit);
        }
    }

    // Completeness
    int nHitsTrueAll{0}, nHitsTrueShower{0}, nHitsTrueTrack{0};
    double completenessSumAll{0.}, completenessSumShower{0.}, completenessSumTrack{0.};
    int nTrueClustersShower{0}, nTrueClustersTrack{0};
    for (const auto &[pMC, hits] : trueClusters)
    {
        int maxIntersectionCard{0};
        for (const auto &[pCluster, hitsInner] : recoClusters)
        {
            std::set<const CaloHit *> intersection;
            std::set_intersection(
                hits.begin(), hits.end(), hitsInner.begin(), hitsInner.end(), std::inserter(intersection, intersection.begin()));
            const int intersectionCard{static_cast<int>(intersection.size())};
            if (intersectionCard > maxIntersectionCard)
            {
                maxIntersectionCard = intersectionCard;
            }
        }

        nHitsTrueAll += static_cast<int>(hits.size());
        double summand{static_cast<double>(maxIntersectionCard)};
        if (!m_hitWeightedPurityCompleteness) {
            summand /= static_cast<double>(hits.size());
        }
        completenessSumAll += summand;

        if (this->IsEM(pMC))
        {
            nHitsTrueShower += static_cast<int>(hits.size());
            nTrueClustersShower++;
            completenessSumShower += summand;
        }
        else
        {
            nHitsTrueTrack += static_cast<int>(hits.size());
            nTrueClustersTrack++;
            completenessSumTrack += summand;
        }
    }

    // Purity
    int nHitsMatchedAll{0}, nHitsMatchedShower{0}, nHitsMatchedTrack{0};
    double puritySumAll{0.}, puritySumShower{0.}, puritySumTrack{0.};
    int nRecoClustersShower{0}, nRecoClustersTrack{0};
    for (const auto &[pCluster, hits] : recoClusters)
    {
        const MCParticle *const pMC{recoClusterToMainMC.at(pCluster)};

        int maxIntersectionCard{0};
        for (const auto &[pMCInner, hitsInner] : trueClusters)
        {
            std::set<const CaloHit *> intersection;
            std::set_intersection(
                hits.begin(), hits.end(), hitsInner.begin(), hitsInner.end(), std::inserter(intersection, intersection.begin()));
            const int intersectionCard{static_cast<int>(intersection.size())};
            if (intersectionCard > maxIntersectionCard)
            {
                maxIntersectionCard = intersectionCard;
            }
        }

        nHitsMatchedAll += static_cast<int>(hits.size());
        double summand{static_cast<double>(maxIntersectionCard)};
        if (!m_hitWeightedPurityCompleteness) {
            summand /= static_cast<double>(hits.size());
        }
        puritySumAll += summand;

        if (this->IsEM(pMC))
        {
            nHitsMatchedShower += static_cast<int>(hits.size());
            nRecoClustersShower++;
            puritySumShower += summand;
        }
        else
        {
            nHitsMatchedTrack += static_cast<int>(hits.size());
            nRecoClustersTrack++;
            puritySumTrack += summand;
        }
    }

    auto computeMetric = [&](const double sum, const int nHits, const int nClusters) {
        int denom{m_hitWeightedPurityCompleteness ? nHits : nClusters};
        return denom > 0 ? sum / static_cast<double>(denom) : -1.;
    };
    metrics.m_completeness = computeMetric(completenessSumAll, nHitsTrueAll, trueClusters.size());
    metrics.m_showerCompleteness = computeMetric(completenessSumShower, nHitsTrueShower, nTrueClustersShower);
    metrics.m_trackCompleteness = computeMetric(completenessSumTrack, nHitsTrueTrack, nTrueClustersTrack);
    metrics.m_purity = computeMetric(puritySumAll, nHitsMatchedAll, recoClusters.size());
    metrics.m_showerPurity = computeMetric(puritySumShower, nHitsMatchedShower, nRecoClustersShower);
    metrics.m_trackPurity = computeMetric(puritySumTrack, nHitsMatchedTrack, nRecoClustersTrack);

    metrics.m_nShowerTrueHits = nHitsTrueShower;
    metrics.m_nTrackTrueHits = nHitsTrueTrack;
    metrics.m_nRecoClusters = recoClusters.size();
    metrics.m_nShowerRecoClusters = nRecoClustersShower;
    metrics.m_nTrackRecoClusters = nRecoClustersTrack;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventClusterValidationAlgorithm::GetMatchedParticleMetrics(
    const std::map<const CaloHit *const, CaloHitParents> &hitParents, MatchedParticleMetrics &metrics) const
{
    std::map<const MCParticle *const, const Cluster *> mcMatchedCluster;
    std::map<const MCParticle *const, int> mcMatchedClusterCorrectHits, mcMatchedClusterTotalHits, mcNTrueHits;
    for (const auto &[pCaloHit, parents] : hitParents)
    {
        if (mcMatchedCluster.find(parents.m_pMainMC) == mcMatchedCluster.end())
        {
            mcMatchedCluster.insert({parents.m_pMainMC, nullptr});
            mcMatchedClusterCorrectHits.insert({parents.m_pMainMC, 0});
            mcMatchedClusterTotalHits.insert({parents.m_pMainMC, 0});
            mcNTrueHits.insert({parents.m_pMainMC, 0});
        }
        mcNTrueHits.at(parents.m_pMainMC)++;
    }

    std::map<const Cluster *const, std::map<const MCParticle *const, int>> clusterMCNHits;
    std::map<const Cluster *const, int> clusterNHits;
    for (const auto &[pCaloHit, parents] : hitParents)
    {
        if (parents.m_pCluster)
        {
            clusterMCNHits[parents.m_pCluster][parents.m_pMainMC]++;
            clusterNHits[parents.m_pCluster]++;
        }
    }

    auto isBetterMatch =
        [&mcMatchedCluster, &mcMatchedClusterCorrectHits, &mcMatchedClusterTotalHits]
        (const MCParticle *const pMC, const int nCorrectHits, const int nTotalHits, const Cluster *const pCluster)
    {
        return
            !mcMatchedCluster.at(pMC) ||                                                // No competitor
            nCorrectHits > mcMatchedClusterCorrectHits.at(pMC) ||                       // More matched hits
            (nCorrectHits == mcMatchedClusterCorrectHits.at(pMC) &&                     // Need to do a tie-breaker
                (nTotalHits < mcMatchedClusterTotalHits.at(pMC) ||                      // Purity tie-breaker
                LArClusterHelper::SortByPosition(pCluster, mcMatchedCluster.at(pMC)))); // Arbitrary tie-breaker
    };

    if (!m_maximalMatching)
    {
        ClusterList seenClusters;
        for (const auto &[pCaloHit, parents] : hitParents)
        {
            const Cluster *const pCluster{parents.m_pCluster};
            const MCParticle *const pMatchedMC{parents.m_pClusterMainMC};

            if (!pCluster)
            {
                continue; // Not letting the null cluster be matched
            }

            if (std::find(seenClusters.begin(), seenClusters.end(), pCluster) != seenClusters.end())
            {
                continue;
            }
            seenClusters.emplace_back(pCluster);

            int nTotalHits{clusterNHits.at(pCluster)};
            int nCorrectHits{clusterMCNHits.at(pCluster).at(pMatchedMC)};

            if (isBetterMatch(pMatchedMC, nCorrectHits, nTotalHits, pCluster))
            {
                mcMatchedCluster.at(pMatchedMC) = pCluster;
                mcMatchedClusterCorrectHits.at(pMatchedMC) = nCorrectHits;
                mcMatchedClusterTotalHits.at(pMatchedMC) = nTotalHits;
            }
        }
    }
    else
    {
        std::map<const Cluster *const, MCParticleList> clusterOrderedMCs;
        for (const auto &[pCluster, mcNHits] : clusterMCNHits)
        {
            MCParticleList orderedMCs;
            for (const auto &[pMC, nHits] : mcNHits)
            {
                orderedMCs.emplace_back(pMC);
            }
            orderedMCs.sort([&mcNHits](const MCParticle *const pA, const MCParticle *const pB) { return mcNHits.at(pA) > mcNHits.at(pB); });
            clusterOrderedMCs.insert({pCluster, orderedMCs});
        }

        std::set<const Cluster *> matchedClusters;
        bool newMatch{true};
        while (newMatch)
        {
            newMatch = false;
            for (auto it = clusterOrderedMCs.begin(); it != clusterOrderedMCs.end(); )
            {
                const Cluster *const pCluster{it->first};
                // Exhausted this cluster's MC contributions, or the cluster was already matched in the previous loop
                if (it->second.empty() || matchedClusters.find(pCluster) != matchedClusters.end())
                {
                    it = clusterOrderedMCs.erase(it);
                    continue;
                }
                const MCParticle *const pMC{it->second.front()}; it->second.pop_front();
                const int nCorrectHits{clusterMCNHits.at(pCluster).at(pMC)};
                const int nTotalHits{clusterNHits.at(pCluster)};

                if (isBetterMatch(pMC, nCorrectHits, nTotalHits, pCluster))
                {
                    newMatch = true;
                    matchedClusters.erase(mcMatchedCluster.at(pMC));
                    matchedClusters.insert(pCluster);
                    mcMatchedCluster.at(pMC) = pCluster;
                    mcMatchedClusterCorrectHits.at(pMC) = nCorrectHits;
                    mcMatchedClusterTotalHits.at(pMC) = nTotalHits;
                }

                it++;
            }
        }
    }

    for (const auto &[pMC, pCluster] : mcMatchedCluster)
    {
        metrics.m_pdg.emplace_back(pMC->GetParticleId());
        metrics.m_causesShower.emplace_back(this->CausesShower(pMC, 0));
        metrics.m_isPrimary.emplace_back(pMC->GetParentList().front()->IsRootParticle());
        metrics.m_trueEnergy.emplace_back(pMC->GetEnergy());
        metrics.m_nTrueHits.emplace_back(mcNTrueHits.at(pMC));
        metrics.m_nMatchedCorrectHits.emplace_back(mcMatchedClusterCorrectHits.at(pMC));
        metrics.m_nMatchedTotalHits.emplace_back(mcMatchedClusterTotalHits.at(pMC));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventClusterValidationAlgorithm::CalcRandIndex(
    std::map<const CaloHit *const, CaloHitParents> &hitParents, ClusterMetrics &metrics) const
{
    // Do a perfect shower growing of the reco clusters
    std::map<const CaloHit *const, const Cluster *const> hitMergeTarget;
    if (m_mergeShowerClustersForRandIndex)
    {
        std::map<const MCParticle *const, const Cluster *const> mcMergeTarget;
        for (const auto &[pCaloHit, parents] : hitParents)
        {
            if (!this->IsEM(parents.m_pClusterMainMC))
            {
                continue;
            }
            // The first cluster we come across will used as the target for the merging within the shower
            if (mcMergeTarget.find(parents.m_pClusterMainMC) == mcMergeTarget.end())
            {
                mcMergeTarget.insert({parents.m_pClusterMainMC, parents.m_pCluster});
            }
            hitMergeTarget.insert({pCaloHit, mcMergeTarget.at(parents.m_pClusterMainMC)});
        }
    }

    // Fill contingency table(s)
    ContingencyTable<const Cluster *const, const MCParticle *const> cTable;
    ContingencyTable<const Cluster *const, const MCParticle *const> cTableShower;
    ContingencyTable<const Cluster *const, const MCParticle *const> cTableTrack;

    std::set<const Cluster *> recoClusters, trackRecoClusters, showerRecoClusters;
    std::set<const MCParticle *> trueClusters, trackTrueClusters, showerTrueClusters;

    for (const auto &[pCaloHit, parents] : hitParents)
    {
        const MCParticle *const pMC = parents.m_pMainMC;
        const Cluster *pCluster = parents.m_pCluster;

        if (hitMergeTarget.find(pCaloHit) != hitMergeTarget.end())
        {
            pCluster = hitMergeTarget.at(pCaloHit);
        }

        recoClusters.insert(pCluster);
        trueClusters.insert(pMC);
        cTable[pCluster][pMC]++;

        if (m_onlyRandIndex)
        {
            continue;
        }

        if (this->IsEM(pMC))
        {
            showerRecoClusters.insert(pCluster);
            showerTrueClusters.insert(pMC);
            cTableShower[pCluster][pMC]++;
        }
        else
        {
            trackRecoClusters.insert(pCluster);
            trackTrueClusters.insert(pMC);
            cTableTrack[pCluster][pMC]++;
        }
    }

    if (m_visualize && m_mergeShowerClustersForRandIndex) // This is only worth seeing if the cheated merging has occurred
    {
        std::cout << "---- Drawing reco clusters used in Rand Index calculation...\n";
        this->VisualizeRandIndexRecoClusters(hitParents, hitMergeTarget);
    }

    metrics.m_ari = LArMonitoringHelper::CalcRandIndex(cTable);
    metrics.m_nAriRecoClusters = static_cast<int>(recoClusters.size());
    metrics.m_nTrueClusters = static_cast<int>(trueClusters.size());
    if (m_onlyRandIndex)
    {
        return;
    }
    metrics.m_showerAri = LArMonitoringHelper::CalcRandIndex(cTableShower);
    metrics.m_nShowerAriRecoClusters = static_cast<int>(showerRecoClusters.size());
    metrics.m_nShowerTrueClusters = static_cast<int>(showerTrueClusters.size());
    metrics.m_trackAri = LArMonitoringHelper::CalcRandIndex(cTableTrack);
    metrics.m_nTrackAriRecoClusters = static_cast<int>(trackRecoClusters.size());
    metrics.m_nTrackTrueClusters = static_cast<int>(trackTrueClusters.size());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventClusterValidationAlgorithm::SetBranches(
    [[maybe_unused]] const ClusterMetrics &clusterMetrics,
    [[maybe_unused]] const MatchedParticleMetrics &matchedParticleMetrics,
    [[maybe_unused]] const int view) const
{
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "event", m_eventNumber - 1));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "view", view));

    // Just the all-hits ARI with some values that might be needed for cuts
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "n_hits", clusterMetrics.m_nHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "n_hits_null", clusterMetrics.m_nHitsNullCluster));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "adjusted_rand_idx", clusterMetrics.m_ari));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "n_true_clusters", clusterMetrics.m_nTrueClusters));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "n_ari_reco_clusters", clusterMetrics.m_nAriRecoClusters));

    if (m_onlyRandIndex)
    {
        PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName));
        return;
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "shower_n_hits", clusterMetrics.m_nShowerTrueHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "track_n_hits", clusterMetrics.m_nTrackTrueHits));
    
    // The shower/track ARI with values needed for cuts
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "shower_adjusted_rand_idx", clusterMetrics.m_showerAri));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "shower_n_true_clusters", clusterMetrics.m_nShowerTrueClusters));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "shower_n_ari_reco_clusters", clusterMetrics.m_nShowerAriRecoClusters));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "track_adjusted_rand_idx", clusterMetrics.m_trackAri));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "track_n_true_clusters", clusterMetrics.m_nTrackTrueClusters));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "track_n_ari_reco_clusters", clusterMetrics.m_nTrackAriRecoClusters));

    // The cluster purity/completenesses with values needed for cuts (already have the number of true clusters from ARI branches)
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "purity", clusterMetrics.m_purity));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "completeness", clusterMetrics.m_completeness));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "n_reco_clusters", clusterMetrics.m_nRecoClusters));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "shower_purity", clusterMetrics.m_showerPurity));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "shower_completeness", clusterMetrics.m_showerCompleteness));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "shower_n_reco_clusters", clusterMetrics.m_nShowerRecoClusters));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "track_purity", clusterMetrics.m_trackPurity));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "track_completeness", clusterMetrics.m_trackCompleteness));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "track_n_reco_clusters", clusterMetrics.m_nTrackRecoClusters));

    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName));

    // Matched particle metrics
    if (m_matchedParticleMetrics)
    {
        for (int i = 0; i < static_cast<int>(matchedParticleMetrics.m_pdg.size()); i++)
        {
            PANDORA_MONITORING_API(SetTreeVariable(
                this->GetPandora(), m_treeName + "_matching", "event", m_eventNumber - 1));
            PANDORA_MONITORING_API(SetTreeVariable(
                this->GetPandora(), m_treeName + "_matching", "view", static_cast<int>(view)));
            PANDORA_MONITORING_API(SetTreeVariable(
                this->GetPandora(), m_treeName + "_matching", "pdg", matchedParticleMetrics.m_pdg.at(i)));
            PANDORA_MONITORING_API(SetTreeVariable(
                this->GetPandora(), m_treeName + "_matching", "causes_shower", matchedParticleMetrics.m_causesShower.at(i)));
            PANDORA_MONITORING_API(SetTreeVariable(
                this->GetPandora(), m_treeName + "_matching", "is_primary", matchedParticleMetrics.m_isPrimary.at(i)));
            PANDORA_MONITORING_API(SetTreeVariable(
                this->GetPandora(), m_treeName + "_matching", "energy", matchedParticleMetrics.m_trueEnergy.at(i)));
            PANDORA_MONITORING_API(SetTreeVariable(
                this->GetPandora(), m_treeName + "_matching", "n_true_hits", matchedParticleMetrics.m_nTrueHits.at(i)));
            PANDORA_MONITORING_API(SetTreeVariable(
                this->GetPandora(), m_treeName + "_matching", "n_correct_matched_hits", matchedParticleMetrics.m_nMatchedCorrectHits.at(i)));
            PANDORA_MONITORING_API(SetTreeVariable(
                this->GetPandora(), m_treeName + "_matching", "n_total_matched_hits", matchedParticleMetrics.m_nMatchedTotalHits.at(i)));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName + "_matching"));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventClusterValidationAlgorithm::VisualizeTargetClusters(std::map<const CaloHit *const, CaloHitParents> &hitParents) const
{
    std::map<const MCParticle *const, CaloHitList> mcToCaloHits;
    for (const auto &[pCaloHit, parents] : hitParents)
    {
        mcToCaloHits[parents.m_pMainMC].push_back(pCaloHit);
    }

    int color {1}; // 0 is white
    for ([[maybe_unused]] const auto &[pMC, caloHits] : mcToCaloHits)
    {
        PANDORA_MONITORING_API(
            VisualizeCaloHits(this->GetPandora(),
                &caloHits,
                "From MC: " + std::to_string(pMC->GetParticleId()) + " - " + std::to_string(pMC->GetEnergy()),
                static_cast<Color>(color++)));
        if (color > 30) // start getting into the AUTO colors territory
        {
            color = 1;
        }
    }
    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventClusterValidationAlgorithm::VisualizeRandIndexRecoClusters(
    std::map<const CaloHit *const, CaloHitParents> &hitParents,
    std::map<const CaloHit *const, const Cluster *const> &hitMergeTargets) const
{
    std::map<const Cluster *const, CaloHitList> clusterToCaloHits;
    std::map<const Cluster *const, const MCParticle *const> clusterToMainMC;
    for (const auto &[pCaloHit, parents] : hitParents)
    {
        if (clusterToMainMC.find(parents.m_pCluster) == clusterToMainMC.end())
        {
            clusterToMainMC.insert({parents.m_pCluster, parents.m_pClusterMainMC});
        }

        if (hitMergeTargets.find(pCaloHit) == hitMergeTargets.end())
        {
            clusterToCaloHits[parents.m_pCluster].push_back(pCaloHit);
        }
        else
        {
            clusterToCaloHits[hitMergeTargets.at(pCaloHit)].push_back(pCaloHit);
        }
    }

    int color {1}; // 0 is white
    for (const auto &[pCluster, caloHits] : clusterToCaloHits)
    {
        [[maybe_unused]] const MCParticle *const pMC{clusterToMainMC.at(pCluster)};
        PANDORA_MONITORING_API(
            VisualizeCaloHits(
                this->GetPandora(),
                &caloHits,
                "From Reco: " + std::to_string(pMC->GetParticleId()) + " - " + std::to_string(pMC->GetEnergy()),
                static_cast<Color>(color++)));
        if (color > 30) // start getting into the AUTO colors territory
        {
            color = 1;
        }
    }
    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventClusterValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "FileName", m_fileName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TreeName", m_treeName));

    std::vector<std::string> caloHitListNames;
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "CaloHitListNames", caloHitListNames));
    if (!caloHitListNames.empty())
    {
        m_caloHitListNames = caloHitListNames;
    }
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "ClusterListNames", m_clusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinMCHitsPerView", m_minMCHitsPerView));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "OnlyRandIndex", m_onlyRandIndex));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FoldShowers", m_foldShowers));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "HandleDeltaRays", m_handleDeltaRays));
    if (m_handleDeltaRays)
    {
        m_deltaRayLengthThresholdSquared = {
            { TPC_VIEW_U, static_cast<float>(pow(LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_U), 2.)) },
            { TPC_VIEW_V, static_cast<float>(pow(LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_V), 2.)) },
            { TPC_VIEW_W, static_cast<float>(pow(LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_W), 2.)) }
        };
    }
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MergeShowerClustersForRandIndex", m_mergeShowerClustersForRandIndex));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MatchedParticleMetrics", m_matchedParticleMetrics));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "DropNullClusterHits", m_dropNullClusterHits));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "HitWeightedPurityCompleteness", m_hitWeightedPurityCompleteness));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaximalMatching", m_maximalMatching));
    if (m_maximalMatching && !m_matchedParticleMetrics)
    {
        return STATUS_CODE_INVALID_PARAMETER;
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualize", m_visualize));


    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
