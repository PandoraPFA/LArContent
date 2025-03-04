/**
 *  @file   larpandoracontent/LArHelpers/LArMonitoringHelper.h
 *
 *  @brief  Header file for the lar monitoring helper helper class.
 *
 *  $Log: $
 */
#ifndef LAR_MONITORING_HELPER_H
#define LAR_MONITORING_HELPER_H 1

#include "Pandora/PandoraInternal.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content
{

/**
 *  @brief  LArMonitoringHelper class
 */
class LArMonitoringHelper
{
public:
    /**
     *  @brief  Count the number of calo hits, in a provided list, of a specified type
     *
     *  @param  hitType the hit type
     *  @param  caloHitList the calo hit list
     *
     *  @return the number of calo hits of the specified type
     */
    static unsigned int CountHitsByType(const pandora::HitType hitType, const pandora::CaloHitList &caloHitList);

    /**
     *  @brief  Order input MCParticles by their number of hits.
     *
     *  @param  selectedMCParticleToGoodHitsMaps the input vector of mappings from selected reconstructable MCParticles to their good hits
     *  @param  orderedMCParticleVector the output vector of ordered MCParticles
     */
    static void GetOrderedMCParticleVector(const LArMCParticleHelper::MCContributionMapVector &selectedMCParticleToGoodHitsMaps,
        pandora::MCParticleVector &orderedMCParticleVector);

    /**
     *  @brief  Order input Pfos by their number of hits.
     *
     *  @param  pfoToReconstructable2DHitsMap the input vector of mappings from Pfos to their reconstructable hits
     *  @param  orderedPfoVector the output vector of ordered Pfos
     */
    static void GetOrderedPfoVector(const LArMCParticleHelper::PfoContributionMap &pfoToReconstructable2DHitsMap, pandora::PfoVector &orderedPfoVector);

    /**
     *  @brief  Print details of selected MCParticles to the terminal in a table.
     *
     *  @param  selectedMCParticleToGoodHitsMap the input mapping from selected reconstructable MCParticles to their good hits
     *  @param  orderedMCParticleVector the input vector of ordered MCParticles
     */
    static void PrintMCParticleTable(const LArMCParticleHelper::MCContributionMap &selectedMCParticleToGoodHitsMaps,
        const pandora::MCParticleVector &orderedMCParticleVector);

    /**
     *  @brief  Print details of input Pfos to the terminal in a table.
     *
     *  @param  pfoToReconstructable2DHitsMap the input vector of mappings from Pfos to their reconstructable hits
     *  @param  orderedPfoVector the input vector of ordered Pfos
     */
    static void PrintPfoTable(const LArMCParticleHelper::PfoContributionMap &pfoToReconstructable2DHitsMap, const pandora::PfoVector &orderedPfoVector);

    /**
     *  @brief  Print the shared good hits between all Pfos and MCParticles
     *
     *  @param  orderedPfoVector the input vector of ordered Pfos
     *  @param  orderedMCParticleVector the input vector of ordered MCParticles
     *  @param  mcParticleToPfoHitSharingMap the output mapping from selected reconstructable MCParticles to Pfos and the number hits shared
     *  @param  nMatches the maximum number of Pfo matches to show
     */
    static void PrintMatchingTable(const pandora::PfoVector &orderedPfoVector, const pandora::MCParticleVector &orderedMCParticleVector,
        const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcParticleToPfoHitSharingMap, const unsigned int nMatches);

    // -- Rand index helpers

    /**
     *  @brief  Calculate the adjusted Rand Index for a given contingency table that summarises two clusterings A and B.
     *          Adjusted Rand Index is a measure of the similarity of two clusterings that is bounded above by 1 (identical)
     *          and is 0 when the two clusterings are as similar as two random clusterings with the same number of clusters.
     *          Reference - https://link.springer.com/article/10.1007/BF01908075
     *
     *  @param[in] cTable Contingency table as a map of the object that defines the clustering A (e.g. pandora::Cluster)
     *                    to a map of the object that defines the clustering B (e.g. pandora::MCParticle)
     *                    to the value of the table at this entry (which is the intersection of the two clusters)
     */
    template<typename Ti, typename Tj>
    static float CalcRandIndex(const std::map<const Ti, std::map<const Tj, int>> &cTable);

    /**
     *  @brief  Calculate the adjusted Rand Index for the clustering defined by MCPartices and by CaloHits.
     *          Adjusted Rand Index is a measure of the similarity of two clusterings that is bounded above by 1 (identical)
     *          and is 0 when the two clusterings are as similar as two random clusterings with the same number of clusters.
     *          Reference - https://link.springer.com/article/10.1007/BF01908075
     *
     *  @param[in] caloHits List of hits
     *  @param[in] clusters List of clusters
     */
    static float CalcRandIndex(const pandora::CaloHitList &caloHits, const pandora::ClusterList &clusters);

    /**
     *  @brief  Fill the contingency table for a set of CaloHits partitioned by Cluster and parent MCParticle.
     *
     *  @param[in]  caloHits List of hits to be considered
     *  @param[in]  clusters List of clusters that contain the hits
     *  @param[out] cTable   Contingency table as a map of Cluster to a map of MCParticle to the value of the table at this entry
     *                       which is the intersection of clusterings defined by each object
     */
    static void FillContingencyTable(const pandora::CaloHitList &caloHits,
                                     const pandora::ClusterList &clusters,
                                     std::map<const pandora::Cluster *const, std::map<const pandora::MCParticle *const, int>> &cTable);

    // --
};

} // namespace lar_content

#endif // #ifndef LAR_MONITORING_HELPER_H
