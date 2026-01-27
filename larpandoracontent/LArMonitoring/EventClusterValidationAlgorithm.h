/**
 *  @file   larpandoracontent/LArMonitoring/EventClusterValidationAlgorithm.h
 *
 *  @brief  Header file of the event-level cluster validation. Calculate metrics that aim to quantify the quality of 2D clusters at
 *          the level of a single view.
 */
#ifndef LAR_EVENT_CLUSTER_VALIDATION_ALGORITHM_H
#define LAR_EVENT_CLUSTER_VALIDATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  EventClusterValidationAlgorithm class
 */
class EventClusterValidationAlgorithm : public pandora::Algorithm
{
private:
    struct ClusterMetrics
    {
        ClusterMetrics();

        double m_purity;
        double m_showerPurity;
        double m_trackPurity;
        double m_completeness;
        double m_showerCompleteness;
        double m_trackCompleteness;
        double m_ari;
        double m_showerAri;
        double m_trackAri;
        int m_nHits;
        int m_nHitsNullCluster;
        int m_nShowerTrueHits;
        int m_nTrackTrueHits;
        int m_nRecoClusters;
        int m_nShowerRecoClusters;
        int m_nTrackRecoClusters;
        int m_nTrueClusters;
        int m_nShowerTrueClusters;
        int m_nTrackTrueClusters;
        int m_nAriRecoClusters;
        int m_nShowerAriRecoClusters;
        int m_nTrackAriRecoClusters;
    };

    struct MatchedParticleMetrics
    {
        MatchedParticleMetrics();

        std::vector<int> m_pdg;
        std::vector<int> m_causesShower;
        std::vector<int> m_isPrimary;
        std::vector<float> m_trueEnergy;
        std::vector<int> m_nTrueHits;
        std::vector<int> m_nMatchedCorrectHits;
        std::vector<int> m_nMatchedTotalHits;
    };

public:
    struct CaloHitParents
    {
        CaloHitParents();

        const pandora::MCParticle *m_pMainMC;
        const pandora::Cluster *m_pCluster;
        const pandora::MCParticle *m_pClusterMainMC;
    };

    /**
    *  @brief  Default constructor
    */
    EventClusterValidationAlgorithm();

    /**
     *  @brief  Destructor saves the validation TTree along with making and saving a metadata TTree
     */
    ~EventClusterValidationAlgorithm();

private:
    template <typename Ti, typename Tj>
    using ContingencyTable = std::map<Ti, std::map<Tj, int>>;

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief Calculate purity and completeness in this view for all/track/shower reco and true clusters.
     *         Completeness(purity) sum over true(reco) clusters the maximum intersection with any of the reco(true) clusters.
     *         According to m_hitWeightedPurityCompleteness, this is either normalised by the number of hits in the true(reco) cluster and
     *         averaged over the total true(reco) clusters, or not normalised and averaged over the total hits.
     *
     *  @param[in]  hitParents Map of hits being considered to the cluster/MC particle they belong to
     *  @param[out] metrics    Output Metrics for the clustering in this view
     */
    void GetClusterMetrics(const std::map<const pandora::CaloHit *const, CaloHitParents> &hitParents, ClusterMetrics &metrics) const;

    /**
     *  @brief Retrieve the metrics for every MC Particle in a view by matching clusters to a true MC particle
     *
     *  @param[in]  hitParents Map of hits being considered to the cluster/MC particle they belong to
     *  @param[out] metrics    Output metrics for the matched MC particles in this view
     */
    void GetMatchedParticleMetrics(const std::map<const pandora::CaloHit *const, CaloHitParents> &hitParents, MatchedParticleMetrics &metrics) const;

    /**
     *  @brief Calculate Rand Index for the reco clusters with the true clusters for all/track/shower hits.
     *         According to m_mergeShowerClustersForRandIndex, the shower clusters may be merged using cheating prior to the calculation.
     *         (Ref. for Adjusted Rand Index: https://link.springer.com/article/10.1007/BF01908075).
     *
     *  @param[in]  hitParents Map of hits being considered to the cluster/MC particle they belong to
     *  @param[out] metrics    Output metrics for the clustering in this view
     */
    void CalcRandIndex(std::map<const pandora::CaloHit *const, CaloHitParents> &hitParents, ClusterMetrics &metrics) const;

    /**
     *  @brief Find the cluster and MC particle each hit belongs to
     *
     *  @param[in]  caloHits   List of hits
     *  @param[in]  clusters   List of clusters
     *  @param[out] hitParents Map of hits to the Cluster/MC particle they belong to
     */
    void GetHitParents(const pandora::CaloHitList &caloHits, const pandora::ClusterList &clusters,
        std::map<const pandora::CaloHit *const, CaloHitParents> &hitParents) const;

    /**
     *  @brief Find the main MC particle for each cluster, this is the MC particle that contributes the most hits
     *
     *  @param[in,out] hitParents Map of hits to the cluster/MC particle they belong to
     */
    void GetClusterMainMC(std::map<const pandora::CaloHit *const, CaloHitParents> &hitParents) const;

    /**
     *  @brief Tries to identify and deal with impossible-to-cluster-correctly delta ray hits by assigning the hit to the
     *         parent particle's cluster if some conditions are met.
     *
     *  @param[in] pCaloHit The hit
     *  @param[in] pMC      The main MC particle of the hit
     *
     *  @return The MC particle the hit should be assigned to, this will either be the inputted MC particle or the parent MC particle
     */
    const pandora::MCParticle *FoldPotentialDeltaRayTo(const pandora::CaloHit *const pCaloHit, const pandora::MCParticle *const pMC) const;

    /**
     *  @brief Finds the ancestor that a child MC particle should be associated with to roll-up an EM shower
     *
     *  @param[in] pMC The MC particle
     *
     *  @return The ancestor MC particle, this will be the inputted MC particle for an MC particle that should not be rolled-up
     */
    const pandora::MCParticle *FoldMCTo(const pandora::MCParticle *const pMC) const;

    /**
     *  @brief Recursive function to check descendent particles for signature of a shower (e -> gamma -> e).
     *         Used to identify delta-rays that shower and photons that only undergo compton scatters leaving diffuse hits.
     *
     *  @param[in] pMC                  The MC particle to check descendents of
     *  @param[in] nDescendentElectrons The number of electrons seen while descending from the original photon MC particle
     *
     *  @return Flag to indicate if more than 1 electron was seen while descending any of the descendent particle association paths
     */
    bool CausesShower(const pandora::MCParticle *const pMC, int nDescendentElectrons) const;

    /**
     *  @brief Check if an MC particle is electron or photon,
     *
     *  @param[in] pMC he MC particle
     *
     *  @return Flag to indicate if the MC particle if EM
     */
    bool IsEM(const pandora::MCParticle *const pMC) const;

    /**
     *  @brief Draw the true clusters being compared to
     *
     *  @param[in] hitParents Map of hits to the cluster/MC particle they belong to
     */
    void VisualizeTargetClusters(std::map<const pandora::CaloHit *const, CaloHitParents> &hitParents) const;

    /**
     *  @brief Draw the reco clusters used in the Rand Index calculation. These will be different from the clusters inputted to the
     *         algorithm if MergeShowerClustersForRandIndex is true.
     *
     *  @param[in] hitParents Map of hits to the cluster/MC particle they belong to
     *  @param[in] hitParents Map of hits to the cluster they merge with for the Rand Index Calculation.
     */
    void VisualizeRandIndexRecoClusters(std::map<const pandora::CaloHit *const, CaloHitParents> &hitParents,
        std::map<const pandora::CaloHit *const, const pandora::Cluster *const> &hitMergeTargets) const;

    /**
     *  @brief Erase hits associated with an MC particle that does meet a minimum number of hits in the view
     *
     *  @param[in,out] hitParents Map of hits to the cluster/MC particle they belong to
     */
    void ApplyMCParticleMinSumHits(std::map<const pandora::CaloHit *const, CaloHitParents> &hitParents) const;

    /**
     *  @brief Handles hits that do not have a reco cluster, this is usually happens when left behind by the 3D hit creation.
     *         Will count the number of hits with no reco cluster and, according to m_dropNullClusterHits, may erase them.
     *
     *  @param[in,out] hitParents Map of hits to the cluster/MC particle they belong to
     */
    int HandleNullClusterHits(std::map<const pandora::CaloHit *const, CaloHitParents> &hitParents) const;

    /**
     *  @brief Update the branches of the TTree for this entry
     *
     *  @param[in] clusterMetrics         Output metrics for clusters in the view
     *  @param[in] matchedParticleMetrics Output metrics each particle matched to a cluster in the view
     *  @param[in] view                   The view
     */
    void SetBranches(const ClusterMetrics &clusterMetrics, const MatchedParticleMetrics &matchedParticleMetrics, const int view) const;

    // members
    int m_eventNumber;                                                  ///< To track the current event number
    std::map<pandora::HitType, float> m_deltaRayLengthThresholdSquared; ///< Threshold for defining small delta rays that will be folded to the parent particle
    float m_deltaRayParentWeightThreshold; ///< Threshold for weight contribution of parent particle for it take the delta ray's hit

    // members that may be set from xml
    std::string m_fileName;                      ///< The filename of the ROOT output file
    std::string m_treeName;                      ///< The name of the ROOT tree
    std::vector<std::string> m_caloHitListNames; ///< The names of the hit lists containing all 2D hits
    std::vector<std::string> m_clusterListNames; ///< The names of the lists of 2D clusters to process
    int m_minMCHitsPerView; ///< Threshold on total main MC particle hits in each view for consideration in metric calculations
    bool m_onlyRandIndex;   ///< Flag to only calculate the adjusted rand index over all particles
    bool m_foldShowers;     ///< Flag to fold shower MC particles to their leading shower MC particle
    bool m_handleDeltaRays; ///< Flag to fold short delta rays + ignore contributions from daughter electrons that overlap with their parent
    bool m_mergeShowerClustersForRandIndex; ///< Flag to merge shower-matched clusters into leading shower MC particle for rand index calculation, note this is only makes any sense for showers folded in the simulation or with m_foldShowers
    bool m_visualize;                       ///< Flag to display the target clustering derived from MC particles
    bool m_matchedParticleMetrics;          ///< Flag to calculate a set of high level clustering metrics by truth-matching clusters
    bool m_dropNullClusterHits;             ///< Flag to ignore any hits associated with the null cluster
    bool m_hitWeightedPurityCompleteness;   ///< Flag for cluster purity and completeness to be averaged over hits rather than clusters
    bool m_maximalMatching; ///< Flag for cluster truth-matching to try matching every MC particle to a cluster, rather than requiring clusters to a hit majority from the matched MC particle
};

} // namespace lar_content

#endif // LAR_EVENT_CLUSTER_VALIDATION_ALGORITHM_H
