/**
 *  @file   larpandoracontent/LArMonitoring/EventClusterValidationAlgorithm.h
 *
 *  @brief  Header file of the event-level cluster validation.
 *
 *  $Log: $
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

        std::vector<int> m_nRecoHits;
        std::vector<float> m_purities;
        std::vector<float> m_completenesses;
        int m_nHits;
        int m_nClusters;
        int m_nMainMCs;
    };

    enum ValidationType
    {
        ALL,
        SHOWER,
        TRACK
    };

public:
    struct CaloHitParents
    {
        CaloHitParents();

        const pandora::MCParticle *m_pMainMC;
        const pandora::Cluster *m_pCluster;
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
    typedef std::map<const pandora::Cluster *const, std::map<const pandora::MCParticle *const, int>> ContingencyTable;

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief Retrieve the metrics of every cluster in a view
     *
     *  @param[in]  hitParents Map of hits being considered to the Cluster/MCParticle they belong to
     *  @param[out] metrics    Metrics for the clusters in this view
     */
    void GetMetrics(const std::map<const pandora::CaloHit *const, CaloHitParents> &hitParents, ClusterMetrics &metrics) const;

    /**
     *  @brief Calculate Rand Index for the clusters with the clusters of true MCParticles.
     *          (Ref. for Adjusted Rand Index: https://link.springer.com/article/10.1007/BF01908075)
     *
     *  @param[in] hitParents Map of hits being considered to the Cluster/MCParticle they belong to
     */
    float CalcRandIndex(std::map<const pandora::CaloHit *const, CaloHitParents> &hitParents) const;

    /**
     *  @brief Find the cluster and MCParticle each hit belongs to
     *
     *  @param[in]  caloHits   List of hits
     *  @param[in]  clusters   List of clusters
     *  @param[out] hitParents Map of hits to the Cluster/MCParticle they belong to
     */
    void GetHitParents(const pandora::CaloHitList &caloHits, const pandora::ClusterList &clusters,
        std::map<const pandora::CaloHit *const, CaloHitParents> &hitParents) const;

    /**
     *  @brief Erase hits associated to an MCParticle that does meet a minimum number of hits in the view
     *
     *  @param[in,out] hitParents Map of hits to the Cluster/MCParticle they belong to
     */
    void ApplyMCParticleMinSumHits(std::map<const pandora::CaloHit *const, CaloHitParents> &hitParents) const;

    /**
     *  @brief Erase hits not associated with an MCParticle PDG incompatible with track/shower/all
     *
     *  @param[in] hitParents Map of hits to the Cluster/MCParticle they belong to
     *  @param[in] valType    Enum for track/shower/all
     */
    std::map<const pandora::CaloHit *const, CaloHitParents> ApplyPDGCut(
        std::map<const pandora::CaloHit *const, CaloHitParents> &hitParents, const ValidationType &valType) const;

    /**
     *  @brief Update the branches of the TTree for this entry
     *
     *  @param[in] metrics      Metrics for clusters in a view
     *  @param[in] randIfx      Adjusted Rand index
     *  @param[in] branchPrefix Prefix for the branch name
     */
    void SetBranches(ClusterMetrics &metrics, float randIdx, std::string branchPrefix) const;

    int m_eventNumber;                           ///< To track the current event number
    std::string m_fileName;                      ///< The filename of the ROOT output file
    std::string m_treeName;                      ///< The name of the ROOT tree
    std::string m_caloHitListName;               ///< The name of the hit list containing all 2D hits
    std::vector<std::string> m_clusterListNames; ///< The names of the lists of 2D clusters to process
    int m_minMCHitsPerView; ///< Threshold on total main MCParticle hits in each view for consideration in metric calculations
};

} // namespace lar_content

#endif // LAR_EVENT_CLUSTER_VALIDATION_ALGORITHM_H
