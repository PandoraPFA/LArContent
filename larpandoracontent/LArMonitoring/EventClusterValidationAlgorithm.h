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

        std::vector<const pandora::MCParticle *> m_pMainMCParticles;
        std::vector<int> m_nContributions;
        std::vector<int> m_nRecoHits;
        std::vector<float> m_purities;
        std::vector<float> m_completenesses;
        std::vector<float> m_fragmentationFractions;
    };

public:
    /**
    *  @brief  Default constructor
    */
    EventClusterValidationAlgorithm();

    virtual ~EventClusterValidationAlgorithm();

private:
    typedef std::map<pandora::Uid, float> MCSumHitWeightMap;

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Retrieve the metrics of every cluster in a view
     *
     *  @param viewClusters The clusters for this view
     *  @param viewCaloHits The calo hits for this view
     *  @param metricsMap The output metrics for the clusters in this view
     */
    void GetMetrics(
        const pandora::ClusterList &viewClusters,
        const pandora::CaloHitList &viewCaloHits,
        ClusterMetrics &metrics
    ) const;

    /**
     *  @brief  Collect the sum of hit weights for each MCParticle in each view
     *
     *  @param pCaloHitList The input calo hit list
     *  @param mcParticleSumWeight The output map of view to MCParticle to sum of all hit weights
     */
    void GetMCParticleViewSumHits(
        const pandora::CaloHitList *const pCaloHitList,
        std::map<pandora::HitType, MCSumHitWeightMap> &mcIDSumHitWeight
    ) const;

    int m_eventNumber;                           ///< To track the current event number
    std::string m_fileName;                      ///< The filename of the ROOT output file
    std::string m_treeName;                      ///< The name of the ROOT tree
    std::string m_caloHitListName;               ///< The name of the hit list containing all hits in all views
    std::vector<std::string> m_clusterListNames; ///< The names of the lists of clusters to process
    float m_minMCHitsPerView;                    ///< Threshold on MCParticle hits in each view
                                                 ///  for consideration in metric calculations
};

} // namespace lar_content

#endif // LAR_EVENT_CLUSTER_VALIDATION_ALGORITHM_H

