/**
 *  @file   LArContent/include/LArTwoDReco/LArCosmicRay/CosmicRayShowerMatchingAlgorithm.h
 * 
 *  @brief  Header file for the cosmic ray shower matching algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_COSMIC_RAY_SHOWER_MATCHING_ALGORITHM_H
#define LAR_COSMIC_RAY_SHOWER_MATCHING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  CosmicRayShowerMatchingAlgorithm class
 */
class CosmicRayShowerMatchingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::map<const pandora::Cluster *const, const pandora::ParticleFlowObject *const> ctopmap_t;
    typedef std::multimap<const pandora::ParticleFlowObject *const, const pandora::Cluster *const> ptocmultimap_t;
    typedef ptocmultimap_t::iterator ptocIter;
    typedef ctopmap_t::iterator ctopIter;
    typedef std::map<const pandora::Cluster *const, const pandora::ParticleFlowObject *const>::value_type ctopValType;
    typedef std::multimap<const pandora::ParticleFlowObject *const, const pandora::Cluster *const>::value_type ptocValType;

    /**
     *  @brief  Add associated clusters to the cosmic ray PFOs
     * 
     *  @param  map of clusters to be associated to pfos
     *  @param  map of pfos pointing back to associated clusters
     */
    pandora::StatusCode CosmicRayShowerMatching(ctopmap_t &clusterToPfoMap, ptocmultimap_t &pfoAssociatedClusterMap) const;

    /**
     *  @brief  Sort pfos by number of constituent hits
     * 
     *  @param  pLhs address of first pfo
     *  @param  pRhs address of second pfo
     */
     static bool SortPfosByNHits(const pandora::ParticleFlowObject *const pLhs, const pandora::ParticleFlowObject *const pRhs);

    /**
     *  @brief  Top level steering of the cluster to PFO assocation
     *
     *  @param  map of clusters to be associated to pfos
     *  @param  map of pfos pointing back to associated clusters
     */
    pandora::StatusCode CosmicRay3DShowerMatching(ctopmap_t &clusterToPfoMap, ptocmultimap_t &pfoAssociatedClusterMap) const;

    /**
     *  @brief  Look at consistency of a UVW combination
     *
     *  @param  pointer to U view cluster
     *  @param  pointer to V view cluster
     *  @param  pointer to W view cluster
     *  @param  pseudo chi2 of consistency
     */
    pandora::StatusCode CompareClusterTriplet(const pandora::Cluster *const pClusterU, const pandora::Cluster *const pClusterV,
        const pandora::Cluster *const pClusterW, float &pseudoChi2) const;

    /**
     *  @brief  Find best PFO to associate a UVW triplet
     *
     *  @param  pointer to U view cluster
     *  @param  pointer to V view cluster
     *  @param  pointer to W view cluster
     *  @param  pointer to best PFO
     *  @param  distance measure to best PFO
     */
    pandora::StatusCode FindBestCosmicPFO(const pandora::Cluster *const pClusterU, const pandora::Cluster *const pClusterV,
      const pandora::Cluster *const pClusterW, pandora::ParticleFlowObject* &pBestPFO, float &distanceToBestPFO) const;

    /**
     *  @brief  Find best PFO to associate a UV or UW or VW pair
     *
     *  @param  pointer to first view cluster
     *  @param  pointer to second view cluster
     *  @param  pointer to best PFO
     *  @param  distance measure to best PFO
     */
    pandora::StatusCode FindBestCosmicPFO(const pandora::Cluster *const pClusterView1, const pandora::Cluster *const pClusterView2, pandora::ParticleFlowObject* &pBestPFO, float &distanceToBestPFO) const;

    /**
     *  @brief  Find best PFO to associate a single unpaired cluster
     *
     *  @param  pointer to first view cluster
     *  @param  pointer to best PFO
     *  @param  distance measure to best PFO
     */
    pandora::StatusCode FindBestCosmicPFO(const pandora::Cluster *const pClusterView1, pandora::ParticleFlowObject* &pBestPFO, float &distanceToBestPFO) const;

    /**
     *  @brief  Visulise the cluster matches
     *
     *  @param  map of clusters to be associated to pfos
     *  @param  map of pfos pointing back to associated clusters
     */
    pandora::StatusCode VisualiseMatches(ctopmap_t &clusterToPfoMap, ptocmultimap_t &pfoAssociatedClusterMap) const;

    /**
     *  @brief  Make association of cluster to pfo in internal maps
     *
     *  @param  map of clusters to be associated to pfos
     *  @param  map of pfos pointing back to associated clusters
     *  @param  pointer to cluster
     *  @param  pounter to pfo
     */
    pandora::StatusCode AssociateClusterWithPfo(ctopmap_t &clusterToPfoMap, ptocmultimap_t &pfoAssociatedClusterMap, const pandora::Cluster *const pC1, const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Make association of cluster pair to pfo in internal maps
     *
     *  @param  map of clusters to be associated to pfos
     *  @param  map of pfos pointing back to associated clusters
     *  @param  pointer to cluster in view 1
     *  @param  pointer to cluster in view 2
     *  @param  pounter to pfo
     */
    pandora::StatusCode AssociateClusterWithPfo(ctopmap_t &clusterToPfoMap, ptocmultimap_t &pfoAssociatedClusterMap, const pandora::Cluster *const pC1, const pandora::Cluster *const pC2, const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Make association of cluster triplet to pfo in internal maps
     *
     *  @param  map of clusters to be associated to pfos
     *  @param  map of pfos pointing back to associated clusters
     *  @param  pointer to cluster in U view
     *  @param  pointer to cluster in V view
     *  @param  pointer to cluster in W view
     *  @param  pounter to pfo
     */
    pandora::StatusCode AssociateClusterWithPfo(ctopmap_t &clusterToPfoMap, ptocmultimap_t &pfoAssociatedClusterMap, const pandora::Cluster *const pC1, const pandora::Cluster *const pC2, const pandora::Cluster *const pC3, const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Return indicative coordinates of cluster at a position x using range indicated
     *
     *  @param  pointer to cluster
     *  @param  x (time) value at which coordinate is returned
     *  @param  minimum x value used to evaluate coordinate
     *  @param  maximum x value used to evaluate coordinate
     *  @param  span in number of hits
     */
    float GetCoordinateAtX(const pandora::Cluster *const pCluster, const float x, const float xmin, const float xmax, const int span)const;

    /**
     *  @brief  Whether a cluster is associated with a MC neutrino
     *
     *  @param  pCluster address of the cluster
     * 
     *  @return boolean
     */
    bool IsNeutrinoCluster(const pandora::Cluster *const pCluster) const;

    std::string             m_inputPfoListName;           ///< The input pfo list name (all PFOs)
    pandora::StringVector   m_inputClusterListNamesU;     ///< The input cluster list names for the u view
    pandora::StringVector   m_inputClusterListNamesV;     ///< The input cluster list names for the v view
    pandora::StringVector   m_inputClusterListNamesW;     ///< The input cluster list names for the w view
    bool                    m_visualiseAllMatches;        ///< Visualise all matches
    bool                    m_visualiseBadMatches;        ///< Visualise bad matches
    bool                    m_debugPrintingOn;            ///< Debug printing
    float                   m_distanceFor1ViewMatching;   ///< Distance cut for single view matching
    float                   m_distanceFor2ViewMatching;   ///< Distance cut for two view matching
    float                   m_distanceFor3ViewMatching;   ///< Distance cut for three view matching
    float                   m_chi2For3ViewMatching;       ///< Pseudo chi2 cut for three view matching
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CosmicRayShowerMatchingAlgorithm::Factory::CreateAlgorithm() const
{
    return new CosmicRayShowerMatchingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_COSMIC_RAY_SHOWER_MATCHING_ALGORITHM_H
