/**
 *  @file   LArContent/include/LArContentFast/CrossedTrackSplittingAlgorithmFast.h
 *
 *  @brief  Header file for the crossed track splitting algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CROSSED_TRACK_SPLITTING_ALGORITHM_FAST_H
#define LAR_CROSSED_TRACK_SPLITTING_ALGORITHM_FAST_H 1

#include "Pandora/Algorithm.h"

#include "LArTwoDReco/LArClusterSplitting/TwoDSlidingFitSplittingAndSwitchingAlgorithm.h"

namespace lar_content_fast
{

template<typename, unsigned int> class KDTreeLinkerAlgo;
template<typename, unsigned int> class KDTreeNodeInfoT;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  CrossedTrackSplittingAlgorithm class
 */
class CrossedTrackSplittingAlgorithm : public lar_content::TwoDSlidingFitSplittingAndSwitchingAlgorithm
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

    /**
     *  @brief  Default constructor
     */
    CrossedTrackSplittingAlgorithm();

private:
    typedef KDTreeLinkerAlgo<const pandora::CaloHit*, 2> HitKDTree2D;
    typedef KDTreeNodeInfoT<const pandora::CaloHit*, 2> HitKDNode2D;
    typedef std::vector<HitKDNode2D> HitKDNode2DList;

    typedef std::unordered_map<const pandora::Cluster*, pandora::ClusterList> ClusterToClustersMap;
    typedef std::unordered_map<const pandora::CaloHit*, const pandora::Cluster*> HitToClusterMap;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    pandora::StatusCode PreparationStep(const pandora::ClusterVector &clusterVector);
    pandora::StatusCode TidyUpStep();
    pandora::StatusCode FindBestSplitPosition(const lar_content::TwoDSlidingFitResult &slidingFit1, const lar_content::TwoDSlidingFitResult &slidingFit2,
        pandora::CartesianVector &splitPosition, pandora::CartesianVector &direction1, pandora::CartesianVector &direction2) const;

    /**
     *  @brief Find average positions of pairs of hits within a maximum separation
     *
     *  @param pCluster1 the first cluster
     *  @param pCluster2 the second cluster
     *  @param candidateList the output list to receive the average positions
     */
    void FindCandidateSplitPositions(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2,
        pandora::CartesianPointList &candidateList) const;

    float                   m_maxClusterSeparation;             ///< maximum separation of two clusters
    float                   m_maxClusterSeparationSquared;      ///< maximum separation of two clusters (squared)
    float                   m_minCosRelativeAngle;              ///< maximum relative angle between tracks after un-crossing

    float                   m_searchRegion1D;                   ///< Search region, applied to each dimension, for look-up from kd-trees
    ClusterToClustersMap    m_nearbyClusters;                   ///< The nearby clusters map
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CrossedTrackSplittingAlgorithm::Factory::CreateAlgorithm() const
{
    return new CrossedTrackSplittingAlgorithm();
}

} // namespace lar_content_fast

#endif // #ifndef LAR_CROSSED_TRACK_SPLITTING_ALGORITHM_FAST_H
