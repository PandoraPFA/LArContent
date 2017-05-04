/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/CrossedTrackSplittingAlgorithm.h
 *
 *  @brief  Header file for the crossed track splitting algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CROSSED_TRACK_SPLITTING_ALGORITHM_H
#define LAR_CROSSED_TRACK_SPLITTING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/TwoDSlidingFitSplittingAndSwitchingAlgorithm.h"

namespace lar_content
{

template<typename, unsigned int> class KDTreeLinkerAlgo;
template<typename, unsigned int> class KDTreeNodeInfoT;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  CrossedTrackSplittingAlgorithm class
 */
class CrossedTrackSplittingAlgorithm : public TwoDSlidingFitSplittingAndSwitchingAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CrossedTrackSplittingAlgorithm();

private:
    typedef KDTreeLinkerAlgo<const pandora::CaloHit*, 2> HitKDTree2D;
    typedef KDTreeNodeInfoT<const pandora::CaloHit*, 2> HitKDNode2D;
    typedef std::vector<HitKDNode2D> HitKDNode2DList;

    typedef std::unordered_map<const pandora::Cluster*, pandora::ClusterSet> ClusterToClustersMap;
    typedef std::unordered_map<const pandora::CaloHit*, const pandora::Cluster*> HitToClusterMap;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    pandora::StatusCode PreparationStep(const pandora::ClusterVector &clusterVector);
    pandora::StatusCode TidyUpStep();
    pandora::StatusCode FindBestSplitPosition(const TwoDSlidingFitResult &slidingFit1, const TwoDSlidingFitResult &slidingFit2,
        pandora::CartesianVector &splitPosition, pandora::CartesianVector &direction1, pandora::CartesianVector &direction2) const;

    /**
     *  @brief Find average positions of pairs of hits within a maximum separation
     *
     *  @param pCluster1 the first cluster
     *  @param pCluster2 the second cluster
     *  @param candidateVector to receive the average positions
     */
    void FindCandidateSplitPositions(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2,
        pandora::CartesianPointVector &candidateVector) const;

    float                   m_maxClusterSeparation;             ///< maximum separation of two clusters
    float                   m_maxClusterSeparationSquared;      ///< maximum separation of two clusters (squared)
    float                   m_minCosRelativeAngle;              ///< maximum relative angle between tracks after un-crossing

    float                   m_searchRegion1D;                   ///< Search region, applied to each dimension, for look-up from kd-trees
    ClusterToClustersMap    m_nearbyClusters;                   ///< The nearby clusters map
};

} // namespace lar_content

#endif // #ifndef LAR_CROSSED_TRACK_SPLITTING_ALGORITHM_H
