/**
 *  @file   LArContent/include/LArTwoDReco/LArClusterAssociation/CrossGapsAssociationAlgorithm.h
 *
 *  @brief  Header file for the cross gaps association algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CROSS_GAPS_ASSOCIATION_ALGORITHM_H
#define LAR_CROSS_GAPS_ASSOCIATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArObjects/LArTwoDSlidingFitResult.h"

#include "LArTwoDReco/LArClusterAssociation/ClusterAssociationAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  CrossGapsAssociationAlgorithm class
 */
class CrossGapsAssociationAlgorithm : public ClusterAssociationAlgorithm
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
    CrossGapsAssociationAlgorithm();

private:
    void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;
    void PopulateClusterAssociationMap(const pandora::ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const;
    bool IsExtremalCluster(const bool isForward, const pandora::Cluster *const pCurrentCluster, const pandora::Cluster *const pTestCluster) const;

    /**
     *  @brief  Determine whether two clusters are associated
     *
     *  @param  innerFitResult two dimensional sliding fit result for the inner cluster
     *  @param  outerFitResult two dimensional sliding fit result for the outer cluster
     *
     *  @return boolean
     */
    bool AreClustersAssociated(const TwoDSlidingFitResult &innerFitResult, const TwoDSlidingFitResult &outerFitResult) const;

    /**
     *  @brief  Sample points along the extrapolation from a starting position to a target fit result to declare cluster association
     *
     *  @param  startPosition the start position
     *  @param  startDirection the start direction
     *  @param  targetFitResult the target fit result
     *
     *  @return boolean
     */
    bool IsAssociated(const pandora::CartesianVector &startPosition, const pandora::CartesianVector &startDirection,
	const TwoDSlidingFitResult &targetFitResult) const;

    /**
     *  @brief  Whether a sampling point lies near a target 2d sliding fit result
     *
     *  @param  samplingPoint the sampling point
     *  @param  targetFitResult the target fit result
     *
     *  @return boolean
     */
    bool IsNearCluster(const pandora::CartesianVector &samplingPoint, const TwoDSlidingFitResult &targetFitResult) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int    m_minClusterHits;               ///< The minimum allowed number of hits in a clean cluster
    unsigned int    m_minClusterLayers;             ///< The minimum allowed number of layers for a clean cluster
    unsigned int    m_slidingFitWindow;             ///< The layer window for the sliding linear fits
    unsigned int    m_maxSamplingPoints;            ///< The maximum number of extension sampling points considered per association check
    float           m_sampleStepSize;               ///< The sampling step size used in association checks, units cm
    unsigned int    m_maxUnmatchedSampleRun;        ///< The maximum run of unmatched (and non-gap) samples to consider before stopping
    float           m_maxOnClusterDistance;         ///< The maximum distance between a sampling point and sliding fit to target cluster
    unsigned int    m_minMatchedSamplingPoints;     ///< Minimum number of matched sampling points to declare association
    float           m_minMatchedSamplingFraction;   ///< Minimum ratio between matched sampling points and expectation to declare association
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CrossGapsAssociationAlgorithm::Factory::CreateAlgorithm() const
{
    return new CrossGapsAssociationAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_CROSS_GAPS_ASSOCIATION_ALGORITHM_H
