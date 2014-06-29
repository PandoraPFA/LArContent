/**
 *  @file   LArContent/include/LArHelpers/LArClusterHelper.h
 *
 *  @brief  Header file for the cluster helper class.
 *
 *  $Log: $
 */
#ifndef LAR_CLUSTER_HELPER_H
#define LAR_CLUSTER_HELPER_H 1

#include "LArObjects/LArTwoDSlidingFitResult.h"

#include "Objects/Cluster.h"

namespace lar
{

/**
 *  @brief  TransverseDirection enum
 */
enum TransverseDirection
{
    POSITIVE_IN_X,
    NEGATIVE_IN_X,
    UNCHANGED_IN_X,
    UNKNOWN
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ShowerEdge enum
 */
enum ShowerEdge
{
    POSITIVE_SHOWER_EDGE,
    NEGATIVE_SHOWER_EDGE
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  LArClusterHelper class
 */
class LArClusterHelper
{
public:
    /**
     *  @brief  Get the hit type associated with a two dimensional cluster
     * 
     *  @param  pCluster the address of the cluster
     * 
     *  @return the cluster hit type
     */
    static pandora::HitType GetClusterHitType(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Perform two dimensional sliding fit, using a three dimensional fit to the cluster to define primary axis
     *
     *  @param  pCluster address of the cluster
     *  @param  layerFitHalfWindow the layer fit half window
     *  @param  twoDSlidingFitResult to receive the fit result
     */
    static void LArTwoDSlidingFit(const pandora::Cluster *const pCluster, const unsigned int layerFitHalfWindow, 
        TwoDSlidingFitResult &twoDSlidingFitResult);

    /**
     *  @brief  Perform two dimensional sliding fit, using z axis as primary axis, fitting x coordinates
     *
     *  @param  pCluster address of the cluster
     *  @param  layerFitHalfWindow the layer fit half window
     *  @param  twoDSlidingFitResult to receive the fit result
     */
    static void LArTwoDSlidingXZFit(const pandora::Cluster *const pCluster, const unsigned int layerFitHalfWindow, 
        TwoDSlidingFitResult &twoDSlidingFitResult);

    /**
     *  @brief  Perform two dimensional sliding fit, using the specified primary axis
     *
     *  @param  pCluster address of the cluster
     *  @param  layerFitHalfWindow the layer fit half window
     *  @param  axisIntercept the axis intercept position
     *  @param  axisDirection the axis direction vector
     *  @param  twoDSlidingFitResult to receive the fit result
     */
    static void LArTwoDSlidingFit(const pandora::Cluster *const pCluster, const unsigned int layerFitHalfWindow, 
        const pandora::CartesianVector &axisIntercept, const pandora::CartesianVector &axisDirection, 
        TwoDSlidingFitResult &twoDSlidingFitResult);

    /**
     *  @brief  Perform two dimensional sliding fit to shower edge, using specified primary axis
     *
     *  @param  pCluster address of the cluster
     *  @param  layerFitHalfWindow the layer fit half window
     *  @param  axisIntercept the axis intercept position
     *  @param  axisDirection the axis direction vector
     *  @param  showerEdge the shower edge
     *  @param  twoDSlidingFitResult to receive the fit result
     */
    static void LArTwoDShowerEdgeFit(const pandora::Cluster *const pCluster, const unsigned int layerFitHalfWindow, 
        const pandora::CartesianVector &axisIntercept, const pandora::CartesianVector &axisDirection, 
        const ShowerEdge showerEdge, TwoDSlidingFitResult &twoDSlidingFitResult);

    /**
     *  @brief  Whether fit results are multivalued in x
     *
     *  @return boolean
     */
    static bool IsMultivaluedInX(const TwoDSlidingFitResult &twoDSlidingFitResult);

    /**
     *  @brief  Get the sliding fit width
     *
     *  @return the sliding fit width
     */
    static float GetSlidingFitWidth(const TwoDSlidingFitResult &twoDSlidingFitResult);

    /**
     *  @brief  Measure width of cluster using multiple straight line fits
     *
     *  @param  pCluster address of the cluster
     *
     *  @return float
     */
    static float LArTrackWidth(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Get length squared of cluster
     *
     *  @param  pCluster address of the cluster
     *
     *  @return the length squared
     */
    static float GetLengthSquared(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Get length of cluster
     *
     *  @param  pCluster address of the cluster
     *
     *  @return the length
     */
    static float GetLength(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Get energy of cluster, based on length
     *
     *  @param  pCluster address of the cluster
     *
     *  @return the energy
     */
    static float GetEnergyFromLength(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Get number of layers spanned by cluster (1+Last-First)
     *
     *  @param  pCluster address of the cluster
     *
     *  @return the layer span
     */

    static unsigned int GetLayerSpan(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Fraction of occupied layers in cluster
     *
     *  @param  pCluster address of the cluster
     *
     *  @return float
     */
    static float GetLayerOccupancy(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Fraction of occupied layers in a pair of clusters
     *
     *  @param  pCluster1 address of the first cluster
     *  @param  pCluster2 address of the second cluster
     *
     *  @return float
     */
    static float GetLayerOccupancy(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2);

    /**
     *  @brief  Get closest distance between clusters in a pair of cluster vectors
     *
     *  @param  clusterVector1 the first cluster vector
     *  @param  clusterVector2 the second cluster vector
     *
     *  @return the closest distance
     */
    static float GetClosestDistance(const pandora::ClusterVector &clusterVector1, const pandora::ClusterVector &clusterVector2);

    /**
     *  @brief  Get closest distance between a specified cluster and vector of clusters
     *
     *  @param  pCluster address of the input cluster
     *  @param  clusterVector vector of input clusters
     *
     *  @return the closest distance
     */
    static float GetClosestDistance(const pandora::Cluster *const pCluster, const pandora::ClusterVector &clusterVector);

    /**
     *  @brief  Get closest distance between a pair of clusters
     *
     *  @param  pCluster1 address of the first cluster
     *  @param  pCluster2 address of the second cluster
     *
     *  @return the closest distance
     */
    static float GetClosestDistance(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2);

    /**
     *  @brief  Get closest distance between a specified position vector and the layer centroids of a specified cluster
     *
     *  @param  position the position vector
     *  @param  pCluster address of the cluster
     *
     *  @return the closest distance
     */
    static float GetClosestDistance(const pandora::CartesianVector &position, const pandora::Cluster *const pCluster);

    /**
     *  @brief  Get closest position on a cluster to a specified input position vector
     *
     *  @param  position the position vector
     *  @param  pCluster address of the cluster
     *
     *  @return the closest position
     */
    static pandora::CartesianVector GetClosestPosition(const pandora::CartesianVector &position, const pandora::Cluster *const pCluster);

    /**
     *  @brief  Get positions of the two most distant calo hits in a 2D cluster (ordered by Z)
     *
     *  @param  pCluster address of the cluster
     *  @param  the inner extremal position
     *  @param  the outer extremal position
     */
    static void GetExtremalCoordinatesXZ(const pandora::Cluster *const pCluster, pandora::CartesianVector &innerCoordinate, 
        pandora::CartesianVector &outerCoordinate);

    /**
     *  @brief  Get minimum and maximum X, Y and Z positions of the calo hits in a 2D cluster
     *
     *  @param  pCluster address of the cluster
     *  @param  the minimum positions (x,y,z)
     *  @param  the maximum positions (x,y,z)
     */
    static void GetClusterSpanXZ(const pandora::Cluster *const pCluster, pandora::CartesianVector &minimumCoordinate, 
        pandora::CartesianVector &maximumCoordinate);

    /**
     *  @brief  Get minimum and maximum x of the calo hits in a 2D cluster
     *
     *  @param  pCluster address of the cluster
     *  @param  the minimum position of x
     *  @param  the maximum position of x
     */
    static void GetClusterSpanX(const pandora::Cluster *const pCluster, float &xmin, float &xmax);

    /**
     *  @brief  Get upper and lower Z positions of the calo hits in a 2D cluster in range xmin to xmax
     *
     *  @param  pCluster address of the cluster
     *  @param  xmin for range in x
     *  @param  xmax for range in x
     *  @param  zmin the lower z for this range of x
     *  @param  zmax the upper z for this range in x
     */

    static void GetClusterSpanZ(const pandora::Cluster *const pCluster, const float xmin, const float xmax, float &zmin, float &zmax);

    /**
     *  @brief  Get average Z positions of the calo hits in a 2D cluster in range xmin to xmax
     *
     *  @param  pCluster address of the cluster
     *  @param  xmin for range in x
     *  @param  xmax for range in x
     */
    static float GetAverageZ(const pandora::Cluster *const pCluster, const float xmin, const float xmax);

    /**
     *  @brief  Sort clusters by inner layer (then use SortByNOccupiedLayers method in event of a tie)
     *
     *  @param  pLhs address of first cluster
     *  @param  pRhs address of second cluster
     */
    static bool SortByInnerLayer(const pandora::Cluster *const pLhs, const pandora::Cluster *const pRhs);

    /**
     *  @brief  Sort clusters by number of occupied layers, and by inner layer, then energy in event of a tie
     *
     *  @param  pLhs address of first cluster
     *  @param  pRhs address of second cluster
     */
    static bool SortByNOccupiedLayers(const pandora::Cluster *const pLhs, const pandora::Cluster *const pRhs);

    /**
     *  @brief  Sort clusters by number of hits and by layer span, then energy in event of a tie
     *
     *  @param  pLhs address of first cluster
     *  @param  pRhs address of second cluster
     */
    static bool SortByNHits(const pandora::Cluster *const pLhs, const pandora::Cluster *const pRhs);

    /**
     *  @brief  Read the vertex helper settings
     *
     *  @param  xmlHandle the relevant xml handle
     */
    static pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

private:
    /**
     *  @brief  Perform two dimensional sliding fit, using the information stored in the sliding fit result object
     *
     *  @param  twoDSlidingFitResult to receive the fit result
     */
    static void StoreSlidingFitResults(TwoDSlidingFitResult &twoDSlidingFitResult);

    /**
     *  @brief  Calculate a fit segment list for a given sliding fit result
     *
     *  @param  twoDSlidingFitResult the sliding fit result
     *  @param  fitSegmentList to receive the fit segment list
     */
    static void CalculateSlidingFitSegments(TwoDSlidingFitResult &twoDSlidingFitResult);

    static unsigned int             m_layerFitHalfWindow;           ///< The layer fit half window for sliding 2d x-z fits
    static float                    m_multiValuedTanThetaCut;       ///< Tan theta cut for finding sliding fits multivalued in x
    static float                    m_multiValuedStepFractionCut;   ///< Step fraction cut for finding sliding fits multivalued in x
    static float                    m_trackResidualQuantile;        ///< Track residual quantile, used for calculating sliding track width
};

} // namespace lar

#endif // #ifndef LAR_CLUSTER_HELPER_H
