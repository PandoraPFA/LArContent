/**
 *  @file   LArContent/include/LArHelpers/LArClusterHelper.h
 * 
 *  @brief  Header file for the cluster helper class.
 * 
 *  $Log: $
 */
#ifndef LAR_CLUSTER_HELPER_H
#define LAR_CLUSTER_HELPER_H 1

#include "Objects/Cluster.h"

namespace lar
{

/**
 *  @brief  ClusterQuality enum
 */
enum ClusterQuality
{
    METHOD_A = 0,  // NEED NAME
    METHOD_B,      // NEED NAME
    METHOD_C,      // NEED NAME
    METHOD_D       // NEED NAME
};

/**
 *  @brief  LArClusterHelper class
 */
class LArClusterHelper
{
public:
    /**
     *  @brief  TwoDSlidingFitResult class
     */
    class TwoDSlidingFitResult
    {
    public:
        /**
         *  @brief  Constructor
         */
        TwoDSlidingFitResult();

        /**
         *  @brief  class LayerFitResult
         */
        class LayerFitResult
        {
        public:
            /**
             *  @brief  Constructor
             * 
             *  @param  l the l coordinate
             *  @param  fitT the fitted t coordinate
             *  @param  gradient the fitted gradient dt/dl
             *  @param  rms the rms of the fit residuals
             */
            LayerFitResult(const double l, const double fitT, const double gradient, const double rms);

            /**
             *  @brief  Get the l coordinate
             * 
             *  @return the l coordinate
             */
            double GetL() const;

            /**
             *  @brief  Get the fitted t coordinate
             * 
             *  @return the fitted t coordinate
             */
            double GetFitT() const;

            /**
             *  @brief  Get the fitted gradient dt/dz
             * 
             *  @return the fitted gradient dt/dl
             */
            double GetGradient() const;

            /**
             *  @brief  Get the rms of the fit residuals
             * 
             *  @return the rms of the fit residuals
             */
            double GetRms() const;

            private:
            double                  m_l;                            ///< The l coordinate
            double                  m_fitT;                         ///< The fitted t coordinate
            double                  m_gradient;                     ///< The fitted gradient dt/dl
            double                  m_rms;                          ///< The rms of the fit residuals
        };

        typedef std::map<int, LayerFitResult> LayerFitResultMap;

        /**
         *  @brief  LayerFitContribution class
         */
        class LayerFitContribution
        {
        public:
            /**
             *  @brief  Default constructor
             */
            LayerFitContribution();

            /**
             *  @brief  Add point to layer fit
             * 
             *  @param  l the longitudinal coordinate
             *  @param  t the transverse coordinate
             */
            void AddPoint(const float l, const float t);

            /**
             *  @brief  Get the sum t
             * 
             *  @return the sum t
             */
            double GetSumT() const;

            /**
             *  @brief  Get the sum l
             * 
             *  @return the sum l
             */
            double GetSumL() const;

            /**
             *  @brief  Get the sum t * t
             * 
             *  @return the sum t * t
             */
            double GetSumTT() const;

            /**
             *  @brief  Get the sum l * t
             * 
             *  @return the sum l * t
             */
            double GetSumLT() const;

            /**
             *  @brief  Get the sum l * l
             * 
             *  @return the sum z * z
             */
            double GetSumLL() const;

            /**
             *  @brief  Get the number of points used
             * 
             *  @return the number of points used
             */
            unsigned int GetNPoints() const;

        private:
            double                  m_sumT;                ///< The sum t
            double                  m_sumL;                ///< The sum l
            double                  m_sumTT;               ///< The sum t * t
            double                  m_sumLT;               ///< The sum l * t
            double                  m_sumLL;               ///< The sum l * l
            unsigned int            m_nPoints;             ///< The number of points used
        };

        typedef std::map<int, LayerFitContribution> LayerFitContributionMap;

        /**
         *  @brief  Get the address of the cluster
         * 
         *  @return the address of the cluster
         */
        const pandora::Cluster *GetCluster() const;

        /**
         *  @brief  Get the layer fit half window
         * 
         *  @return the layer fit half window
         */
        unsigned int GetLayerFitHalfWindow() const;

        /**
         *  @brief  Get the axis intercept position
         * 
         *  @return the axis intercept position
         */
        const pandora::CartesianVector &GetAxisIntercept() const;

        /**
         *  @brief  Get the axis direction vector
         * 
         *  @return the axis direction vector
         */
        const pandora::CartesianVector &GetAxisDirection() const;

        /**
         *  @brief  Get sliding linear fit coordinates for a given cartesian vector
         * 
         *  @param  position the position cartesian vector
         *  @param  rL to receive the longitudinal coordinate
         *  @param  rT to receive the transverse coordinate
         */
        void GetLocalCoordinates(const pandora::CartesianVector &position, float &rL, float &rT) const;

        /**
         *  @brief  Get global coordinates for given sliding linear fit coordinates
         * 
         *  @param  rL the longitudinal coordinate
         *  @param  rT the transverse coordinate
         *  @param  position to receive the position cartesian vector
         */
        void GetGlobalCoordinates(const float rL, const float rT, pandora::CartesianVector &position) const;

        /**
         *  @brief  Get layer number for given sliding linear fit longitudinal coordinate
         * 
         *  @param  rL the longitudinal coordinate
         */
        int GetLayer(const float rL) const;

        /**
         *  @brief  Get sliding linear fit coordinates for a given x coordinate
         * 
         *  @param  x the x coordinate
         *  @param  rL to receive the longitudinal coordinate
         *  @param  rT to receive the transverse coordinate
         *  @param  layer to receive the layer
         */
        void GetLocalFitCoordinates(const float x, float &rL, float &rT, int &layer) const;

        /**
         *  @brief  Get global fit coordinates for a given x coordinate
         * 
         *  @param  x the x coordinate
         *  @param  position to receive the position cartesian vector
         */
        void GetGlobalFitCoordinates(const float x, pandora::CartesianVector &position) const;

        /**
         *  @brief  Get global position corresponding to the fit result in minimum fit layer
         * 
         *  @return the position
         */
        pandora::CartesianVector GetGlobalMinLayerPosition() const;

        /**
         *  @brief  Get global position corresponding to the fit result in maximum fit layer
         * 
         *  @return the position
         */
        pandora::CartesianVector GetGlobalMaxLayerPosition() const;

        /**
         *  @brief  Whether fit results are multivalued in x
         * 
         *  @return boolean
         */
        bool IsMultivaluedInX() const;

        /**
         *  @brief  Get the sliding fit width
         * 
         *  @return the sliding fit width
         */
        float GetSlidingFitWidth() const;

        /**
         *  @brief  Find the largest scatter in the cluster, if above a threshold value
         * 
         *  @param  largestScatterLayer to receive the layer corresponding to the largest scatter
         * 
         *  @return STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND
         */
        pandora::StatusCode FindLargestScatter(unsigned int &largestScatterLayer) const;

        /**
         *  @brief  Get the layer fit result map
         * 
         *  @return the layer fit result map
         */
        const LayerFitResultMap &GetLayerFitResultMap() const;

        /**
         *  @brief  Get the layer fit contribution map
         * 
         *  @return the layer fit contribution map
         */
        const LayerFitContributionMap &GetLayerFitContributionMap() const;

    private:
        const pandora::Cluster     *m_pCluster;                 ///< The address of the cluster
        unsigned int                m_layerFitHalfWindow;       ///< The layer fit half window
        pandora::CartesianVector    m_axisIntercept;            ///< The axis intercept position
        pandora::CartesianVector    m_axisDirection;            ///< The axis direction vector
        LayerFitResultMap           m_layerFitResultMap;        ///< The layer fit result map
        LayerFitContributionMap     m_layerFitContributionMap;  ///< The layer fit contribution map

        friend class LArClusterHelper;
    };

    /**
     *  @brief  Perform two dimensional sliding fit, using a three dimensional fit to the cluster to define primary axis
     * 
     *  @param  pCluster address of the cluster
     *  @param  layerFitHalfWindow the layer fit half window
     *  @param  twoDSlidingFitResult to receive the fit result
     */
    static void LArTwoDSlidingFit(const pandora::Cluster *const pCluster, const unsigned int layerFitHalfWindow, TwoDSlidingFitResult &twoDSlidingFitResult);

    /**
     *  @brief  Perform two dimensional sliding fit, using z axis as primary axis, fitting x coordinates
     * 
     *  @param  pCluster address of the cluster
     *  @param  layerFitHalfWindow the layer fit half window
     *  @param  twoDSlidingFitResult to receive the fit result
     */
    static void LArTwoDSlidingXZFit(const pandora::Cluster *const pCluster, const unsigned int layerFitHalfWindow, TwoDSlidingFitResult &twoDSlidingFitResult);

    /**
     *  @brief  Perform two dimensional sliding fit, using the specified primary axis
     * 
     *  @param  pCluster address of the cluster
     *  @param  layerFitHalfWindow the layer fit half window
     *  @param  axisIntercept the axis intercept position
     *  @param  axisDirection the axis direction vector
     *  @param  twoDSlidingFitResult to receive the fit result
     */
    static void LArTwoDSlidingFit(const pandora::Cluster *const pCluster, const unsigned int layerFitHalfWindow, const pandora::CartesianVector &axisIntercept,
        const pandora::CartesianVector &axisDirection, TwoDSlidingFitResult &twoDSlidingFitResult);

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
     *  @param  pCluster address of the first cluster
     * 
     *  @return the length squared
     */
    static float GetLengthSquared( const pandora::Cluster* const pCluster );

    /**
     *  @brief  Get length of cluster
     * 
     *  @param  pCluster address of the first cluster
     * 
     *  @return the length
     */
    static float GetLength( const pandora::Cluster* const pCluster ); 

    /**
     *  @brief  Get energy of cluster, based on length
     * 
     *  @param  pCluster address of the first cluster
     * 
     *  @return the energy
     */
    static float GetEnergyFromLength( const pandora::Cluster* const pCluster );

    /**
     *  @brief  Get number of layers spanned by cluster (1+Last-First)
     * 
     *  @param  pCluster address of the first cluster
     * 
     *  @return the layer span
     */

    static unsigned int GetLayerSpan( const pandora::Cluster* const pCluster );

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
     *  @brief  Get closest distance between the layer centroids of a pair of clusters
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
     * @brief Populate cluster vector with subset of cluster list, containing clusters judged to be clean
     *
     *  @param  method the type of selection 
     *  @param  pClusterList address of the cluster list
     *  @param  clusterVector to receive the populated cluster vector
     */
    static void GetListOfCleanClusters(const ClusterQuality method, const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector);

    /** 
     *  @brief Determine if cluster is clean
     *
     *  @param  method the type of selection
     *  @param  pCluster address of the cluster
     *
     *  @return boolean
     */
    static bool IsCleanCluster(const ClusterQuality method, const pandora::Cluster *const pCluster);

    /** 
     *  @brief Determine if cluster is clean (method A)
     *
     *  @param  pCluster address of the cluster
     *
     *  @return boolean
     */
    static bool IsCleanCluster_MethodA(const pandora::Cluster *const pCluster);

    /** 
     *  @brief Determine if cluster is clean (method B)
     *
     *  @param  pCluster address of the cluster 
     * 
     *  @return boolean 
     */
    static bool IsCleanCluster_MethodB(const pandora::Cluster *const pCluster);

    /** 
     *  @brief Determine if cluster is clean (method C)
     *
     *  @param  pCluster address of the cluster
     *  
     *  @return boolean
     */
    static bool IsCleanCluster_MethodC(const pandora::Cluster *const pCluster);

    /** 
     *  @brief Determine if cluster is clean (method D)
     *
     *  @param  pCluster address of the cluster 
     *
     *  @return boolean
     */
    static bool IsCleanCluster_MethodD(const pandora::Cluster *const pCluster);

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
    static unsigned int             m_layerFitHalfWindow;       ///< The layer fit half window for sliding 2d x-z fits
    static float                    m_trackFitMaxRms;           ///< Max RMS of track segment to be considered for kink splitting
    static float                    m_minCosScatteringAngle;    ///< Min kink angle at which to enable kink splitting
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Cluster *LArClusterHelper::TwoDSlidingFitResult::GetCluster() const
{
    return m_pCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int LArClusterHelper::TwoDSlidingFitResult::GetLayerFitHalfWindow() const
{
    return m_layerFitHalfWindow;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &LArClusterHelper::TwoDSlidingFitResult::GetAxisIntercept() const
{
    return m_axisIntercept;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &LArClusterHelper::TwoDSlidingFitResult::GetAxisDirection() const
{
    return m_axisDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArClusterHelper::TwoDSlidingFitResult::LayerFitResultMap &LArClusterHelper::TwoDSlidingFitResult::GetLayerFitResultMap() const
{
    return m_layerFitResultMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArClusterHelper::TwoDSlidingFitResult::LayerFitContributionMap &LArClusterHelper::TwoDSlidingFitResult::GetLayerFitContributionMap() const
{
    return m_layerFitContributionMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArClusterHelper::TwoDSlidingFitResult::LayerFitResult::GetL() const
{
    return m_l;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArClusterHelper::TwoDSlidingFitResult::LayerFitResult::GetFitT() const
{
    return m_fitT;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArClusterHelper::TwoDSlidingFitResult::LayerFitResult::GetGradient() const
{
    return m_gradient;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArClusterHelper::TwoDSlidingFitResult::LayerFitResult::GetRms() const
{
    return m_rms;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArClusterHelper::TwoDSlidingFitResult::LayerFitContribution::GetSumT() const
{
    return m_sumT;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArClusterHelper::TwoDSlidingFitResult::LayerFitContribution::GetSumL() const
{
    return m_sumL;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArClusterHelper::TwoDSlidingFitResult::LayerFitContribution::GetSumLT() const
{
    return m_sumLT;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArClusterHelper::TwoDSlidingFitResult::LayerFitContribution::GetSumLL() const
{
    return m_sumLL;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArClusterHelper::TwoDSlidingFitResult::LayerFitContribution::GetSumTT() const
{
    return m_sumTT;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int LArClusterHelper::TwoDSlidingFitResult::LayerFitContribution::GetNPoints() const
{
    return m_nPoints;
}

} // namespace lar

#endif // #ifndef LAR_CLUSTER_HELPER_H
