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
     *  @brief  TwoDSlidingXZFitResult class
     */
    class TwoDSlidingXZFitResult
    {
    public:
        /**
         *  @brief  class LayerFitResult
         */
        class LayerFitResult
        {
        public:
            /**
             *  @brief  Constructor
             * 
             *  @param  z the z coordinate
             *  @param  fitX the fitted x coordinate
             *  @param  gradient the fitted gradient dx/dz
             *  @param  rms the rms of the fit residuals
             */
            LayerFitResult(const double z, const double fitX, const double gradient, const double rms);

            /**
             *  @brief  Get the z coordinate
             * 
             *  @return the z coordinate
             */
            double GetZ() const;

            /**
             *  @brief  Get the fitted x coordinate
             * 
             *  @return the fitted x coordinate
             */
            double GetFitX() const;

            /**
             *  @brief  Get the fitted gradient dx/dz
             * 
             *  @return the fitted gradient dx/dz
             */
            double GetGradient() const;

            /**
             *  @brief  Get the rms of the fit residuals
             * 
             *  @return the rms of the fit residuals
             */
            double GetRms() const;

            private:
            double                  m_z;                            ///< The z coordinate
            double                  m_fitX;                         ///< The fitted x coordinate
            double                  m_gradient;                     ///< The fitted gradient dx/dz
            double                  m_rms;                          ///< The rms of the fit residuals
        };

        typedef std::map<unsigned int, LayerFitResult> LayerFitResultMap;

        /**
         *  @brief  LayerFitContribution class
         */
        class LayerFitContribution
        {
        public:
            /**
             *  @brief  Constructor
             * 
             *  @param  pCaloHitList address of the calo hit list
             */
            LayerFitContribution(const pandora::CaloHitList *const pCaloHitList);

            /**
             *  @brief  Get the sum x
             * 
             *  @return the sum x
             */
            double GetSumX() const;

            /**
             *  @brief  Get the sum z
             * 
             *  @return the sum z
             */
            double GetSumZ() const;

            /**
             *  @brief  Get the sum x * x
             * 
             *  @return the sum x * x
             */
            double GetSumXX() const;

            /**
             *  @brief  Get the sum z * x
             * 
             *  @return the sum z * x
             */
            double GetSumZX() const;

            /**
             *  @brief  Get the sum z * z
             * 
             *  @return the sum z * z
             */
            double GetSumZZ() const;

            /**
             *  @brief  Get the number of points used
             * 
             *  @return the number of points used
             */
            unsigned int GetNPoints() const;

        private:
            double                  m_sumX;                ///< The sum x
            double                  m_sumZ;                ///< The sum z
            double                  m_sumXX;               ///< The sum x * x
            double                  m_sumZX;               ///< The sum z * x
            double                  m_sumZZ;               ///< The sum z * z
            unsigned int            m_nPoints;             ///< The number of points used
        };

        typedef std::map<unsigned int, LayerFitContribution> LayerFitContributionMap;

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
         *  @brief  Find the largest scatter in the cluster, if above a threshold value
         * 
         *  @param  largestScatterLayer to receive the layer corresponding to the largest scatter
         * 
         *  @return STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND
         */
        pandora::StatusCode FindLargestScatter(unsigned int &largestScatterLayer) const;

        /**
         *  @brief  Get the sliding fit width
         * 
         *  @return the sliding fit width
         */
        float GetSlidingFitWidth() const;

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
        LayerFitResultMap           m_layerFitResultMap;        ///< The layer fit result map
        LayerFitContributionMap     m_layerFitContributionMap;  ///< The layer fit contribution map

        friend class LArClusterHelper;
    };

    typedef TwoDSlidingXZFitResult::LayerFitResultMap LayerFitResultMap;
    typedef TwoDSlidingXZFitResult::LayerFitContributionMap LayerFitContributionMap;

    /**
     *  @brief  Perform two dimensional sliding x-z fit
     * 
     *  @param  pCluster address of the cluster
     *  @param  twoDSlidingXZFitResult to receive the fit result
     */
    static void LArTwoDSlidingXZFit(const pandora::Cluster *const pCluster, TwoDSlidingXZFitResult &twoDSlidingXZFitResult);

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

inline const pandora::Cluster *LArClusterHelper::TwoDSlidingXZFitResult::GetCluster() const
{
    return m_pCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int LArClusterHelper::TwoDSlidingXZFitResult::GetLayerFitHalfWindow() const
{
    return m_layerFitHalfWindow;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArClusterHelper::TwoDSlidingXZFitResult::LayerFitResultMap &LArClusterHelper::TwoDSlidingXZFitResult::GetLayerFitResultMap() const
{
    return m_layerFitResultMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArClusterHelper::TwoDSlidingXZFitResult::LayerFitContributionMap &LArClusterHelper::TwoDSlidingXZFitResult::GetLayerFitContributionMap() const
{
    return m_layerFitContributionMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArClusterHelper::TwoDSlidingXZFitResult::LayerFitResult::GetZ() const
{
    return m_z;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArClusterHelper::TwoDSlidingXZFitResult::LayerFitResult::GetFitX() const
{
    return m_fitX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArClusterHelper::TwoDSlidingXZFitResult::LayerFitResult::GetGradient() const
{
    return m_gradient;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArClusterHelper::TwoDSlidingXZFitResult::LayerFitResult::GetRms() const
{
    return m_rms;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArClusterHelper::TwoDSlidingXZFitResult::LayerFitContribution::GetSumX() const
{
    return m_sumX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArClusterHelper::TwoDSlidingXZFitResult::LayerFitContribution::GetSumZ() const
{
    return m_sumZ;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArClusterHelper::TwoDSlidingXZFitResult::LayerFitContribution::GetSumZX() const
{
    return m_sumZX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArClusterHelper::TwoDSlidingXZFitResult::LayerFitContribution::GetSumZZ() const
{
    return m_sumZZ;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double LArClusterHelper::TwoDSlidingXZFitResult::LayerFitContribution::GetSumXX() const
{
    return m_sumXX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int LArClusterHelper::TwoDSlidingXZFitResult::LayerFitContribution::GetNPoints() const
{
    return m_nPoints;
}

} // namespace lar

#endif // #ifndef LAR_CLUSTER_HELPER_H
