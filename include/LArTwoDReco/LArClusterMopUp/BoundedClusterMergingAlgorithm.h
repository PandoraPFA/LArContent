/**
 *  @file   LArContent/include/LArTwoDReco/LArClusterMopUp/BoundedClusterMergingAlgorithm.h
 * 
 *  @brief  Header file for the bounded cluster merging algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_BOUNDED_CLUSTER_MERGING_ALGORITHM_H
#define LAR_BOUNDED_CLUSTER_MERGING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArTwoDReco/LArClusterMopUp/ClusterMopUpAlgorithm.h"

namespace lar
{

/**
 *  @brief  BoundedClusterMergingAlgorithm class
 */

class BoundedClusterMergingAlgorithm : public ClusterMopUpAlgorithm
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
    /**
     *  @brief  XSampling class
     */
    class XSampling
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  fitResult the sliding fit result
         */
        XSampling(const TwoDSlidingFitResult &fitResult);

        float       m_minX;          ///< The min x value
        float       m_maxX;          ///< The max x value
        float       m_xPitch;        ///< The x sampling pitch to be used
    };

    /**
     *  @brief  ShowerEdge
     */
    class ShowerEdge
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  xCoordinate the x coordinate
         *  @param  highEdgeZ the shower high edge z coordinate
         *  @param  lowEdgeZ the shower low edge z coordinate
         */
        ShowerEdge(const float xCoordinate, const float highEdgeZ, const float lowEdgeZ);

        /**
         *  @param  Get the x coordinate
         * 
         *  @return the x coordinate
         */
        float GetXCoordinate() const;

        /**
         *  @param  Get the shower high edge z coordinate
         * 
         *  @return the shower high edge z coordinate
         */
        float GetHighEdgeZ() const;

        /**
         *  @param  Get the shower low edge z coordinate
         * 
         *  @return the shower low edge z coordinate
         */
        float GetLowEdgeZ() const;

    private:
        float       m_xCoordinate;   ///< The x coordinate
        float       m_highEdgeZ;     ///< The shower high edge z coordinate
        float       m_lowEdgeZ;      ///< The shower low edge z coordinate
    };

    typedef std::map<int, ShowerEdge> ShowerPositionMap;

    void ClusterMopUp(const pandora::ClusterList &pfoClusters, const pandora::ClusterList &remnantClusters, ClusterAssociationMap &clusterAssociationMap) const;

    /**
     *  @brief  Get the shower position map containing high and low edge z positions in bins of x
     * 
     *  @param  negativeEdgeFitResult the negative edge sliding fit result
     *  @param  positiveEdgeFitResult the positive edge sliding fit result
     *  @param  xSampling the x sampling details
     *  @param  showerPositionMap to receive the shower position map
     */
    void GetShowerPositionMap(const TwoDSlidingFitResult &negativeEdgeFitResult, const TwoDSlidingFitResult &positiveEdgeFitResult,
        const XSampling &xSampling, ShowerPositionMap &showerPositionMap) const;

    /**
     *  @brief  Get the fraction of hits in a cluster bounded by a specified shower position map
     * 
     *  @param  pCluster address of the cluster
     *  @param  xSampling the x sampling details
     *  @param  showerPositionMap the shower position map
     * 
     *  @return the fraction of bounded hits
     */
    float GetBoundedFraction(const pandora::Cluster *const pCluster, const XSampling &xSampling, const ShowerPositionMap &showerPositionMap) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int            m_slidingFitWindow;             ///< The layer window for the sliding linear fits
    float                   m_minBoundedFraction;           ///< The minimum cluster bounded fraction for merging
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *BoundedClusterMergingAlgorithm::Factory::CreateAlgorithm() const
{
    return new BoundedClusterMergingAlgorithm();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline BoundedClusterMergingAlgorithm::ShowerEdge::ShowerEdge(const float xCoordinate, const float edge1, const float edge2) :
    m_xCoordinate(xCoordinate),
    m_highEdgeZ(std::max(edge1, edge2)),
    m_lowEdgeZ(std::min(edge1, edge2))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float BoundedClusterMergingAlgorithm::ShowerEdge::GetXCoordinate() const
{
    return m_xCoordinate;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float BoundedClusterMergingAlgorithm::ShowerEdge::GetHighEdgeZ() const
{
    return m_highEdgeZ;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float BoundedClusterMergingAlgorithm::ShowerEdge::GetLowEdgeZ() const
{
    return m_lowEdgeZ;
}

} // namespace lar

#endif // #ifndef LAR_BOUNDED_CLUSTER_MERGING_ALGORITHM_H
