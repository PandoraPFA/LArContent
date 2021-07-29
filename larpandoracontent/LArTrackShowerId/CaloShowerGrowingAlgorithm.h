/**
 *  @file   larpandoracontent/LArTrackShowerId/CaloShowerGrowingAlgorithm.h
 *
 *  @brief  Header file for the calorimetric shower growing algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CALO_SHOWER_GROWING_ALGORITHM_H
#define LAR_CALO_SHOWER_GROWING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  CaloShowerGrowingAlgorithm class
 */
class CaloShowerGrowingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CaloShowerGrowingAlgorithm();

private:
    class Bounds
    {
    public:
        /**
         *  Constructor. It doesn't really matter if the vertices correspond to top-left etc, as long as the coordinates are specified in a
         *  way that makes them cycle in either a counter-clockwise or clockwise direction about the chain tl -> bl -> br -> tr -> tl
         *
         *  @param  tl top left vertex
         *  @param  bl bottom left vertex
         *  @param  br bottom right bertex
         *  @param  tr top right vertex
         *  @param  pAlgorithm the parent algorithm
         */
        Bounds(const pandora::CartesianVector &tl, const pandora::CartesianVector &bl, const pandora::CartesianVector &br, const pandora::CartesianVector &tr);

        /*  @brief  Determines if a point is contained within a the bounding region
         *
         *  @return true if the point is within the bounds, false otherwise
         */
        bool Contains(const pandora::CartesianVector &point) const;

        const pandora::CartesianVector m_a;     ///< The top-left corner of the bounding region
        const pandora::CartesianVector m_b;     ///< The bottom-left corner of the bounding region
        const pandora::CartesianVector m_c;     ///< The bottom-right corner of the bounding region
        const pandora::CartesianVector m_d;     ///< The top-right corner of the bounding region

    private:
        /*  @brief  Determines if a point is contained within a triangle using the Barycentric technique, as described in
         *          Real-Time Collision Detection, C. Ericson (Morgan Kaufmann, 2005).
         *
         *  @param  ab the vector describing the edge ab of the triangle abc
         *  @param  ac the vector describing the edge ac of the triangle abc
         *  @param  ap the vector describing the displacement of a point p from vertex a of the triangle abc
         *
         *  @return true if the point p is within the bounds of the triangle abc, false otherwise
         */
        bool Contains(const pandora::CartesianVector &ab, const pandora::CartesianVector &ac, const pandora::CartesianVector &ap) const;

        // ATTN: The edges ab and ac describe one triangle of the bounding region, ac and ad describe the other
        const pandora::CartesianVector m_ab;    ///< Vector from tl to bl
        const pandora::CartesianVector m_ac;    ///< Vector from tl to br
        const pandora::CartesianVector m_ad;    ///< Vector from tl to tr
    };

    typedef std::map<const pandora::Cluster *, pandora::ClusterVector> ClusterAssociationMap;

    pandora::StatusCode Run();

    /**
     *  @brief  Grow the showers in a cluster list
     *
     *  @param  clusterList the list of available seed candidates
     */
    void GrowShowers(const pandora::ClusterList &clusterList) const;

    /**
     *  @brief  Get the seed clusters
     *
     *  @param  clusterList the list of available seed candidates
     *  @param  seedClusterList the output list of seed clusters
     */
    void GetSeedClusters(const pandora::ClusterList &clusterList, pandora::ClusterList &seedClusterList) const;

    /**
     *  @brief  Get the bounding box for a seed cluster
     *
     *  @param  pSeed the seed whose bounds should be calculated
     *  
     *  @return the bounding region
     */
    Bounds GetSeedBounds(const pandora::Cluster *pSeed) const;

    /**
     *  @brief  Determine if two clusters are associated
     *
     *  @param  pClusterSeed the seed cluster
     *  @param  pCluster the cluster to consider for association
     *
     *  @return Whether or not the clusters are associated
     */
    bool AreClustersAssociated(const pandora::Cluster *const pClusterSeed, const pandora::Cluster *const pCluster) const;

    /**
     *  @brief  Get a figure of merit representing the consistency of the provided cluster association map
     *
     *  @param  clusterAssociationMap the cluster association map
     *
     *  @return the figure of merit
     */
    float GetFigureOfMerit(const ClusterAssociationMap &clusterAssociationMap) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector m_inputClusterListNames; ///< The names of the input cluster lists
    unsigned int m_minCaloHitsForSeed; ///< The minimum number of calo hits per seed cluster
    bool m_visualize; ///< Whether or not to visualize the algorithm steps
};

} // namespace lar_content

#endif // #ifndef LAR_CALO_SHOWER_GROWING_ALGORITHM_H
