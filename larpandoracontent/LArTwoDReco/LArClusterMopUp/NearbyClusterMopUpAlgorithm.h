/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterMopUp/NearbyClusterMopUpAlgorithm.h
 *
 *  @brief  Header file for the nearby cluster mop up algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_NEARBY_CLUSTER_MOP_UP_ALGORITHM_H
#define LAR_NEARBY_CLUSTER_MOP_UP_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArTwoDReco/LArClusterMopUp/ClusterMopUpBaseAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  NearbyClusterMopUpAlgorithm class
 */

class NearbyClusterMopUpAlgorithm : public ClusterMopUpBaseAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    NearbyClusterMopUpAlgorithm();

private:
    void ClusterMopUp(const pandora::ClusterList &pfoClusters, const pandora::ClusterList &remnantClusters) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int    m_minHitsInCluster;         ///< Minimum number of hits in order to consider a cluster
    float           m_vertexProximity;          ///< Distance between cluster inner/outer centroid and vtx to declare cluster vtx associated
    float           m_minClusterSeparation;     ///< Minimum distance between parent and daughter clusters to declare clusters associated
    float           m_touchingDistance;         ///< Threshold (small) distance below which parent and daughter clusters are declated touching
};

} // namespace lar_content

#endif // #ifndef LAR_NEARBY_CLUSTER_MOP_UP_ALGORITHM_H
