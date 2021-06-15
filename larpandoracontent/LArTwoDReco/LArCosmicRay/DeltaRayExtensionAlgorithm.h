/**
 *  @file   larpandoracontent/LArTwoDReco/LArCosmicRay/DeltaRayExtensionAlgorithm.h
 *
 *  @brief  Header file for the delta ray extension algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_DELTA_RAY_EXTENSION_ALGORITHM_H
#define LAR_DELTA_RAY_EXTENSION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/ClusterExtensionAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  DeltaRayExtensionAlgorithm class
 */
class DeltaRayExtensionAlgorithm : public ClusterExtensionAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    DeltaRayExtensionAlgorithm();

private:
    void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;
    void FillClusterAssociationMatrix(const pandora::ClusterVector &clusterVector, ClusterAssociationMatrix &clusterAssociationMatrix) const;
    void FillClusterMergeMap(const ClusterAssociationMatrix &clusterAssociationMatrix, ClusterMergeMap &clusterMergeMap) const;

    typedef std::unordered_map<const pandora::Cluster *, pandora::CartesianVector> ClusterToCoordinateMap;

    /**
     *  @brief  Reduce number of extremal coordinates calculations by caching results when they are first obtained
     *
     *  @param  pParentCluster the cluster
     *  @param  innerCoordinateMap the map from cluster to inner extremal coordinate
     *  @param  outerCoordinateMap the map from cluster to outer extremal coordinate
     *  @param  innerCoordinate to receive the inner coordinate
     *  @param  outerCoordinate to receive the outer coordinate
     */
    void GetExtremalCoordinatesFromCache(const pandora::Cluster *const pCluster, ClusterToCoordinateMap &innerCoordinateMap,
        ClusterToCoordinateMap &outerCoordinateMap, pandora::CartesianVector &innerCoordinate, pandora::CartesianVector &outerCoordinate) const;

    /**
     *  @brief  Form association between two clusters
     *
     *  @param  pParentCluster the parent cluster
     *  @param  pDaughterCluster the daughter cluster
     *  @param  innerCoordinateMap the map from cluster to inner extremal coordinate
     *  @param  outerCoordinateMap the map from cluster to outer extremal coordinate
     *  @param  clusterAssociationMatrix the matrix of cluster associations
     */
    void FillClusterAssociationMatrix(const pandora::Cluster *const pParentCluster, const pandora::Cluster *const pDaughterCluster,
        ClusterToCoordinateMap &innerCoordinateMap, ClusterToCoordinateMap &outerCoordinateMap, ClusterAssociationMatrix &clusterAssociationMatrix) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_minClusterLength; ///<
    float m_maxClusterLength; ///<

    float m_maxLongitudinalDisplacement; ///<
    float m_maxTransverseDisplacement;   ///<
};

} // namespace lar_content

#endif // #ifndef LAR_DELTA_RAY_EXTENSION_ALGORITHM_H
