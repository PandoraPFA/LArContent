/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/ThreeDAssociationAlgorithm.h
 *
 *  @brief  Header file for the 3D association algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_THREED_ASSOCIATION_ALGORITHM_H
#define LAR_THREED_ASSOCIATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "Helpers/ClusterFitHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/ClusterAssociationAlgorithm.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"

#include <string_view>

namespace lar_content
{

/**
 *  @brief  ThreeDAssociationAlgorithm class
 */
class ThreeDAssociationAlgorithm : public ClusterAssociationAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ThreeDAssociationAlgorithm();

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;
    void PopulateClusterAssociationMap(const pandora::ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const;
    bool IsExtremalCluster(const bool isForward, const pandora::Cluster *const pCurrentCluster, const pandora::Cluster *const pTestCluster) const;

    /**
     *  @brief  Determine whether two clusters are associated
     *
     *  @param  pInnerCluster address of the inner cluster
     *  @param  pOuterCluster address of the outer cluster
     *
     *  @return whether the clusters are associated
     */
    bool AreClustersAssociated(const pandora::Cluster *const pInnerCluster, const pandora::Cluster *const pOuterCluster) const;

    /**
     *  @brief  Determine whether two clusters are associated
     *
     *  @param  innerClusterEnd inner cluster end position
     *  @param  outerClusterStart outer cluster start position
     *  @param  innerFit inner cluster fit result
     *  @param  outerFit outer cluster fit result
     *
     *  @return whether the clusters are associated
     */
    bool AreClustersAssociated(const pandora::CaloHitList &pInnerCluster, const pandora::CaloHitList &pOuterCluster,
        const pandora::CartesianVector &innerCentroid, const pandora::CartesianVector &outerCentroid,
        const pandora::CartesianVector &innerDirection, const pandora::CartesianVector &outerDirection) const;

    void PCAFit(const pandora::Cluster *pCluster) const;

    void PCAFit(const pandora::CaloHitList &mergedClusterCaloHit,
        pandora::CartesianVector &centroid, pandora::CartesianVector &direction ) const;
    
    // For debugging purposes
    // TODO: Remove
    void VisualizeClusters(const pandora::Cluster *const pInner, const pandora::Cluster *const pOuter, 
        const std::string_view str) const;

    unsigned int m_minClusterLayers;     ///< minimum allowed number of layers for a clean cluster
    unsigned int m_maxGapLayers;         ///< maximum allowed number of layers between associated clusters
    unsigned int m_fitLayers;            ///< number of layers to fit at start and end of cluster
    float m_maxGapDistanceSquared;       ///< maximum allowed distance (squared) between associated clusters
    float m_minCosRelativeAngle;         ///< maximum allowed relative angle between associated clusters
    float m_minClusterLength;            ///< minimum cluster length
    float m_maxTransverseDisplacement;   ///< maximum allowed transverse displacement after extrapolation (normalised to cell size)
    float m_maxLongitudinalDisplacement; ///< maximum allowed longitudinal displacement after extrapolation (normalised to cell size)
    float m_hitSizeZ;                    ///< estimated hit size in z (wire number) dimension, units cm
    float m_hitSizeX;                    ///< estimated hit size in x (drift time) dimension, units cm
    mutable pandora::HitType m_view;     ///< The view to which the hits under consideration belong

    // struct to cache PCA attributes
    struct PCAAttr {
        pandora::CartesianVector centroid;
        pandora::CartesianVector direction;
    };

    typedef std::unordered_map<const pandora::Cluster*, PCAAttr> ClusterPCAMap;
    mutable ClusterPCAMap clusterToPCAMap;

    typedef std::map<const pandora::Cluster*, pandora::CaloHitList> ClusterToHitListMap;
    mutable ClusterToHitListMap clusterToHitListMap;
};

} // namespace lar_content

#endif // #ifndef LAR_LONGITUDINAL_ASSOCIATION_ALGORITHM_H