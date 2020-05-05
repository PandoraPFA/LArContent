/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/HitWidthClusterMergingAlgorithm.h
 *
 *  @brief  Header file for the hit width cluster merging algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_HIT_WIDTH_CLUSTER_MERGING_ALGORITHM_H
#define LAR_HIT_WIDTH_CLUSTER_MERGING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/ClusterAssociationAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArHitWidthHelper.h"

namespace lar_content
{

/**
 *  @brief HitWidthClusterMergingAlgorithm class
 */
class HitWidthClusterMergingAlgorithm : public ClusterAssociationAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    HitWidthClusterMergingAlgorithm();

private:
    void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;
    void PopulateClusterAssociationMap(const pandora::ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const;
    bool IsExtremalCluster(const bool isForward, const pandora::Cluster *const pCurrentCluster,  const pandora::Cluster *const pTestCluster) const;

    /**
     *  @brief  Determine whether two clusters are associated
     *
     *  @param  currentClusterParameters parameters defining the current cluster
     *  @param  testClusterParameters parameters defining the test cluster
     *
     *  @return boolean whether the clusters are associated
     */
    bool AreClustersAssociated(const LArHitWidthHelper::ClusterParameters &currentClusterParameters, const LArHitWidthHelper::ClusterParameters &testClusterParameters) const;

    /**
     *  @brief  Determine the position of the constituent hit that lies closest to a specified position
     *
     *  @param  position the point to which the consituent hits will be compared
     *  @param  constituentHitVector the input vector of constituent hits
     *  @param  closestPoint the position of the closest constituent hit
     *
     */
    void FindClosestPointToPosition(const pandora::CartesianVector &position, const LArHitWidthHelper::ConstituentHitVector &constituentHitVector,
        pandora::CartesianVector &closestPoint) const;

    /**
     *  @brief  Determine the cluster direction at a reference point by performing a weighted least squared fit to the input consitutent hit positions
     *          The fit is performed in a rotated cartesian coordinate system defined by a fitting axis determined by eigen
     *          The function is composed of two main loops that calculate the:
     *          1. rL & rT weighted means
     *          2. gradient of the fit
     *          The fit is performed using a subset of points that are closest to the fitReferencePoint
     *
     *  @param  constituentHitVector the input vector of constituent hits
     *  @param  direction the fitted cluster direction
     *  @param  fitReferencePoint the hits closest to this point are included in the fit
     *  @param  fittingWeight the weight that is considered in the fit
     *
     */
    void GetClusterDirection(const LArHitWidthHelper::ConstituentHitVector &constituentHitVector, pandora::CartesianVector &direction,
        const pandora::CartesianVector &fitReferencePoint, const float fittingWeight) const;

    /**
     *  @brief  Obtain a vector of the minimum number of hits closest to a reference point that exceed a given weight
     *
     *  @param  constituentHitVector the input vector of constituent hits
     *  @param  fitReferencePoint the reference point
     *  @param  fittingWeight the specified cumulative hit weight
     *  @param  constituentHitSubsetVector the subset of constituent hits
     *
     */
    void GetConstituentHitSubsetVector(const LArHitWidthHelper::ConstituentHitVector &constituentHitVector, const pandora::CartesianVector &fitReferencePoint,
        const float fittingWeight, LArHitWidthHelper::ConstituentHitVector &constituentHitSubsetVector) const;

    /**
     *  @brief  Obtain the axes of the fitting frame
     *
     *  @param  constituentHitSubsetVector the input vector of constituent hits
     *  @param  axisDirection the fitting 'x-axis'
     *  @param  orthoDirection the fitting 'z-axis'
     *
     */
    void GetFittingAxes(const LArHitWidthHelper::ConstituentHitVector &constituentHitSubsetVector, pandora::CartesianVector &axisDirection,
        pandora::CartesianVector &orthoDirection) const;

    /**
     *  @brief  Translate from (x, y, z) coordinates to (rL, rT) coordinates
     *
     *  @param  axisDirection the fitting 'x-axis'
     *  @param  constituentHitPosition the (x, y, z) position of a constituent hit
     *  @param  rL the fitting 'x' coordinate
     *  @param  rT the fitting 'z' coordinate
     *
     */
    void GetFittingCoordinates(const pandora::CartesianVector &axisDirection, const pandora::CartesianVector &constituentHitPosition, float &rL, float &rT) const;

    /**
     *  @brief  Translate a gradient in the fitting coordinate frame to a direction vector in the detector frame
     *
     *  @param  axisDirection the fitting 'x-axis'
     *  @param  gradient the gradient dT/dL
     *  @param  globalDirection the direction vector in the detector frame
     *
     */
    void GetGlobalDirection(const pandora::CartesianVector &axisDirection, const float gradient, pandora::CartesianVector &globalDirection) const;

    /**
     *  @brief  Remove 'shortcut' associations from the cluster association map
     *
     *  @param  clusterVector the vector of clusters considered in the merging process
     *  @param  clusterAssociationMap the mapping of clusters to forward/backward associations
     *
     */
    void RemoveShortcutAssociations(const pandora::ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_maxConstituentHitWidth;           ///< The maximum hit width of a constituent hit of broken up hit, units cm
    float m_hitWidthScalingFactor;            ///< The scaling factor of the hit widths
    float m_fittingWeight;                    ///< The maximum hit weight considered in the least squared fit
    float m_minClusterWeight;                 ///< The threshold hit weight of the original, unscaled cluster to be considered in the merging process
    float m_maxXMergeDistance;                ///< The maximum x distance between merging points of associated clusters, units cm
    float m_maxZMergeDistance;                ///< The maximum z distance between merging points of associated clusters, units cm
    float m_minMergeCosOpeningAngle;          ///< The minimum cosine opening angle of the directions of associated clusters
    float m_minDirectionDeviationCosAngle;    ///< The minimum cosine opening angle of the direction of and associated cluster before and after merge
    float m_minClusterSparseness;             ///< The threshold sparseness of a cluster to be considered in the merging process

    // ATTN Dangling pointers emerge during cluster merging, here explicitly not dereferenced
    mutable LArHitWidthHelper::ClusterToParametersMap m_clusterToParametersMap;   ///< The map [cluster -> cluster parameters]
};

} //namespace lar_content

#endif //LAR_HIT_WIDTH_CLUSTER_MERGING_ALGORITHM_H
