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
     *  @brief  Return the cluster direction at a reference point
     *
     *  @param  constituentHitVector the input vector of constituent hits
     *  @param  fitReferencePoint the hits closest to this point are included in the fit
     *
     *  @return CartesianVector the cluster direction at the reference point
     */
    pandora::CartesianVector GetClusterDirection(const LArHitWidthHelper::ConstituentHitVector &constituentHitVector, const pandora::CartesianVector &fitReferencePoint) const;

    /**
     *  @brief  Perform a weighted least squared fit to the input consitutent hit positions by minimising the longitudinal (as traditional) or transverse distance in the fit
     *          Function composed of three main loops to calculate the:
     *          1. x & z weighted means
     *          2. gradient and x/z intercept of the fit
     *          3. chi-squared of the fit
     *          The gradient & intercept are transformed such that the output is that of a z = mx + c fit
     *          The fit is performed using a subset of points that are closest to the fitReferencePoint
     *
     *  @param  unsortedConstituentHitVector the input vector of constituent hits to which the fit is applied
     *  @param  isLongitudinal whether to minimise the longitudinal distance (true) or transverse distance (false) in the fit
     *  @param  direction the fitted cluster direction
     *  @param  zIntercept the z intercept of the least squared fit
     *  @param  chiSquared the chi squared of the fit
     *  @param  fitReferencePoint the hits closest to this point are included in the fit
     */
    void GetWeightedGradient(const LArHitWidthHelper::ConstituentHitVector &unsortedConstituentHitVector, const bool isLongitudinal, pandora::CartesianVector &direction,
        pandora::CartesianVector &zIntercept, float &chiSquared, const pandora::CartesianVector &fitReferencePoint) const;


    void GetConstituentHitSubsetVector(const LArHitWidthHelper::ConstituentHitVector &constituentHitVector, const pandora::CartesianVector &fitReferencePoint,
        const float fittingWeight, LArHitWidthHelper::ConstituentHitVector &constituentHitSubsetVector) const;

    void GetFittingAxes(const LArHitWidthHelper::ConstituentHitVector &constituentHitSubsetVector, pandora::CartesianVector &axisDirection,
        pandora::CartesianVector &orthoDirection) const;

    void GetFittingCoordinates(const pandora::CartesianVector &axisDirection, const pandora::CartesianVector &constituentHitPosition, float &rL, float &rT) const;

    void GetGlobalDirection(const pandora::CartesianVector &axisDirection, const float gradient, pandora::CartesianVector &globalDirection) const;

    void GetWeightedGradient(const LArHitWidthHelper::ConstituentHitVector &constituentHitVector, pandora::CartesianVector &direction,
        const pandora::CartesianVector &fitReferencePoint, const float fittingWeight) const;

    void FindClosestPointToPosition(const pandora::CartesianVector &position, const LArHitWidthHelper::ConstituentHitVector &constituentHitVector,
        pandora::CartesianVector &closestPoint) const;
    
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
    float m_fittingWeight;                    ///< The maximum hit weight of considered in the least squared fit
    float m_minClusterWeight;                 ///< The threshold hit weight of the original, unscaled cluster to be considered in the merging process
    float m_maxXMergeDistance;                ///< The maximum x distance between merging points of associated clusters, units cm
    float m_maxZMergeDistance;                ///< The maximum z distance between merging points of associated clusters, units cm
    float m_minMergeCosOpeningAngle;          ///< The minimum cosine opening angle of the directions of associated clusters
    float m_minDirectionDeviationCosAngle;    ///< The minimum cosine opening angle of the direction of and associated cluster before and after merge
    float m_minClusterSparseness;
    bool  m_useOldDirectionMethod;
    bool  m_useClosestMergePoint;
    bool  m_doubleFittingWeight;


    // ATTN Dangling pointers emerge during cluster merging, here explicitly not dereferenced
    mutable LArHitWidthHelper::ClusterToParametersMap m_clusterToParametersMap;   ///< The map [cluster -> cluster parameters]
};

} //namespace lar_content

#endif //LAR_HIT_WIDTH_CLUSTER_MERGING_ALGORITHM_H
