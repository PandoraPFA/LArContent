/**
 *  @file   HitWidthClusterMergingAlgorithm.h
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

    void RemoveShortcutAssociations(const pandora::ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const;

    bool AreClustersAssociated(const LArHitWidthHelper::ClusterParameters &currentFitParameters, const LArHitWidthHelper::ClusterParameters &testFitParameters) const;

    bool IsExtremalCluster(const bool isForward, const pandora::Cluster *const pCurrentCluster,  const pandora::Cluster *const pTestCluster) const;

    pandora::CartesianVector GetClusterDirection(const LArHitWidthHelper::ClusterParameters &clusterFitParameters, const pandora::CartesianVector &fitReferencePoint) const;

    void GetWeightedGradient(const LArHitWidthHelper::ClusterParameters &clusterFitParameters, bool isTransverse, pandora::CartesianVector &direction, pandora::CartesianVector &zIntercept, float &chiSquared, const pandora::CartesianVector &fitReferencePoint) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    
    std::string m_clusterListName;
    float m_maxConstituentHitWidth;
    float m_hitWidthScalingFactor;
    float m_fittingWeight;
    float m_minClusterWeight;
    float m_maxXMergeDistance;         //Distance either side of point
    float m_maxZMergeDistance;         //Distance either side of point
    float m_minMergeCosOpeningAngle;
    float m_minDirectionDeviationCosAngle;
};

} //namespace lar_content

#endif //LAR_HIT_WIDTH_CLUSTER_MERGING_ALGORITHM_H
