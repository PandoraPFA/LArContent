/**
 *  @file   larpandoracontent/LArVertex/DirectionClusterSplittingAlgorithm.h
 * 
 *  @brief  Header file for the candidate vertex creation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_DIRECTION_CLUSTER_SPLITTING_ALGORITHM_H
#define LAR_DIRECTION_CLUSTER_SPLITTING_ALGORITHM_H 1

#include "larpandoracontent/LArControlFlow/TrackDirectionTool.h"

#include "Pandora/Algorithm.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  DirectionClusterSplittingAlgorithm::Algorithm class
 */
class DirectionClusterSplittingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    DirectionClusterSplittingAlgorithm();

private:
    pandora::StatusCode Run();

    void CountClusters();

    void SelectClusters(pandora::ClusterVector &clusterVectorU, pandora::ClusterVector &clusterVectorV, pandora::ClusterVector &clusterVectorW);

    const pandora::Cluster* GetTargetCluster(pandora::ClusterVector &clusterVectorW);

    std::pair<const pandora::CaloHit*, const pandora::Cluster*> RetrieveSplitCaloHitClusterPair(const pandora::Cluster* pCluster, TrackDirectionTool::DirectionFitObject &fitResult);

    std::pair<const pandora::CaloHit*, const pandora::Cluster*> FindMatchingCaloHitClusterPair(TrackDirectionTool::DirectionFitObject &fitResult, pandora::ClusterVector &clusterVector);

    pandora::StatusCode SplitCluster(std::pair<const pandora::CaloHit*, const pandora::Cluster*> hitClusterPair, pandora::HitType listHitType);

    void DivideCaloHits(const pandora::Cluster* pCluster, const pandora::CaloHit* pTargetCaloHit, pandora::CaloHitList &caloHitList1, pandora::CaloHitList &caloHitList2);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector   m_inputClusterListNames;            ///< The list of cluster list names

    unsigned int            m_minClusterCaloHits;               ///< The min number of hits in base cluster selection method
    float                   m_minClusterLengthSquared;          ///< The min length (squared) in base cluster selection method
    bool                    m_enableDirection;                  ///< Flag for running batch jobs: whether this algorithm should do anything

    TrackDirectionTool      *m_pTrackDirectionTool;
};

} // namespace lar_content

#endif // #ifndef LAR_DIRECTION_CLUSTER_SPLITTING_ALGORITHM_H
