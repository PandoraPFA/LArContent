/**
 *  @file   EMTrackAlgorithm.h
 *
 *  @brief  Header file for the em track algorithm class
 *
 *  $Log: $
 */
#ifndef LAR_EM_TRACK_ALGORITHM_H
#define LAR_EM_TRACK_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

namespace lar_content
{


/**
 *  @brief  EMTrackAlgorithm class
 */
class EMTrackAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    EMTrackAlgorithm();

    typedef std::unordered_map<const pandora::Cluster*, const pandora::Cluster*> ClusterToAssociatedClusterMap;
    
 private:
    pandora::StatusCode Run();
    void SelectCleanClusters(const pandora::ClusterList *pClusterList, pandora::ClusterVector &clusterVector);
    void FillClusterToAssociatedClusterMap(const pandora::ClusterVector &clusterVector, ClusterToAssociatedClusterMap &clusterToAssociatedClusterMap);
    bool AreClustersAssociated(const TwoDSlidingFitResult &currentClusterFit, const TwoDSlidingFitResult &testClusterFit);
    void GetClusterMergingCoordinates(const TwoDSlidingFitResult &cluster1FitResult, pandora::CartesianVector &cluster1MergePoint, pandora::CartesianVector &cluster1Direction, const TwoDSlidingFitResult &cluster2FitResult, pandora::CartesianVector &cluster2MergePoint, pandora::CartesianVector &cluster2Direction);
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_minCaloHits; // minimum number of calo hits for cluster to be considered

};

} // namespace lar_content


#endif // #ifndef LAR_EM_TRACK_ALGORITHM_H
