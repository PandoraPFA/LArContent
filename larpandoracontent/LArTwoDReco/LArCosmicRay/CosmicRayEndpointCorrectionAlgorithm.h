/**
 *  @file   larpandoracontent/LArTwoDReco/LArCosmicRay/CosmicRayEndpointCorrectionAlgorithm.h
 *
 *  @brief  Header file for the cosmic ray endpoint correction class
 *
 *  $Log: $
 */
#ifndef LAR_COSMIC_RAY_ENDPOINT_CORRECTION_ALGORITHM_H
#define LAR_COSMIC_RAY_ENDPOINT_CORRECTION_ALGORITHM_H 1

//#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArTwoDReco/LArCosmicRay/CosmicRayTrackRefinementBaseAlgorithm.h"

namespace lar_content
{

class CosmicRayEndpointCorrectionAlgorithm : public CosmicRayTrackRefinementBaseAlgorithm
{
public:
    
    CosmicRayEndpointCorrectionAlgorithm();
    
private:

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    class ClusterEndpointAssociation : public ClusterAssociation
    {
    public:
        /**
         *  @brief  Default constructor
         */
        ClusterEndpointAssociation(const pandora::CartesianVector &upstreamMergePoint, const pandora::CartesianVector &upstreamMergeDirection,
            const pandora::CartesianVector &downstreamMergePoint, const pandora::CartesianVector &downstreamMergeDirection, const pandora::Cluster *pMainTrackCluster, const bool isEndUpstream);

        /**
         *  @brief  Returns the upstream cluster address
         *
         *  @return  Cluster the address of the upstream cluster
         */
        const pandora::Cluster *GetMainTrackCluster() const;

        /**
         *  @brief  Returns the upstream cluster address
         *
         *  @return  Cluster the address of the upstream cluster
         */
        const bool IsEndUpstream() const;        

        //void SetMainTrackCluster(const pandora::Cluster *pMainTrackCluster);

    private:
        const pandora::Cluster    *m_pMainTrackCluster;
        bool                       m_isEndUpstream;
    };


    /**
     *  @brief  Select clusters to be considered in algorithm
     *
     *  @param  pClusterList the input cluster list
     *  @param  clusterVector the output cluster vector
     */
    void SelectCleanClusters(const pandora::ClusterList *pClusterList, pandora::ClusterVector &clusterVector) const;

    void FindBestClusterAssociation(const pandora::ClusterVector &clusterVector, const SlidingFitResultMapPair &slidingFitResultMapPair,
        ClusterAssociationVector &clusterAssociationVector);
    
    bool IsDeltaRay(const pandora::Cluster *const pCluster, const pandora::CartesianVector &clusterMergePoint, const pandora::CartesianVector &clusterMergeDirection,
        const bool isEndUpstream) const;   

    void CreateMainTrack(const ClusterAssociationCaloHitOwnershipMap &clusterAssociationCaloHitOwnershipMap, const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector, SlidingFitResultMapPair &slidingFitResultMapPair) const;

    void UpdateAfterMainModification(const pandora::Cluster *const pDeletedCluster, const pandora::Cluster *const pNewCluster, pandora::ClusterVector &clusterVector, SlidingFitResultMapPair &slidingFitResultMapPair, ClusterAssociationCaloHitOwnershipMap &clusterAssociationCaloHitOwnershipMap) const;
    
    int m_minCaloHits;
    float m_maxDistanceFromTPC;
    float m_curveThreshold;
    float m_minScaledZOffset;
    float m_thresholdAngleDeviation;
    float m_thresholdAngleDeviationBetweenLayers;
    int m_maxAnomalousPoints;
    float m_thresholdMaxAngleDeviation;
    
};

//------------------------------------------------------------------------------------------------------------------------------------------    

inline const pandora::Cluster *CosmicRayEndpointCorrectionAlgorithm::ClusterEndpointAssociation::GetMainTrackCluster() const
{
    return m_pMainTrackCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------    

inline bool CosmicRayEndpointCorrectionAlgorithm::ClusterEndpointAssociation::IsEndUpstream() const
{
    return m_isEndUpstream;
}
        
    

} // namespace lar_content

#endif // #ifndef LAR_COSMIC_RAY_ENDPOINT_CORRECTION_ALGORITHM_H
