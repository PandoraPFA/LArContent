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

        void SetMainTrackCluster(const pandora::Cluster *const pMainTrackCluster);
        
        bool IsEndUpstream() const;        


    private:
        const pandora::Cluster    *m_pMainTrackCluster;
        bool                       m_isEndUpstream;
    };




    
class CosmicRayEndpointCorrectionAlgorithm : public CosmicRayTrackRefinementBaseAlgorithm<ClusterEndpointAssociation>
{
public:
    
    CosmicRayEndpointCorrectionAlgorithm();
    
private:

    //pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);


    void FindBestClusterAssociation(const pandora::ClusterVector &clusterVector, const SlidingFitResultMapPair &slidingFitResultMapPair,
        ClusterAssociationVector &clusterAssociationVector);
    
    bool IsDeltaRay(const pandora::Cluster *const pCluster, const pandora::CartesianVector &clusterMergePoint, const pandora::CartesianVector &clusterMergeDirection,
        const bool isEndUpstream) const;   

    void CreateMainTrack(ClusterEndpointAssociation &clusterEndpointAssociation, const ClusterToCaloHitListMap &clusterToCaloHitListMap, const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector, SlidingFitResultMapPair &slidingFitResultMapPair) const;

    void UpdateAfterMainTrackModification(const pandora::Cluster *const pMainTrackCluster, ClusterEndpointAssociation &clusterEndpointAssociation, pandora::ClusterVector &clusterVector, SlidingFitResultMapPair &slidingFitResultMapPai) const;   

    
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

inline const pandora::Cluster *ClusterEndpointAssociation::GetMainTrackCluster() const
{
    return m_pMainTrackCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------    

inline bool ClusterEndpointAssociation::IsEndUpstream() const
{
    return m_isEndUpstream;
}

//------------------------------------------------------------------------------------------------------------------------------------------        

inline void ClusterEndpointAssociation::SetMainTrackCluster(const pandora::Cluster *const pMainTrackCluster)
{
    m_pMainTrackCluster = pMainTrackCluster;
}    

} // namespace lar_content

#endif // #ifndef LAR_COSMIC_RAY_ENDPOINT_CORRECTION_ALGORITHM_H
