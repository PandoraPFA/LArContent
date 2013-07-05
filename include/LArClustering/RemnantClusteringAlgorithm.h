/**
 *  @file   LArContent/include/Clustering/RemnantClusteringAlgorithm.h
 * 
 *  @brief  Header file for the remnant clustering algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_REMNANT_CLUSTERING_ALGORITHM_H
#define LAR_REMNANT_CLUSTERING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  RemnantClusteringAlgorithm class
 */
class RemnantClusteringAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

private:
    pandora::StatusCode Run();
    pandora::StatusCode FindClusters(const pandora::CaloHitList* theHitList);
    pandora::StatusCode BuildClusters(const pandora::CaloHitList* theHitList);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string    m_nonSeedClusterListName;    ///< The non seed cluster list name

    float          m_minPulseHeight;            ///< Minimum pulse height
    float          m_clusterRadius;             ///< Clustering window, radius
    float          m_clusterRadiusSquared;      ///< Clustering window, radius squared)
    int            m_clusterLayers;             ///< Clustering window, layers
    int            m_clusterMinSize;            ///< Minimum cluster size

private:

    bool IsAssociated( pandora::CaloHit* CaloHitI, pandora::CaloHit* CaloHitJ );
    void MakeAssociation( pandora::CaloHit* CaloHitI, pandora::CaloHit* CaloHitJ );  

    int GetClusterID( pandora::CaloHit* pCaloHit );
    void SetClusterID( pandora::CaloHit* pCaloHit, int clusterID );
    void ResetClusterID( pandora::CaloHit* pCaloHitI, pandora::CaloHit* pCaloHitJ );

    typedef std::map<pandora::CaloHit*,int> HitAssociationMap;
    HitAssociationMap fHitAssociations;

    typedef std::map<int,int> ClusterAssociationMap;
    ClusterAssociationMap fClusterAssociations;
    ClusterAssociationMap fClusterSizes;

    typedef std::map< int, pandora::Cluster* > ClusterMap;
    ClusterMap fClusterMap;

    int fClusterID;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *RemnantClusteringAlgorithm::Factory::CreateAlgorithm() const
{
    return new RemnantClusteringAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_REMNANT_CLUSTERING_ALGORITHM_H
