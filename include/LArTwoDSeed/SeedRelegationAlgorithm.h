/**
 *  @file   LArContent/include/LArTwoDSeed/SeedRelegationAlgorithm.h
 * 
 *  @brief  Header file for the seed relegation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_SEED_RELEGATION_ALGORITHM_H
#define LAR_SEED_RELEGATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  SeedRelegationAlgorithm class
 */
class SeedRelegationAlgorithm : public pandora::Algorithm
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
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Selection requirement on length of cluster
     * 
     *  @return boolean
     */
    bool IsSeedClusterLongEnough( const pandora::Cluster* const pCluster ) const;

    /**
     *  @brief  Check whether a cluster is associated with the reconstructed vertex
     * 
     *  @return boolean
     */
    bool IsSeedClusterConnectedToVertex( const pandora::Cluster* const pCluster ) const;

    /**
     *  @brief  Check whether a cluster is associated with another, larger, cluster
     * 
     *  @return boolean
     */
    bool IsSeedClusterConnectedToLargerCluster( const pandora::Cluster* const pClusterI, const pandora::Cluster* const pClusterJ );

    /**
     *  @brief  Check whether a cluster subtends a small angle to another cluster
     * 
     *  @return boolean
     */
    bool IsSeedClusterAtSmallAngleToLargerCluster( const pandora::CartesianVector positionI, const pandora::CartesianVector directionI, const pandora::CartesianVector positionJ, const pandora::CartesianVector directionJ );


    std::string         m_seedClusterListName;      ///< The seed cluster list name
    std::string         m_nonSeedClusterListName;   ///< The non seed cluster list name

    float               m_minClusterLength;
    float               m_minClusterDotProduct;
    float               m_minClusterSeparation;
    float               m_minClusterRadialSeparation;
    float               m_minClusterVertexSeparation;

    float               m_vertexInnerRadialDistance;
    float               m_vertexOuterRadialDistance;
    float               m_vertexImpactParameter;

    float               m_minClusterLengthSquared;
    float               m_minClusterSeparationSquared;
    float               m_minClusterRadialSeparationSquared;
    float               m_minClusterVertexSeparationSquared;
    float               m_vertexInnerRadialDistanceSquared;
    float               m_vertexOuterRadialDistanceSquared;

};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *SeedRelegationAlgorithm::Factory::CreateAlgorithm() const
{
    return new SeedRelegationAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_SEED_RELEGATION_ALGORITHM_H
