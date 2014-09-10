/**
 *  @file   LArContent/include/LArTwoDReco/LArClusterAssociation/LongitudinalAssociationAlgorithm.h
 * 
 *  @brief  Header file for the longitudinal association algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_LONGITUDINAL_ASSOCIATION_ALGORITHM_H
#define LAR_LONGITUDINAL_ASSOCIATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "Helpers/ClusterFitHelper.h"

#include "LArTwoDReco/LArClusterAssociation/ClusterAssociationAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  LongitudinalAssociationAlgorithm class
 */
class LongitudinalAssociationAlgorithm : public ClusterAssociationAlgorithm
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
     *  @param  hitSizeX hit size x
     *  @param  hitSizeZ hit size z
     *  @param  innerFit inner cluster fit result
     *  @param  outerFit outer cluster fit result
     * 
     *  @return whether the clusters are associated
     */
    bool AreClustersAssociated(const pandora::CartesianVector &innerClusterEnd, const pandora::CartesianVector &outerClusterStart, const float hitSizeX,
        const float hitSizeZ, const pandora::ClusterFitResult &innerFit, const pandora::ClusterFitResult &outerFit) const;

    unsigned int m_minClusterLayers;            ///< minimum allowed number of layers for a clean cluster
    unsigned int m_maxGapLayers;                ///< maximum allowed number of layers between associated clusters
    unsigned int m_fitLayers;                   ///< number of layers to fit at start and end of cluster
    float        m_maxGapDistanceSquared;       ///< maximum allowed distance (squared) between associated clusters
    float        m_minCosRelativeAngle;         ///< maximum allowed relative angle between associated clusters
    float        m_maxTransverseDisplacement;   ///< maximum allowed transverse displacement after extrapolation (normalised to cell size)
    float        m_maxLongitudinalDisplacement; ///< maximum allowed longitudinal displacement after extrapolation (normalised to cell size)
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *LongitudinalAssociationAlgorithm::Factory::CreateAlgorithm() const
{
    return new LongitudinalAssociationAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_LONGITUDINAL_ASSOCIATION_ALGORITHM_H
