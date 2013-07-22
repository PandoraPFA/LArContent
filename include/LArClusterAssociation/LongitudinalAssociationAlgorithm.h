/**
 *  @file   LArContent/include/LArClusterAssociation/LongitudinalAssociationAlgorithm.h
 * 
 *  @brief  Header file for the cluster merging algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_LONGITUDINAL_ASSOCIATION_ALGORITHM_H
#define LAR_LONGITUDINAL_ASSOCIATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArClusterAssociation/ClusterAssociationAlgorithm.h"

#include "Helpers/ClusterHelper.h"

namespace lar
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
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Populate cluster vector with subset of cluster list, containing clusters judged to be clean
     * 
     *  @param  pClusterList address of the cluster list
     *  @param  clusterVector to receive the populated cluster vector
     */
    void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;

    /**
     *  @brief  Populate the cluster association map
     * 
     *  @param  clusterVector the cluster vector
     *  @param  clusterAssociationMap to receive the populated cluster association map
     */
    void PopulateClusterAssociationMap(const pandora::ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const;

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
        const float hitSizeZ, const pandora::ClusterHelper::ClusterFitResult &innerFit, const pandora::ClusterHelper::ClusterFitResult &outerFit) const;

    /**
     *  @brief  Determine which of two clusters is extremal
     * 
     *  @param  isForward whether propagation direction is forward
     *  @param  pCurrentCluster current extremal cluster
     *  @param  pTestCluster potential extremal cluster
     * 
     *  @return boolean
     */
    bool IsExtremalCluster(const bool isForward, const pandora::Cluster *const pCurrentCluster, const pandora::Cluster *const pTestCluster) const;

};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *LongitudinalAssociationAlgorithm::Factory::CreateAlgorithm() const
{
    return new LongitudinalAssociationAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_LONGITUDINAL_ASSOCIATION_ALGORITHM_H
