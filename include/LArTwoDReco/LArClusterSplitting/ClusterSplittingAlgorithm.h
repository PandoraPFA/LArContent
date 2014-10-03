/**
 *  @file   LArContent/include/LArTwoDReco/LArClusterSplitting/ClusterSplittingAlgorithm.h
 *
 *  @brief  Header file for the cluster splitting algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CLUSTER_SPLITTING_ALGORITHM_H
#define LAR_CLUSTER_SPLITTING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include <list>

namespace lar_content
{

/**
 *  @brief  ClusterSplittingAlgorithm class
 */
class ClusterSplittingAlgorithm : public pandora::Algorithm
{
protected:
    virtual pandora::StatusCode Run();
    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Divide calo hits in a cluster into two lists, each associated with a separate fragment cluster
     *
     *  @param  pCluster address of the cluster
     *  @param  firstCaloHitList the hits in the first fragment
     *  @param  secondCaloHitList the hits in the second fragment
     */
    virtual pandora::StatusCode DivideCaloHits(const pandora::Cluster *const pCluster, pandora::CaloHitList &firstCaloHitList,
        pandora::CaloHitList &secondCaloHitList) const = 0;

private:
    typedef std::list<pandora::Cluster*> ClusterSplittingList;

    /**
     *  @brief  Split cluster into two fragments
     *
     *  @param  pCluster address of the cluster
     *  @param  clusterSplittingList to receive the two cluster fragments
     */
    pandora::StatusCode SplitCluster(pandora::Cluster *const pCluster, ClusterSplittingList &clusterSplittingList) const;

    std::string     m_inputClusterList;         ///< The name of the input cluster list - if empty, use the current cluster list
};

} // namespace lar_content

#endif // #ifndef LAR_CLUSTER_SPLITTING_ALGORITHM_H
