/**
 *  @file   LArContent/include/LArTwoDReco/LArClusterMopUp/ClusterMopUpAlgorithm.h
 * 
 *  @brief  Header file for the cluster mop up algorithm base class.
 * 
 *  $Log: $
 */
#ifndef LAR_CLUSTER_MOP_UP_ALGORITHM_H
#define LAR_CLUSTER_MOP_UP_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  ClusterMopUpAlgorithm class
 */

class ClusterMopUpAlgorithm : public pandora::Algorithm
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

protected:
    virtual pandora::StatusCode Run();

    /**
     *  @brief  Get the two dimensional clusters contained in the input pfo list, divided into three different views
     * 
     *  @param  clusterListU to receive the list of clusters in the u view
     *  @param  clusterListV to receive the list of clusters in the v view
     *  @param  clusterListW to receive the list of clusters in the w view
     */
    virtual void GetPfoClusterLists(pandora::ClusterList &clusterListU, pandora::ClusterList &clusterListV, pandora::ClusterList &clusterListW) const;

    /**
     *  @brief  Get the two dimensional clusters contained in the input remant cluster lists, divided into three different views
     * 
     *  @param  clusterListU to receive the list of clusters in the u view
     *  @param  clusterListV to receive the list of clusters in the v view
     *  @param  clusterListW to receive the list of clusters in the w view
     */
    virtual void GetRemnantClusterLists(pandora::ClusterList &clusterListU, pandora::ClusterList &clusterListV, pandora::ClusterList &clusterListW) const;

    /**
     *  @brief  Get the two dimensional clusters contained in an input cluster list, divided into three different views
     * 
     *  @param  inputClusterList the input cluster list
     *  @param  availabilityFlag only clusters with matching availability will be considered
     *  @param  clusterListU to receive the list of clusters in the u view
     *  @param  clusterListV to receive the list of clusters in the v view
     *  @param  clusterListW to receive the list of clusters in the w view
     */
    virtual void GetClusterLists(const pandora::ClusterList &inputClusterList, const bool availabilityFlag, pandora::ClusterList &clusterListU,
        pandora::ClusterList &clusterListV, pandora::ClusterList &clusterListW) const;

    typedef std::map<pandora::Cluster*, std::string> ClusterToListNameMap;

    /**
     *  @brief  Get the map from cluster to cluster list name
     * 
     *  @param  clusterToListNameMap to receive the cluster to list name map
     */
    virtual void GetClusterToListNameMap(ClusterToListNameMap &clusterToListNameMap) const;

    /**
     *  @brief  Cluster mop up for a single view. This function is responsible for instructing pandora to make cluster alterations
     * 
     *  @param  pfoClusters the list of pfo clusters
     *  @param  remnantClusters the list of remnant clusters
     *  @param  clusterToListNameMap the cluster list name map
     */
    virtual void ClusterMopUp(const pandora::ClusterList &pfoClusters, const pandora::ClusterList &remnantClusters, const ClusterToListNameMap &clusterToListNameMap) const = 0;

    typedef std::map<pandora::Cluster*, float> AssociationDetails;
    typedef std::map<pandora::Cluster*, AssociationDetails> ClusterAssociationMap;

    /**
     *  @brief  Make the cluster merges specified in the cluster association map, using list name information in the cluster list name map
     * 
     *  @param  clusterAssociationMap the cluster association map
     *  @param  clusterToListNameMap the cluster list name map
     */
    virtual void MakeClusterMerges(const ClusterAssociationMap &clusterAssociationMap, const ClusterToListNameMap &clusterToListNameMap) const;

    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string             m_pfoListName;                  ///< The pfo list name
    pandora::StringVector   m_remnantClusterListNames;      ///< The list of remnant cluster list names
    pandora::StringVector   m_additionalClusterListNames;   ///< The list of additional cluster list names, maybe specifying lists containing pfo clusters
};

} // namespace lar

#endif // #ifndef LAR_CLUSTER_MOP_UP_ALGORITHM_H
