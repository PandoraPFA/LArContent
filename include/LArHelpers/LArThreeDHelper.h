/**
 *  @file   LArContent/include/LArHelpers/LArThreeDHelper.h
 * 
 *  @brief  Header file for the lar three dimensional reconstruction helper class.
 * 
 *  $Log: $
 */
#ifndef LAR_THREE_D_HELPER_H
#define LAR_THREE_D_HELPER_H 1

#include "Pandora/PandoraInternal.h"
#include "Pandora/PandoraInputTypes.h"
#include "Pandora/StatusCodes.h"

#include "Xml/tinyxml.h"

namespace lar
{

/**
 *  @brief  LArThreeDHelper class
 */
class LArThreeDHelper
{
public:
    /**
     *  @brief  Store components of a cluster that has undergone reclustering/fragmentation
     * 
     *  @param  seedComponents the list of seed component clusters
     *  @param  nonSeedComponents the list of non seed component clusters
     */
    static void StoreClusterComponents(const pandora::ClusterList &seedComponents, const pandora::ClusterList &nonSeedComponents);

    /**
     *  @brief  Remove all stored clusters from the helper class
     */
    static void RemoveAllStoredClusters();

    /**
     *  @brief  Get the number of stored clusters
     * 
     *  @return the number of stored clusters
     */
    static unsigned int GetNStoredClusters();

    /**
     *  @brief  Get the cluster id corresponding to a cluster originally containing a specified seed
     * 
     *  @param  pSeedCluster address of a seed component cluster
     * 
     *  @return the cluster id
     */
    static unsigned int GetClusterIdFromSeed(const pandora::Cluster *const pSeedCluster);

    /**
     *  @brief  Get a list of all cluster components corresponding to a cluster originally containing a specified seed
     * 
     *  @param  pSeedCluster address of a seed component cluster
     *  @param  clusterList to receive the component clusters
     */
    static void GetAllClusterComponents(const pandora::Cluster *const pSeedCluster, pandora::ClusterList &clusterList);

    /**
     *  @brief  Get a list of all seed components corresponding to a cluster originally containing a specified seed
     * 
     *  @param  pSeedCluster address of a seed component cluster
     *  @param  clusterList to receive the seed component clusters
     */
    static void GetAllSeedComponents(const pandora::Cluster *const pSeedCluster, pandora::ClusterList &clusterList);

    /**
     *  @brief  Get a list of all non seed components corresponding to a cluster originally containing a specified seed
     * 
     *  @param  pSeedCluster address of a seed component cluster
     *  @param  clusterList to receive the non seed component clusters
     */
    static void GetAllNonSeedComponents(const pandora::Cluster *const pSeedCluster, pandora::ClusterList &clusterList);

    /**
     *  @brief  Read the lar particle id settings
     * 
     *  @param  xmlHandle the relevant xml handle
     */
    static pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

private:
    typedef std::map<const pandora::Cluster*, unsigned int> ClusterToIdMap;
    typedef std::map<unsigned int, pandora::ClusterList> IdToClusterListMap;

    static ClusterToIdMap       m_seedClusterToIdMap;               ///< The seed cluster to id map
    static IdToClusterListMap   m_idToSeedClusterListMap;           ///< The id to seed cluster list map
    static IdToClusterListMap   m_idToNonSeedClusterListMap;        ///< The id to non seed cluster list map
};

} // namespace lar

#endif // #ifndef LAR_THREE_D_HELPER_H
