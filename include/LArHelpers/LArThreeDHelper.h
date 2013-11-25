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
     *  @brief  Get the hit type associated with a two dimensional cluster
     * 
     *  @param  pCluster the address of the cluster
     * 
     *  @return the cluster hit type
     */
    static pandora::HitType GetClusterHitType(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Store components of a cluster that has undergone reclustering/fragmentation
     * 
     *  @param  seedComponents the list of seed component clusters
     *  @param  nonSeedComponents the list of non seed component clusters
     */
    static void StoreClusterComponents(const pandora::ClusterList &seedComponents, const pandora::ClusterList &nonSeedComponents);

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
     *  @brief  Store a lone cluster
     * 
     *  @param  pCluster address of the lone cluster
     */
    static void StoreLoneCluster(pandora::Cluster *pCluster);

    /**
     *  @brief  Get the list of lone clusters stored in the three d helper
     * 
     *  @return the list of lone clusters
     */
    static const pandora::ClusterList &GetLoneClusterList();

    /**
     *  @brief  Remove all stored clusters from the helper class
     */
    static void Reset();

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
    static pandora::ClusterList m_loneClusterList;                  ///< The lone cluster list
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArThreeDHelper::StoreLoneCluster(pandora::Cluster *pCluster)
{
    (void) m_loneClusterList.insert(pCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::ClusterList &LArThreeDHelper::GetLoneClusterList()
{
    return m_loneClusterList;
}

} // namespace lar

#endif // #ifndef LAR_THREE_D_HELPER_H
