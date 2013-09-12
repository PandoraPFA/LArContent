/**
 *  @file   LArContent/src/LArHelpers/LArThreeDHelper.cc
 * 
 *  @brief  Implementation of the lar three dimensional reconstruction helper class.
 * 
 *  $Log: $
 */

#include "Helpers/XmlHelper.h"

#include "Objects/Cluster.h"

#include "LArHelpers/LArThreeDHelper.h"

namespace lar
{

using namespace pandora;

LArThreeDHelper::ClusterToIdMap LArThreeDHelper::m_seedClusterToIdMap;
LArThreeDHelper::IdToClusterListMap LArThreeDHelper::m_idToSeedClusterListMap;
LArThreeDHelper::IdToClusterListMap LArThreeDHelper::m_idToNonSeedClusterListMap;
ClusterList LArThreeDHelper::m_loneClusterList;

//------------------------------------------------------------------------------------------------------------------------------------------

HitType LArThreeDHelper::GetClusterHitType(const Cluster *const pCluster)
{
    HitType hitType(CUSTOM);

    if (pCluster->ContainsHitType(VIEW_U))
    {
        if (CUSTOM == hitType) hitType = VIEW_U;
        else throw StatusCodeException(STATUS_CODE_FAILURE);
    }

    if (pCluster->ContainsHitType(VIEW_V))
    {
        if (CUSTOM == hitType) hitType = VIEW_V;
        else throw StatusCodeException(STATUS_CODE_FAILURE);
    }

    if (pCluster->ContainsHitType(VIEW_W))
    {
        if (CUSTOM == hitType) hitType = VIEW_W;
        else throw StatusCodeException(STATUS_CODE_FAILURE);
    }

    if (CUSTOM == hitType)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    return hitType;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArThreeDHelper::StoreClusterComponents(const ClusterList &seedComponents, const ClusterList &nonSeedComponents)
{
    const unsigned int clusterId(m_idToSeedClusterListMap.size());

    if (m_idToNonSeedClusterListMap.size() != clusterId)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    m_idToSeedClusterListMap[clusterId] = seedComponents;
    m_idToNonSeedClusterListMap[clusterId] = nonSeedComponents;

    for (ClusterList::const_iterator iter = seedComponents.begin(), iterEnd = seedComponents.end(); iter != iterEnd; ++iter)
        m_seedClusterToIdMap[*iter] = clusterId;
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int LArThreeDHelper::GetClusterIdFromSeed(const Cluster *const pSeedCluster)
{
    ClusterToIdMap::const_iterator iter = m_seedClusterToIdMap.find(pSeedCluster);

    if (m_seedClusterToIdMap.end() == iter)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArThreeDHelper::GetAllClusterComponents(const Cluster *const pSeedCluster, ClusterList &clusterList)
{
    LArThreeDHelper::GetAllSeedComponents(pSeedCluster, clusterList);
    LArThreeDHelper::GetAllNonSeedComponents(pSeedCluster, clusterList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArThreeDHelper::GetAllSeedComponents(const Cluster *const pSeedCluster, ClusterList &clusterList)
{
    const unsigned int clusterId(LArThreeDHelper::GetClusterIdFromSeed(pSeedCluster));
    IdToClusterListMap::const_iterator iter = m_idToSeedClusterListMap.find(clusterId);

    if (m_idToSeedClusterListMap.end() == iter)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    clusterList.insert(iter->second.begin(), iter->second.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArThreeDHelper::GetAllNonSeedComponents(const Cluster *const pSeedCluster, ClusterList &clusterList)
{
    const unsigned int clusterId(LArThreeDHelper::GetClusterIdFromSeed(pSeedCluster));
    IdToClusterListMap::const_iterator iter = m_idToNonSeedClusterListMap.find(clusterId);

    if (m_idToNonSeedClusterListMap.end() == iter)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    clusterList.insert(iter->second.begin(), iter->second.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArThreeDHelper::Reset()
{
    m_seedClusterToIdMap.clear();
    m_idToSeedClusterListMap.clear();
    m_idToNonSeedClusterListMap.clear();
    m_loneClusterList.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LArThreeDHelper::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar
