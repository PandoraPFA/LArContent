/**
 *  @file   LArContent/src/LArTwoDReco/LArClusterAssociation/ClusterAssociationAlgorithm.cc
 * 
 *  @brief  Implementation of the cluster association algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArContent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArContent/LArTwoDReco/LArClusterAssociation/ClusterAssociationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ClusterAssociationAlgorithm::ClusterAssociationAlgorithm() :
    m_mergeMade(false),
    m_resolveAmbiguousAssociations(true)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterAssociationAlgorithm::Run()
{
    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

    ClusterVector clusterVector;
    this->GetListOfCleanClusters(pClusterList, clusterVector);

    ClusterAssociationMap clusterAssociationMap;
    this->PopulateClusterAssociationMap(clusterVector, clusterAssociationMap);

    m_mergeMade = true;

    while (m_mergeMade)
    {
        // Unambiguous propagation
        while (m_mergeMade)
        {
            m_mergeMade = false;

            for (const Cluster *const pCluster : clusterVector)
            {
                // ATTN The clusterVector may end up with dangling pointers; only protected by this check against managed cluster list
                if (pClusterList->end() == pClusterList->find(pCluster))
                    continue;

                this->UnambiguousPropagation(pCluster, true,  clusterAssociationMap);
                this->UnambiguousPropagation(pCluster, false, clusterAssociationMap);
            }
        }

        if (!m_resolveAmbiguousAssociations)
            continue;

        // Propagation with ambiguities
        for (const Cluster *const pCluster : clusterVector)
        {
            // ATTN The clusterVector may end up with dangling pointers; only protected by this check against up-to-date association list
            ClusterAssociationMap::const_iterator mapIter = clusterAssociationMap.find(pCluster);

            if (clusterAssociationMap.end() == mapIter)
                continue;

            if (mapIter->second.m_backwardAssociations.empty() && !mapIter->second.m_forwardAssociations.empty())
                this->AmbiguousPropagation(pCluster, true, clusterAssociationMap);

            if (mapIter->second.m_forwardAssociations.empty() && !mapIter->second.m_backwardAssociations.empty())
                this->AmbiguousPropagation(pCluster, false, clusterAssociationMap);
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterAssociationAlgorithm::UnambiguousPropagation(const Cluster *const pCluster, const bool isForward, ClusterAssociationMap &clusterAssociationMap) const
{
    const Cluster *const pClusterToEnlarge = pCluster;
    ClusterAssociationMap::iterator iterEnlarge = clusterAssociationMap.find(pClusterToEnlarge);

    if (clusterAssociationMap.end() == iterEnlarge)
        return;

    ClusterList &clusterListEnlarge(isForward ? iterEnlarge->second.m_forwardAssociations : iterEnlarge->second.m_backwardAssociations);

    if (clusterListEnlarge.size() != 1)
        return;

    const Cluster *const pClusterToDelete = *(clusterListEnlarge.begin());
    ClusterAssociationMap::iterator iterDelete = clusterAssociationMap.find(pClusterToDelete);

    if (clusterAssociationMap.end() == iterDelete)
        return;

    ClusterList &clusterListDelete(isForward ? iterDelete->second.m_backwardAssociations : iterDelete->second.m_forwardAssociations);

    if (clusterListDelete.size() != 1)
        return;

    this->UpdateForUnambiguousMerge(pClusterToEnlarge, pClusterToDelete, isForward, clusterAssociationMap);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pClusterToEnlarge, pClusterToDelete));
    m_mergeMade = true;

    this->UnambiguousPropagation(pClusterToEnlarge, isForward, clusterAssociationMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterAssociationAlgorithm::AmbiguousPropagation(const Cluster *const pCluster, const bool isForward, ClusterAssociationMap &clusterAssociationMap) const
{
    ClusterAssociationMap::iterator cIter = clusterAssociationMap.find(pCluster);

    if (clusterAssociationMap.end() == cIter)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const Cluster *pExtremalCluster = pCluster;

    ClusterList firstClusterList;
    this->NavigateAlongAssociations(clusterAssociationMap, pCluster, isForward, pExtremalCluster, firstClusterList);

    ClusterList secondClusterList;
    this->NavigateAlongAssociations(clusterAssociationMap, pExtremalCluster, !isForward, pExtremalCluster, secondClusterList);

    ClusterVector daughterClusterVector;

    if (pCluster == pExtremalCluster)
    {
        for (ClusterList::const_iterator fIter = firstClusterList.begin(), fIterEnd = firstClusterList.end(); fIter != fIterEnd; ++fIter)
        {
            if ((secondClusterList.end() != secondClusterList.find(*fIter)) && (pCluster != (*fIter)))
                daughterClusterVector.push_back(*fIter);
        }
    }

    std::sort(daughterClusterVector.begin(), daughterClusterVector.end(), LArClusterHelper::SortByNHits);

    for (ClusterVector::iterator dIter = daughterClusterVector.begin(), dIterEnd = daughterClusterVector.end(); dIter != dIterEnd; ++dIter)
    {
        this->UpdateForAmbiguousMerge(pCluster, *dIter, isForward, clusterAssociationMap);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pCluster, *dIter));
        m_mergeMade = true;
        *dIter = NULL;
    } 
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterAssociationAlgorithm::UpdateForUnambiguousMerge(const Cluster *const pClusterToEnlarge, const Cluster *const pClusterToDelete, const bool isForwardMerge,
    ClusterAssociationMap &clusterAssociationMap) const
{
    ClusterAssociationMap::iterator iterEnlarge = clusterAssociationMap.find(pClusterToEnlarge);
    ClusterAssociationMap::iterator iterDelete = clusterAssociationMap.find(pClusterToDelete);

    if ((clusterAssociationMap.end() == iterEnlarge) || (clusterAssociationMap.end() == iterDelete))
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    ClusterList &clusterListToMove(isForwardMerge ? iterDelete->second.m_forwardAssociations : iterDelete->second.m_backwardAssociations);
    ClusterList &clusterListToReplace(isForwardMerge ? iterEnlarge->second.m_forwardAssociations : iterEnlarge->second.m_backwardAssociations);
    clusterListToReplace = clusterListToMove;
    clusterAssociationMap.erase(iterDelete);

    for (ClusterAssociationMap::iterator iter = clusterAssociationMap.begin(), iterEnd = clusterAssociationMap.end(); iter != iterEnd; ++iter)
    {
        ClusterList &forwardClusters = iter->second.m_forwardAssociations;
        ClusterList &backwardClusters = iter->second.m_backwardAssociations;

        ClusterList::iterator forwardIter = forwardClusters.find(pClusterToDelete);
        ClusterList::iterator backwardIter = backwardClusters.find(pClusterToDelete);

        if (forwardClusters.end() != forwardIter)
        {
            forwardClusters.erase(pClusterToDelete);
            forwardClusters.insert(pClusterToEnlarge);
        }

        if (backwardClusters.end() != backwardIter)
        {
            backwardClusters.erase(pClusterToDelete);
            backwardClusters.insert(pClusterToEnlarge);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterAssociationAlgorithm::UpdateForAmbiguousMerge(const Cluster *const pClusterToEnlarge, const Cluster *const pClusterToDelete, const bool isForwardMerge,
    ClusterAssociationMap &clusterAssociationMap) const
{
    ClusterAssociationMap::iterator iterEnlarge = clusterAssociationMap.find(pClusterToEnlarge);
    ClusterAssociationMap::iterator iterDelete = clusterAssociationMap.find(pClusterToDelete);

    if ((clusterAssociationMap.end() == iterEnlarge) || (clusterAssociationMap.end() == iterDelete))
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    ClusterList &clusterListEnlarge(isForwardMerge ? iterEnlarge->second.m_forwardAssociations : iterEnlarge->second.m_backwardAssociations);
    ClusterList &clusterListDelete(isForwardMerge ? iterDelete->second.m_backwardAssociations : iterDelete->second.m_forwardAssociations);

    for (ClusterList::iterator iter = clusterListEnlarge.begin(); iter != clusterListEnlarge.end();)
    {
        if ((*iter) != pClusterToDelete)
        {
            ClusterAssociationMap::iterator iterAssociation = clusterAssociationMap.find(*iter);

            if (clusterAssociationMap.end() == iterAssociation)
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);

            ClusterList &associatedClusterList(isForwardMerge ? iterAssociation->second.m_backwardAssociations : iterAssociation->second.m_forwardAssociations);

            ClusterList::iterator enlargeIter = associatedClusterList.find(pClusterToEnlarge);

            if (associatedClusterList.end() == enlargeIter)
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);

            associatedClusterList.erase(enlargeIter);
            clusterListEnlarge.erase(iter++);
        }
        else
        {
            ++iter;
        }
    }

    for (ClusterList::iterator iter = clusterListDelete.begin(); iter != clusterListDelete.end();)
    {
        if ((*iter) != pClusterToEnlarge)
        {
            ClusterAssociationMap::iterator iterAssociation = clusterAssociationMap.find(*iter);

            if (clusterAssociationMap.end() == iterAssociation)
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);

            ClusterList &associatedClusterList(isForwardMerge ? iterAssociation->second.m_forwardAssociations : iterAssociation->second.m_backwardAssociations);

            ClusterList::iterator deleteIter = associatedClusterList.find(pClusterToDelete);

            if (associatedClusterList.end() == deleteIter)
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);

            associatedClusterList.erase(deleteIter);
            clusterListDelete.erase(iter++);
        }
        else
        {
            ++iter;
        }
    }

    return this->UpdateForUnambiguousMerge(pClusterToEnlarge, pClusterToDelete, isForwardMerge, clusterAssociationMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterAssociationAlgorithm::NavigateAlongAssociations(const ClusterAssociationMap &clusterAssociationMap, const Cluster *const pCluster,
    const bool isForward, const Cluster *&pExtremalCluster, ClusterList &clusterList) const
{
    ClusterAssociationMap::const_iterator iterAssociation = clusterAssociationMap.find(pCluster);

    if (clusterAssociationMap.end() == iterAssociation)
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    clusterList.insert(pCluster);

    if ((pCluster != pExtremalCluster) && this->IsExtremalCluster(isForward, pExtremalCluster, pCluster))
          pExtremalCluster = pCluster;

    const ClusterList &associatedClusterList(isForward ? iterAssociation->second.m_forwardAssociations : iterAssociation->second.m_backwardAssociations);

    for (ClusterList::const_iterator iter = associatedClusterList.begin(), iterEnd = associatedClusterList.end(); iter != iterEnd; ++iter)
    {
        this->NavigateAlongAssociations(clusterAssociationMap, *iter, isForward, pExtremalCluster, clusterList);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterAssociationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ResolveAmbiguousAssociations", m_resolveAmbiguousAssociations));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
