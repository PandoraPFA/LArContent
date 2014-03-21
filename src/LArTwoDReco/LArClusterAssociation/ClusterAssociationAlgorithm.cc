/**
 *  @file   LArContent/src/LArTwoDReco/LArClusterAssociation/ClusterAssociationAlgorithm.cc
 * 
 *  @brief  Implementation of the cluster association algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArTwoDReco/LArClusterAssociation/ClusterAssociationAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode ClusterAssociationAlgorithm::Run()
{
    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

    m_mergeMade = true;
    ClusterAssociationMap clusterAssociationMap;

    while (m_mergeMade)
    {
        // Unambiguous propagation
        while (m_mergeMade)
        {
            m_mergeMade = false;

            ClusterVector clusterVector;
            this->GetListOfCleanClusters(pClusterList, clusterVector);

            clusterAssociationMap.clear();
            this->PopulateClusterAssociationMap(clusterVector, clusterAssociationMap);

            for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
            {
                if (pClusterList->end() == pClusterList->find(*iter))
                    continue;

                this->UnambiguousPropagation(*iter, true,  clusterAssociationMap);
                this->UnambiguousPropagation(*iter, false, clusterAssociationMap);
            }
        }

        if (!m_resolveAmbiguousAssociations)
            continue;

        ClusterVector clusterVector;
        this->GetListOfCleanClusters(pClusterList, clusterVector);
        std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNOccupiedLayers);


        // Propagation with ambiguities
        for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
        {
            Cluster *pCluster = *iter;
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

void ClusterAssociationAlgorithm::UnambiguousPropagation(Cluster *pCluster, const bool isForward, ClusterAssociationMap &clusterAssociationMap) const
{
    Cluster *pClusterToEnlarge = pCluster;
    ClusterAssociationMap::iterator iterEnlarge = clusterAssociationMap.find(pClusterToEnlarge);

    if (clusterAssociationMap.end() == iterEnlarge)
        return;

    ClusterList &clusterListEnlarge(isForward ? iterEnlarge->second.m_forwardAssociations : iterEnlarge->second.m_backwardAssociations);

    if (clusterListEnlarge.size() != 1)
        return;

    Cluster *pClusterToDelete = *(clusterListEnlarge.begin());
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

void ClusterAssociationAlgorithm::AmbiguousPropagation(Cluster *pCluster, const bool isForward, ClusterAssociationMap &clusterAssociationMap) const
{
    ClusterAssociationMap::iterator iter = clusterAssociationMap.find(pCluster);

    if (clusterAssociationMap.end() == iter)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    Cluster *pExtremalCluster = pCluster;

    ClusterList firstClusterList;
    this->NavigateAlongAssociations(clusterAssociationMap, pCluster, isForward, pExtremalCluster, firstClusterList);
// ClusterList tempList1; tempList1.insert(pCluster);
// ClusterList tempList2; tempList2.insert(pExtremalCluster);
// PANDORA_MONITORING_API(SetEveDisplayParameters(false, false, -1, 1));
// PANDORA_MONITORING_API(VisualizeClusters(&firstClusterList, "ForwardProp", GREEN));
// PANDORA_MONITORING_API(VisualizeClusters(&tempList1, "InitialCluster", RED));
// PANDORA_MONITORING_API(VisualizeClusters(&tempList2, "ExtremalCluster", BLUE));
// PANDORA_MONITORING_API(ViewEvent());

    ClusterList secondClusterList;
    this->NavigateAlongAssociations(clusterAssociationMap, pExtremalCluster, !isForward, pExtremalCluster, secondClusterList);
// ClusterList tempList3; tempList3.insert(pExtremalCluster);
// PANDORA_MONITORING_API(SetEveDisplayParameters(false, false, -1, 1));
// PANDORA_MONITORING_API(VisualizeClusters(&secondClusterList, "BackwardProp", GREEN));
// PANDORA_MONITORING_API(VisualizeClusters(&tempList3, "ExtremalCluster", BLUE));
// PANDORA_MONITORING_API(ViewEvent());

    ClusterList daughterClusterList;

    if ( pCluster==pExtremalCluster )
    {
        for (ClusterList::const_iterator iter = firstClusterList.begin(), iterEnd = firstClusterList.end(); iter != iterEnd; ++iter)
        {
            if ((secondClusterList.end() != secondClusterList.find(*iter)) && (pCluster != (*iter)))
                daughterClusterList.insert(*iter);
        }
    }

// if ( daughterClusterList.empty() == false )
// {
// ClusterList tempList; tempList.insert(pCluster);
// PANDORA_MONITORING_API(SetEveDisplayParameters(false, false, -1, 1));
// PANDORA_MONITORING_API(VisualizeClusters(&tempList, "PrimaryCluster", GREEN));
// PANDORA_MONITORING_API(VisualizeClusters(&daughterClusterList, "Daughters", BLUE));
// PANDORA_MONITORING_API(ViewEvent());
// }

    for (ClusterList::const_iterator iter = daughterClusterList.begin(), iterEnd = daughterClusterList.end(); iter != iterEnd; ++iter)
    {
        this->UpdateForAmbiguousMerge(pCluster, *iter, isForward, clusterAssociationMap);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pCluster, *iter));
        m_mergeMade = true;
    } 
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterAssociationAlgorithm::UpdateForUnambiguousMerge(Cluster *pClusterToEnlarge, Cluster *pClusterToDelete, const bool isForwardMerge,
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

void ClusterAssociationAlgorithm::UpdateForAmbiguousMerge(Cluster *pClusterToEnlarge, Cluster *pClusterToDelete, const bool isForwardMerge,
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

void ClusterAssociationAlgorithm::NavigateAlongAssociations(const ClusterAssociationMap &clusterAssociationMap, Cluster *pCluster,
    const bool isForward, Cluster *&pExtremalCluster, ClusterList &clusterList) const
{
    ClusterAssociationMap::const_iterator iterAssociation = clusterAssociationMap.find(pCluster);

    if (clusterAssociationMap.end() == iterAssociation)
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    clusterList.insert(pCluster);

    if (this->IsExtremalCluster(isForward, pExtremalCluster, pCluster))
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
    m_resolveAmbiguousAssociations = true;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ResolveAmbiguousAssociations", m_resolveAmbiguousAssociations));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
