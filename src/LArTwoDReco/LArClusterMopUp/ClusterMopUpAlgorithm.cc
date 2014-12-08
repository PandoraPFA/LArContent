/**
 *  @file   LArContent/src/LArTwoDReco/LArClusterMopUp/ClusterMopUpAlgorithm.cc
 * 
 *  @brief  Implementation of the cluster mop up algorithm base class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArTwoDReco/LArClusterMopUp/ClusterMopUpAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ClusterMopUpAlgorithm::ClusterMopUpAlgorithm() :
    m_excludePfosContainingTracks(true)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterMopUpAlgorithm::Run()
{
    ClusterList pfoClusterListU, pfoClusterListV, pfoClusterListW;
    this->GetPfoClusterLists(pfoClusterListU, pfoClusterListV, pfoClusterListW);

    ClusterList remnantClusterListU, remnantClusterListV, remnantClusterListW;
    this->GetRemnantClusterLists(remnantClusterListU, remnantClusterListV, remnantClusterListW);

    ClusterToListNameMap clusterToListNameMap;
    this->GetClusterToListNameMap(clusterToListNameMap);

    this->ClusterMopUp(pfoClusterListU, remnantClusterListU, clusterToListNameMap);
    this->ClusterMopUp(pfoClusterListV, remnantClusterListV, clusterToListNameMap);
    this->ClusterMopUp(pfoClusterListW, remnantClusterListW, clusterToListNameMap);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterMopUpAlgorithm::GetPfoClusterLists(ClusterList &clusterListU, ClusterList &clusterListV, ClusterList &clusterListW) const
{
    const PfoList *pPfoList = NULL;
    if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, m_pfoListName, pPfoList))
    {
        std::cout << "ClusterMopUpAlgorithm: pfo list " << m_pfoListName << " unavailable." << std::endl;
        return;
    }

    for (PfoList::const_iterator pIter = pPfoList->begin(), pIterEnd = pPfoList->end(); pIter != pIterEnd; ++pIter)
    {
        bool containsFixedTrack(false);
        const ClusterList &clusterList((*pIter)->GetClusterList());

        for (ClusterList::const_iterator cIter = clusterList.begin(), cIterEnd = clusterList.end(); cIter != cIterEnd; ++cIter)
        {
            if (m_excludePfosContainingTracks && (*cIter)->IsFixedMuon())
            {
                containsFixedTrack = true;
                break;
            }
        }

        if (!containsFixedTrack)
            this->GetClusterLists(clusterList, false, clusterListU, clusterListV, clusterListW);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterMopUpAlgorithm::GetRemnantClusterLists(ClusterList &clusterListU, ClusterList &clusterListV, ClusterList &clusterListW) const
{
    for (StringVector::const_iterator sIter = m_remnantClusterListNames.begin(), sIterEnd = m_remnantClusterListNames.end(); sIter != sIterEnd; ++sIter)
    {
        const ClusterList *pClusterList(NULL);
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, *sIter, pClusterList))
            continue;

        this->GetClusterLists(*pClusterList, true, clusterListU, clusterListV, clusterListW);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterMopUpAlgorithm::GetClusterLists(const ClusterList &inputClusterList, const bool availabilityFlag, ClusterList &clusterListU,
    ClusterList &clusterListV, ClusterList &clusterListW) const
{
    for (ClusterList::const_iterator cIter = inputClusterList.begin(), cIterEnd = inputClusterList.end(); cIter != cIterEnd; ++cIter)
    {
        Cluster *pCluster(*cIter);

        if (availabilityFlag != pCluster->IsAvailable())
            continue;

        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));

        if ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType))
            continue;

        ClusterList &target((TPC_VIEW_U == hitType) ? clusterListU : (TPC_VIEW_V == hitType) ? clusterListV : clusterListW);
        target.insert(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterMopUpAlgorithm::MakeClusterMerges(const ClusterAssociationMap &clusterAssociationMap, const ClusterToListNameMap &clusterToListNameMap) const
{
    for (ClusterAssociationMap::const_iterator rIter = clusterAssociationMap.begin(), rIterEnd = clusterAssociationMap.end(); rIter != rIterEnd; ++rIter)
    {
        Cluster *pRemnantCluster(rIter->first);
        const AssociationDetails associationDetails(rIter->second);

        Cluster *pBestPfoCluster(NULL);
        float bestFigureOfMerit(-std::numeric_limits<float>::max());

        for (AssociationDetails::const_iterator pIter = associationDetails.begin(), pIterEnd = associationDetails.end(); pIter != pIterEnd; ++pIter)
        {
            Cluster *pPfoCluster(pIter->first);
            const float figureOfMerit(pIter->second);

            if (figureOfMerit > bestFigureOfMerit)
            {
                pBestPfoCluster = pPfoCluster;
                bestFigureOfMerit = figureOfMerit;
            }
        }

        if (NULL == pBestPfoCluster)
            continue;

        ClusterToListNameMap::const_iterator listIterP = clusterToListNameMap.find(pBestPfoCluster);
        ClusterToListNameMap::const_iterator listIterR = clusterToListNameMap.find(pRemnantCluster);

        if ((clusterToListNameMap.end() == listIterP) || (clusterToListNameMap.end() == listIterR))
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pBestPfoCluster, pRemnantCluster,
            listIterP->second, listIterR->second));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterMopUpAlgorithm::GetClusterToListNameMap(ClusterToListNameMap &clusterToListNameMap) const
{
    StringSet stringSet;
    stringSet.insert(m_remnantClusterListNames.begin(), m_remnantClusterListNames.end());
    stringSet.insert(m_additionalClusterListNames.begin(), m_additionalClusterListNames.end());

    for (StringSet::const_iterator sIter = stringSet.begin(), sIterEnd = stringSet.end(); sIter != sIterEnd; ++sIter)
    {
        const ClusterList *pClusterList(NULL);
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, *sIter, pClusterList))
            continue;

        for (ClusterList::const_iterator cIter = pClusterList->begin(), cIterEnd = pClusterList->end(); cIter != cIterEnd; ++cIter)
        {
            if (!clusterToListNameMap.insert(ClusterToListNameMap::value_type(*cIter, *sIter)).second)
                throw StatusCodeException(STATUS_CODE_FAILURE);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterMopUpAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "RemnantClusterListNames", m_remnantClusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "AdditionalClusterListNames", m_additionalClusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ExcludePfosContainingTracks", m_excludePfosContainingTracks));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
