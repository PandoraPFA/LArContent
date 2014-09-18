/**
 *  @file   LArContent/src/LArTwoDReco/LArClusterSplitting/TwoDSlidingFitConsolidationAlgorithm.cc
 *
 *  @brief  Implementation of the 2D sliding fit consolidation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"

#include "LArPlugins/LArTransformationPlugin.h"

#include "LArTwoDReco/LArClusterSplitting/TwoDSlidingFitConsolidationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

TwoDSlidingFitConsolidationAlgorithm::TwoDSlidingFitConsolidationAlgorithm() :
    m_minTrackLength(7.5f),
    m_maxClusterLength(15.f),
    m_halfWindowLayers(25)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitConsolidationAlgorithm::Run()
{
    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

    // Select tracks and showers for re-clustering (Note: samples not mutually exclusive)
    ClusterVector trackClusters, showerClusters;
    this->SortInputClusters(pClusterList, trackClusters, showerClusters);

    // Build sliding linear fits from track clusters
    TwoDSlidingFitResultList slidingFitResultList;
    this->BuildSlidingLinearFits(trackClusters, slidingFitResultList);

    // Recluster the hits
    ClusterToHitMap clustersToExpand, clustersToContract;
    this->GetReclusteredHits(slidingFitResultList, showerClusters, clustersToExpand, clustersToContract);

    // Consolidate and re-build clusters
    ClusterList unavailableClusters;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RemoveHitsFromClusters(clustersToContract, unavailableClusters));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->AddHitsToClusters(clustersToExpand, unavailableClusters));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RebuildClusters(clustersToContract, unavailableClusters));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitConsolidationAlgorithm::SortInputClusters(const ClusterList *const pClusterList, ClusterVector &trackClusters,
    ClusterVector &showerClusters) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        const float thisLengthSquared(LArClusterHelper::GetLengthSquared(pCluster));

        if (thisLengthSquared < m_maxClusterLength * m_maxClusterLength)
            showerClusters.push_back(pCluster);

        if (thisLengthSquared > m_minTrackLength * m_minTrackLength)
            trackClusters.push_back(pCluster);
    }

    std::sort(trackClusters.begin(), trackClusters.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitConsolidationAlgorithm::BuildSlidingLinearFits(const ClusterVector &trackClusters,
    TwoDSlidingFitResultList &slidingFitResultList) const
{
    const float slidingFitPitch(LArGeometryHelper::GetLArTransformationPlugin(this->GetPandora())->GetWireZPitch());

    for (ClusterVector::const_iterator iter = trackClusters.begin(), iterEnd = trackClusters.end(); iter != iterEnd; ++iter)
    {
        try
        {
            const TwoDSlidingFitResult slidingFitResult(*iter, m_halfWindowLayers, slidingFitPitch);
            slidingFitResultList.push_back(slidingFitResult);
        }
        catch (StatusCodeException &statusCodeException)
        {
            if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                throw statusCodeException;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitConsolidationAlgorithm::RemoveHitsFromClusters(const ClusterToHitMap &clustersToContract,
    ClusterList &unavailableClusters) const
{
    for (ClusterToHitMap::const_iterator iterI = clustersToContract.begin(), iterEndI = clustersToContract.end(); iterI != iterEndI; ++iterI)
    {
        Cluster* pCluster = const_cast<Cluster*>(iterI->first);
        const CaloHitList &caloHitListToRemove = iterI->second;

        if (caloHitListToRemove.empty())
            continue;

        if (unavailableClusters.count(pCluster))
            continue;

        CaloHitList caloHitList, caloHitListToKeep;
        pCluster->GetOrderedCaloHitList().GetCaloHitList(caloHitList);
        for (CaloHitList::const_iterator iterJ = caloHitList.begin(), iterEndJ = caloHitList.end(); iterJ != iterEndJ; ++iterJ)
        {
            CaloHit* pCaloHit = *iterJ;
            if (!caloHitListToRemove.count(pCaloHit))
                caloHitListToKeep.insert(pCaloHit);
        }

        if (caloHitListToKeep.empty())
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete<Cluster>(*this, pCluster));
            unavailableClusters.insert(pCluster);
            continue;
        }

        for (CaloHitList::const_iterator iterJ = caloHitListToRemove.begin(), iterEndJ = caloHitListToRemove.end(); iterJ != iterEndJ; ++iterJ)
        {
            CaloHit* pCaloHit = *iterJ;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromCluster(*this, pCluster, pCaloHit));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitConsolidationAlgorithm::AddHitsToClusters(const ClusterToHitMap &clustersToExpand,
    ClusterList &unavailableClusters) const
{
    for (ClusterToHitMap::const_iterator iterI = clustersToExpand.begin(), iterEndI = clustersToExpand.end(); iterI != iterEndI; ++iterI)
    {
        Cluster* pCluster = const_cast<Cluster*>(iterI->first);
        const CaloHitList &caloHitList = iterI->second;

        if (caloHitList.empty())
            continue;

        if (unavailableClusters.count(pCluster))
            continue;

        unavailableClusters.insert(pCluster);

        for (CaloHitList::const_iterator iterJ = caloHitList.begin(), iterEndJ = caloHitList.end(); iterJ != iterEndJ; ++iterJ)
        {
            CaloHit* pCaloHit = *iterJ;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pCluster, pCaloHit));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitConsolidationAlgorithm::RebuildClusters(const ClusterToHitMap &clustersToRebuild,
    const ClusterList &unavailableClusters) const
{
    if (clustersToRebuild.empty())
        return STATUS_CODE_SUCCESS;

    for (ClusterToHitMap::const_iterator iter = clustersToRebuild.begin(), iterEnd = clustersToRebuild.end(); iter != iterEnd; ++iter)
    {
        const Cluster* pCluster = iter->first;
        const CaloHitList &caloHitList = iter->second;

        Cluster* pClusterToDelete = const_cast<Cluster*>(pCluster);

        if (unavailableClusters.count(pClusterToDelete))
            continue;

        if (caloHitList.empty())
            continue;

        std::string currentClusterListName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*this, currentClusterListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete<Cluster>(*this, pClusterToDelete));

        const ClusterList *pClusterList = NULL;
        std::string newClusterListName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunClusteringAlgorithm(*this, m_reclusteringAlgorithmName,
            pClusterList, newClusterListName));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, newClusterListName, currentClusterListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, currentClusterListName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitConsolidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, xmlHandle,
        "ClusterRebuilding", m_reclusteringAlgorithmName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinTrackLength", m_minTrackLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxClusterLength", m_maxClusterLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitHalfWindow", m_halfWindowLayers));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
