/**
 *  @file   LArContent/src/LArTwoDReco/LArClusterSplitting/TwoDSlidingFitConsolidationAlgorithm.cc
 *
 *  @brief  Implementation of the 2D sliding fit consolidation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/TwoDSlidingFitConsolidationAlgorithm.h"

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
        const Cluster *const pCluster = *iter;

        const float thisLengthSquared(LArClusterHelper::GetLengthSquared(pCluster));

        if (thisLengthSquared < m_maxClusterLength * m_maxClusterLength)
            showerClusters.push_back(pCluster);

        if (thisLengthSquared > m_minTrackLength * m_minTrackLength)
            trackClusters.push_back(pCluster);
    }

    std::sort(trackClusters.begin(), trackClusters.end(), LArClusterHelper::SortByNHits);
    std::sort(showerClusters.begin(), showerClusters.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitConsolidationAlgorithm::BuildSlidingLinearFits(const ClusterVector &trackClusters, TwoDSlidingFitResultList &slidingFitResultList) const
{
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));

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

StatusCode TwoDSlidingFitConsolidationAlgorithm::RemoveHitsFromClusters(const ClusterToHitMap &clustersToContract, ClusterList &unavailableClusters) const
{
    for (const ClusterToHitMap::value_type &mapEntry : clustersToContract)
    {
        const Cluster *const pCluster = mapEntry.first;
        const CaloHitList &caloHitListToRemove = mapEntry.second;

        if (caloHitListToRemove.empty())
            continue;

        if (unavailableClusters.count(pCluster))
            continue;

        CaloHitList caloHitList, caloHitListToKeep;
        pCluster->GetOrderedCaloHitList().GetCaloHitList(caloHitList);

        for (const CaloHit *const pCaloHit : caloHitList)
        {
            if (!caloHitListToRemove.count(pCaloHit))
                caloHitListToKeep.insert(pCaloHit);
        }

        if (caloHitListToKeep.empty())
        {
            // ATTN clustersToContract and unavailable clusters now contain dangling pointers
            unavailableClusters.insert(pCluster);
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete<Cluster>(*this, pCluster));
            continue;
        }

        for (const CaloHit *const pCaloHit : caloHitListToRemove)
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromCluster(*this, pCluster, pCaloHit));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitConsolidationAlgorithm::AddHitsToClusters(const ClusterToHitMap &clustersToExpand, ClusterList &unavailableClusters) const
{
    for (const ClusterToHitMap::value_type &mapEntry : clustersToExpand)
    {
        const Cluster *const pCluster = mapEntry.first;
        const CaloHitList &caloHitList = mapEntry.second;

        if (caloHitList.empty())
            continue;

        if (unavailableClusters.count(pCluster))
            continue;

        unavailableClusters.insert(pCluster);

        for (const CaloHit *const pCaloHit : caloHitList)
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pCluster, pCaloHit));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitConsolidationAlgorithm::RebuildClusters(const ClusterToHitMap &clustersToRebuild, const ClusterList &unavailableClusters) const
{
    if (clustersToRebuild.empty())
        return STATUS_CODE_SUCCESS;

    ClusterVector sortedClusters;
    for (const auto &mapEntry : clustersToRebuild)
    {
        if (!unavailableClusters.count(mapEntry.first))
            sortedClusters.push_back(mapEntry.first);
    }

    std::sort(sortedClusters.begin(), sortedClusters.end(), LArClusterHelper::SortByNHits);

    for (const Cluster *const pCluster : sortedClusters)
    {
        const CaloHitList &caloHitList(clustersToRebuild.at(pCluster));
        const Cluster *const pClusterToDelete(pCluster);

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
