/**
 *  @file   LArContent/src/LArTwoDReco/LArClusterSplitting/TwoDSlidingFitConsolidationAlgorithm.cc
 *
 *  @brief  Implementation of the 2D sliding fit consolidation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArTwoDReco/LArClusterSplitting/TwoDSlidingFitConsolidationAlgorithm.h"

using namespace pandora;

namespace lar
{

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

    // Get Initial configuration of hits
    ClusterToHitMap clustersAtStart;
    this->GetInitialHits(trackClusters, clustersAtStart);
    this->GetInitialHits(showerClusters, clustersAtStart);

    // Recluster the hits
    ClusterToHitMap clustersToExpand, clustersToContract;
    this->GetReclusteredHits(slidingFitResultList, showerClusters, clustersToExpand, clustersToContract);

    // Consolidate and re-build clusters
    ClusterList modifiedShowers, modifiedTracks;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RemoveHitsFromClusters(clustersToContract));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->AddHitsToClusters(clustersToExpand));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RebuildClusters(clustersAtStart));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitConsolidationAlgorithm::SortInputClusters(const ClusterList *const pClusterList, ClusterVector &trackClusters,
    ClusterVector &showerClusters) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster* pCluster = *iter;

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
    for (ClusterVector::const_iterator iter = trackClusters.begin(), iterEnd = trackClusters.end(); iter != iterEnd; ++iter)
    {
         TwoDSlidingFitResult slidingFitResult;
         LArClusterHelper::LArTwoDSlidingFit(*iter, m_halfWindowLayers, slidingFitResult);
         slidingFitResultList.push_back(slidingFitResult);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitConsolidationAlgorithm::GetInitialHits(const ClusterVector &inputClusters, ClusterToHitMap &caloHitsAtStart) const
{
    for (ClusterVector::const_iterator iter = inputClusters.begin(), iterEnd = inputClusters.end(); iter != iterEnd; ++iter)
    {
        const Cluster* pCluster = *iter;

        if (caloHitsAtStart[pCluster].size() > 0)
            continue;

        CaloHitList caloHitList;
        pCluster->GetOrderedCaloHitList().GetCaloHitList(caloHitList);
        caloHitsAtStart[pCluster].insert(caloHitList.begin(), caloHitList.end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitConsolidationAlgorithm::RemoveHitsFromClusters(const ClusterToHitMap &clustersToContract) const
{
    for (ClusterToHitMap::const_iterator iterI = clustersToContract.begin(), iterEndI = clustersToContract.end(); iterI != iterEndI; ++iterI)
    {
        Cluster* pCluster = const_cast<Cluster*>(iterI->first);
        const CaloHitList &caloHitListToRemove = iterI->second;

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

StatusCode TwoDSlidingFitConsolidationAlgorithm::AddHitsToClusters(const ClusterToHitMap &clustersToExpand) const
{
    for (ClusterToHitMap::const_iterator iterI = clustersToExpand.begin(), iterEndI = clustersToExpand.end(); iterI != iterEndI; ++iterI)
    {
        Cluster* pCluster = const_cast<Cluster*>(iterI->first);
        const CaloHitList &caloHitList = iterI->second;

        for (CaloHitList::const_iterator iterJ = caloHitList.begin(), iterEndJ = caloHitList.end(); iterJ != iterEndJ; ++iterJ)
        {
            CaloHit* pCaloHit = *iterJ;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pCluster, pCaloHit));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitConsolidationAlgorithm::RebuildClusters(const ClusterToHitMap &clustersAtStart) const
{
    ClusterToHitMap clustersToRebuild;

    for (ClusterToHitMap::const_iterator iterI = clustersAtStart.begin(), iterEndI = clustersAtStart.end(); iterI != iterEndI; ++iterI)
    {
        Cluster* pCluster = const_cast<Cluster*>(iterI->first);
        const CaloHitList &caloHitList = iterI->second;

        for (CaloHitList::const_iterator iterJ = caloHitList.begin(), iterEndJ = caloHitList.end(); iterJ != iterEndJ; ++iterJ)
        {
            CaloHit* pCaloHit = *iterJ;

            if (PandoraContentApi::IsAvailable(*this, pCaloHit))
                clustersToRebuild[pCluster].insert(pCaloHit);
        }
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->BuildNewClusters(clustersToRebuild));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitConsolidationAlgorithm::BuildNewClusters(const ClusterToHitMap &clustersToRebuild) const
{
    if (clustersToRebuild.empty())
        return STATUS_CODE_SUCCESS;

    const ClusterList *pClusterList = NULL;
    std::string currentClusterListName, newClusterListName;

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*this, currentClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent<ClusterList>(*this, pClusterList, newClusterListName));

    for (ClusterToHitMap::const_iterator iter = clustersToRebuild.begin(), iterEnd = clustersToRebuild.end(); iter != iterEnd; ++iter)
    {
        const CaloHitList &caloHitList = iter->second;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->BuildNewClusters(caloHitList));
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, newClusterListName, currentClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, currentClusterListName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitConsolidationAlgorithm::BuildNewClusters(const CaloHitList &inputCaloHitList) const
{
    if (inputCaloHitList.empty())
        return STATUS_CODE_SUCCESS;

    // Form simple associations between residual hits from deleted cluster
    HitToHitMap hitAssociationMap;

    for (CaloHitList::const_iterator iterI = inputCaloHitList.begin(), iterEndI = inputCaloHitList.end(); iterI != iterEndI; ++iterI)
    {
        CaloHit* pCaloHitI = *iterI;

        for (CaloHitList::const_iterator iterJ = iterI, iterEndJ = iterEndI; iterJ != iterEndJ; ++iterJ)
        {
            CaloHit* pCaloHitJ = *iterJ;

            if (pCaloHitI == pCaloHitJ)
                continue;

            if ((pCaloHitI->GetPositionVector() - pCaloHitJ->GetPositionVector()).GetMagnitudeSquared() < m_reclusteringWindow * m_reclusteringWindow)
            {
                hitAssociationMap[pCaloHitI].insert(pCaloHitJ);
                hitAssociationMap[pCaloHitJ].insert(pCaloHitI);
            }
        }
    }

    // Collect up associations and build new clusters
    CaloHitList vetoList;

    for (CaloHitList::const_iterator iterI = inputCaloHitList.begin(), iterEndI = inputCaloHitList.end(); iterI != iterEndI; ++iterI)
    {
        CaloHit* pSeedCaloHit = *iterI;

        if (vetoList.count(pSeedCaloHit))
            continue;

        CaloHitList mergeList;
        this->CollectAssociatedHits(pSeedCaloHit, pSeedCaloHit, hitAssociationMap, vetoList, mergeList);

        Cluster *pCluster = NULL;
        PandoraContentApi::Cluster::Parameters parameters;
        parameters.m_caloHitList.insert(pSeedCaloHit);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pCluster));
        vetoList.insert(pSeedCaloHit);

        for (CaloHitList::const_iterator iterJ = mergeList.begin(), iterEndJ = mergeList.end(); iterJ != iterEndJ; ++iterJ)
        {
            CaloHit* pAssociatedCaloHit = *iterJ;

            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pCluster, pAssociatedCaloHit));
            vetoList.insert(pAssociatedCaloHit);
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitConsolidationAlgorithm::CollectAssociatedHits(CaloHit *pSeedCaloHit, CaloHit *pCurrentCaloHit,
    const HitToHitMap &hitAssociationMap, const CaloHitList &vetoList, CaloHitList &mergeList) const
{
    if (vetoList.count(pCurrentCaloHit))
        return;

    HitToHitMap::const_iterator iter1 = hitAssociationMap.find(pCurrentCaloHit);
    if (iter1 == hitAssociationMap.end())
        return;

    for (CaloHitList::const_iterator iter2 = iter1->second.begin(), iterEnd2 = iter1->second.end(); iter2 != iterEnd2; ++iter2)
    {
        CaloHit* pAssociatedCaloHit = *iter2;

        if (pAssociatedCaloHit == pSeedCaloHit)
            continue;

        if (!mergeList.insert(pAssociatedCaloHit).second)
            continue;

        this->CollectAssociatedHits(pSeedCaloHit, pAssociatedCaloHit, hitAssociationMap, vetoList, mergeList);
    }

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitConsolidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_minTrackLength = 7.5f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinTrackLength", m_minTrackLength));

    m_maxClusterLength = 15.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxClusterLength", m_maxClusterLength));

    m_halfWindowLayers = 25;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitHalfWindow", m_halfWindowLayers));

    m_reclusteringWindow = 2.5f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ReclusteringWindow", m_reclusteringWindow));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
