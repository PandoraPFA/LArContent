/**
 *  @file   LArContent/src/LArTwoDReco/LArClusterSplitting/TwoDTrackConsolidationAlgorithm.cc
 *
 *  @brief  Implementation of the 2D track consolidation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArTwoDReco/LArClusterSplitting/TwoDTrackConsolidationAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode TwoDTrackConsolidationAlgorithm::Run()
{
    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

    // Select tracks and showers for re-clustering (Note: samples not mutually exclusive)
    ClusterVector trackClusters, showerClusters;
    this->SortInputClusters(pClusterList, trackClusters, showerClusters);

    // Build sliding linear fits from track clusters
    TwoDSlidingFitResultList slidingFitResultList;
    this->BuildSlidingLinearFits(trackClusters, slidingFitResultList);

    // Form associations: decide which hits should be moved from shower clusters into track clusters
    ClusterToHitMap clustersToExpand, clustersToContract;
    this->GetAssociatedHits(slidingFitResultList, showerClusters, clustersToExpand, clustersToContract);

    // Consolidate and re-build clusters
    ClusterList modifiedShowers, modifiedTracks;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RemoveHitsFromShowers(clustersToContract, modifiedShowers));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->AddHitsToTracks(clustersToExpand, modifiedTracks));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RebuildClusters(modifiedTracks, modifiedShowers));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDTrackConsolidationAlgorithm::SortInputClusters(const ClusterList *const pClusterList, ClusterVector &trackClusters,
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

void TwoDTrackConsolidationAlgorithm::BuildSlidingLinearFits(const ClusterVector &trackClusters,
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

void TwoDTrackConsolidationAlgorithm::GetAssociatedHits(const TwoDSlidingFitResultList &slidingFitResultListI,
    const ClusterVector &showerClustersJ, ClusterToHitMap &caloHitsToAddI, ClusterToHitMap &caloHitsToRemoveJ) const
{
    for (TwoDSlidingFitResultList::const_iterator iterI = slidingFitResultListI.begin(), iterEndI = slidingFitResultListI.end(); iterI != iterEndI; ++iterI)
    {
        const TwoDSlidingFitResult slidingFitResultI = *iterI;

        const Cluster* pClusterI = slidingFitResultI.GetCluster();
        const float thisLengthSquaredI(LArClusterHelper::GetLengthSquared(pClusterI));

        for (ClusterVector::const_iterator iterJ = showerClustersJ.begin(), iterEndJ = showerClustersJ.end(); iterJ != iterEndJ; ++iterJ)
        {
            const Cluster* pClusterJ = *iterJ;
            const float thisLengthSquaredJ(LArClusterHelper::GetLengthSquared(pClusterJ));

            if (pClusterI == pClusterJ)
                continue;

            if (2.0 * thisLengthSquaredJ > thisLengthSquaredI)
                continue;

            this->GetAssociatedHits(slidingFitResultI, pClusterJ, caloHitsToAddI, caloHitsToRemoveJ);
        }
    }

// --- EVENT DISPLAY [BEGIN] ---
// ClusterList tempClusters;
// CaloHitList tempCaloHitsToKeep, tempCaloHitsToMove;

// for (ClusterToHitMap::const_iterator iterI = caloHitsToAddI.begin(), iterEndI = caloHitsToAddI.end(); iterI != iterEndI; ++iterI)
// {
// const Cluster* pClusterI = iterI->first;
// const CaloHitList &caloHitListI = iterI->second;

// tempClusters.insert((Cluster*)pClusterI);
// tempCaloHitsToMove.insert(caloHitListI.begin(), caloHitListI.end());
// }

// for (ClusterToHitMap::const_iterator iterJ = caloHitsToRemoveJ.begin(), iterEndJ = caloHitsToRemoveJ.end(); iterJ != iterEndJ; ++iterJ)
// {
// const Cluster* pClusterJ = iterJ->first;
// const CaloHitList &caloHitListMove = iterJ->second;

// CaloHitList caloHitListKeep;
// pClusterJ->GetOrderedCaloHitList().GetCaloHitList(caloHitListKeep);
// for (CaloHitList::const_iterator iterK = caloHitListKeep.begin(), iterEndK = caloHitListKeep.end(); iterK != iterEndK; ++iterK)
// {
// CaloHit* pCaloHit = *iterK;
// if (0 == caloHitListMove.count(pCaloHit))
// tempCaloHitsToKeep.insert((CaloHit*)pCaloHit);
// }
// }

// if (!tempClusters.empty())
// {
// PandoraMonitoringApi::VisualizeClusters(&tempClusters, "Cluster", BLUE);
// PandoraMonitoringApi::VisualizeCaloHits(&tempCaloHitsToMove, "HitsToMove", RED);
// PandoraMonitoringApi::VisualizeCaloHits(&tempCaloHitsToKeep, "HitsToKeep", GREEN);
// PandoraMonitoringApi::ViewEvent();
// }
// --- EVENT DISPLAY [END] ---
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDTrackConsolidationAlgorithm::GetAssociatedHits(const TwoDSlidingFitResult& slidingFitResultI, const Cluster* pClusterJ,
    ClusterToHitMap &caloHitsToAddI, ClusterToHitMap &caloHitsToRemoveJ) const
{
    const Cluster* pClusterI(slidingFitResultI.GetCluster());

    CaloHitList associatedHits, caloHitListJ;
    pClusterJ->GetOrderedCaloHitList().GetCaloHitList(caloHitListJ);

    float minL(std::numeric_limits<float>::max());
    float maxL(std::numeric_limits<float>::max());

    // Loop over hits from shower clusters, and make associations with track clusters
    // (Determine if hits from shower clusters can be used to fill gaps in track cluster)
    //
    // Apply the following selection:
    //  rJ = candidate hit from shower cluster
    //  rI = nearest hit on track cluster
    //  rK = projection of shower hit onto track cluster
    //
    //                  o rJ
    //  o o o o o o - - x - - - - o o o o o o o
    //           rI    rK
    //
    //  Require: rJK < std::min(rCut, rIJ, rKI)

    for (CaloHitList::const_iterator iterJ = caloHitListJ.begin(), iterEndJ = caloHitListJ.end(); iterJ != iterEndJ; ++iterJ)
    {
        CaloHit* pCaloHitJ = *iterJ;

        try
        {
            const CartesianVector positionJ(pCaloHitJ->GetPositionVector());
            const CartesianVector positionI(LArClusterHelper::GetClosestPosition(positionJ, pClusterI));

            float rL(0.f), rT(0.f);
            CartesianVector positionK(0.f, 0.f, 0.f);
            slidingFitResultI.GetGlobalFitProjection(positionJ, positionK);
            slidingFitResultI.GetLocalPosition(positionK, rL, rT);

            const float rsqIJ((positionI - positionJ).GetMagnitudeSquared());
            const float rsqJK((positionJ - positionK).GetMagnitudeSquared());
            const float rsqKI((positionK - positionI).GetMagnitudeSquared());

            if (rsqJK < std::min(m_maxTransverseDisplacement * m_maxTransverseDisplacement, std::min(rsqIJ, rsqKI)))
            {
                if (associatedHits.empty())
                {
                    minL = rL;
                    maxL = rL;
                }
                else
                {
                    minL = std::min(minL, rL);
                    maxL = std::max(maxL, rL);
                }

                associatedHits.insert(pCaloHitJ);
            }
        }
        catch (StatusCodeException &)
        {
        }
    }

    const float associatedSpan(maxL - minL);
    const float associatedFraction(associatedHits.empty() ? 0.f : static_cast<float>(associatedHits.size()) / static_cast<float>(pClusterJ->GetNCaloHits()));

    if (associatedSpan > m_minAssociatedSpan || associatedFraction > m_minAssociatedFraction)
    {
        for (CaloHitList::const_iterator iterK = associatedHits.begin(), iterEndK = associatedHits.end(); iterK != iterEndK; ++iterK)
        {
            CaloHit* pCaloHit = *iterK;

            if (caloHitsToRemoveJ[pClusterJ].count(pCaloHit))
                continue;

            caloHitsToAddI[pClusterI].insert(pCaloHit);
            caloHitsToRemoveJ[pClusterJ].insert(pCaloHit);
        }
    }
}
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDTrackConsolidationAlgorithm::RemoveHitsFromShowers(const ClusterToHitMap &clustersToContract, ClusterList &modifiedShowers) const
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

        modifiedShowers.insert(pCluster);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDTrackConsolidationAlgorithm::AddHitsToTracks(const ClusterToHitMap &clustersToExpand, ClusterList &modifiedTracks) const
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

        modifiedTracks.insert(pCluster);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDTrackConsolidationAlgorithm::RebuildClusters(const ClusterList &tracksToRebuild, const ClusterList &showersToRebuild) const
{
    ClusterToHitMap clustersToRebuild;

    for (ClusterList::const_iterator iter = showersToRebuild.begin(), iterEnd = showersToRebuild.end(); iter != iterEnd; ++iter)
    {
        Cluster* pCluster = *iter;

        if (tracksToRebuild.count(pCluster))
            continue;

        CaloHitList caloHitList;
        pCluster->GetOrderedCaloHitList().GetCaloHitList(caloHitList);
        clustersToRebuild[pCluster].insert(caloHitList.begin(), caloHitList.end());

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete<Cluster>(*this, pCluster));
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->BuildNewClusters(clustersToRebuild));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDTrackConsolidationAlgorithm::BuildNewClusters(const ClusterToHitMap &clustersToRebuild) const
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

StatusCode TwoDTrackConsolidationAlgorithm::BuildNewClusters(const CaloHitList &inputCaloHitList) const
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
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, pSeedCaloHit, pCluster));
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

void TwoDTrackConsolidationAlgorithm::CollectAssociatedHits(CaloHit *pSeedCaloHit, CaloHit *pCurrentCaloHit,
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

StatusCode TwoDTrackConsolidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_minTrackLength = 7.5f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinTrackLength", m_minTrackLength));

    m_maxClusterLength = 15.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxClusterLength", m_maxClusterLength));

     m_halfWindowLayers = 20;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitHalfWindow", m_halfWindowLayers));

    m_maxTransverseDisplacement = 1.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTransverseDisplacement", m_maxTransverseDisplacement));

    m_minAssociatedSpan = 1.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinAssociatedSpan", m_minAssociatedSpan));

    m_minAssociatedFraction = 0.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinAssociatedFraction", m_minAssociatedFraction));

    m_reclusteringWindow = 2.5f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ReclusteringWindow", m_reclusteringWindow));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
