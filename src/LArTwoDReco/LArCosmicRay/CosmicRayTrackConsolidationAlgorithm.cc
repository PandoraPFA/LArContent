/**
 *  @file   LArContent/src/LArTwoDReco/LArCosmicRay/CosmicRayTrackConsolidationAlgorithm.cc
 *
 *  @brief  Implementation of the cosmic ray track consolidation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArTwoDReco/LArCosmicRay/CosmicRayTrackConsolidationAlgorithm.h"



using namespace pandora;

namespace lar
{

StatusCode CosmicRayTrackConsolidationAlgorithm::Run()
{
    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));


    ClusterVector trackClusters, shortClusters;
    this->SortInputClusters(pClusterList, trackClusters, shortClusters);
  

    TwoDSlidingFitResultList slidingFitResultList;
    this->BuildSlidingLinearFits(trackClusters, slidingFitResultList);
   

    ClusterToHitMap caloHitsToAdd, caloHitsToRemove;
    this->GetAssociatedHits(slidingFitResultList, shortClusters, caloHitsToAdd, caloHitsToRemove);



    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RunReclustering(caloHitsToAdd, caloHitsToRemove));
   

    // RE-CLUSTERING GOES HERE


    return STATUS_CODE_SUCCESS;
}
 
//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTrackConsolidationAlgorithm::SortInputClusters(const ClusterList *const pClusterList, ClusterVector &trackClusters, 
    ClusterVector &shortClusters) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster* pCluster = *iter;

        const float thisLengthSquared(LArClusterHelper::GetLengthSquared(pCluster));

        if (thisLengthSquared < m_maxClusterLength * m_maxClusterLength)
            shortClusters.push_back(pCluster);

        if (thisLengthSquared > m_minTrackLength * m_minTrackLength)
            trackClusters.push_back(pCluster);
    }

    std::sort(trackClusters.begin(), trackClusters.end(), LArClusterHelper::SortByNHits);
}
 
//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTrackConsolidationAlgorithm::BuildSlidingLinearFits(const ClusterVector &trackClusters, 
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
   
void CosmicRayTrackConsolidationAlgorithm::GetAssociatedHits(const TwoDSlidingFitResultList &slidingFitResultListI, 
    const ClusterVector &shortClustersJ, ClusterToHitMap &caloHitsToAddI, ClusterToHitMap &caloHitsToRemoveJ) const
{
    for (TwoDSlidingFitResultList::const_iterator iterI = slidingFitResultListI.begin(), iterEndI = slidingFitResultListI.end(); iterI != iterEndI; ++iterI)
    {
        const TwoDSlidingFitResult slidingFitResultI = *iterI;

        const Cluster* pClusterI = slidingFitResultI.GetCluster();
        const float thisLengthSquaredI(LArClusterHelper::GetLengthSquared(pClusterI));

        for (ClusterVector::const_iterator iterJ = shortClustersJ.begin(), iterEndJ = shortClustersJ.end(); iterJ != iterEndJ; ++iterJ)
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

void CosmicRayTrackConsolidationAlgorithm::GetAssociatedHits(const TwoDSlidingFitResult& slidingFitResultI, const Cluster* pClusterJ,
    ClusterToHitMap &caloHitsToAddI, ClusterToHitMap &caloHitsToRemoveJ) const
{

    const Cluster* pClusterI(slidingFitResultI.GetCluster());
   
    CaloHitList associatedHits, caloHitListJ;
    pClusterJ->GetOrderedCaloHitList().GetCaloHitList(caloHitListJ);

    float minL(std::numeric_limits<float>::max());
    float maxL(std::numeric_limits<float>::max());

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

StatusCode CosmicRayTrackConsolidationAlgorithm::RunReclustering(const ClusterToHitMap &clustersToSave, const ClusterToHitMap &clustersToDelete) const
{

    for (ClusterToHitMap::const_iterator iterI = clustersToDelete.begin(), iterEndI = clustersToDelete.end(); iterI != iterEndI; ++iterI)
    {
        Cluster* pClusterToDelete = const_cast<Cluster*>(iterI->first);
        
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete<Cluster>(*this, pClusterToDelete));
    }

    for (ClusterToHitMap::const_iterator iterI = clustersToSave.begin(), iterEndI = clustersToSave.end(); iterI != iterEndI; ++iterI)
    {
        Cluster* pClusterToSave = const_cast<Cluster*>(iterI->first);
        const CaloHitList &caloHitList = iterI->second;

        for (CaloHitList::const_iterator iterJ = caloHitList.begin(), iterEndJ = caloHitList.end(); iterJ != iterEndJ; ++iterJ)
        {
            CaloHit* pCaloHit = *iterJ;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pClusterToSave, pCaloHit));
        }
    }


    // Run the initial cluster formation algorithm
    

    const ClusterList *pClusterList = NULL;
    std::string currentClusterListName, newClusterListName;

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*this, currentClusterListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunClusteringAlgorithm(*this, m_clusteringAlgorithmName,
        pClusterList, newClusterListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, newClusterListName, currentClusterListName));



    /// -----------------------------

    /*
    PRR_IF PandoraContentApi::Delete(*this, pCluster, "ClusterListName");
    
    PRR_IF PandoraContentApi::AddToCluster(*this, pCluster, pCaloHit);

    // Find example, e.g. ClusteringParentAlgorithm to instantiate daughter clustering alg and receive name
    PRR_IF PandoraContentApi::RunClusteringAlgorithm(*this, m_clusteringAlgorithmName);

    // No daughter clustering alg?
    PRR_IF PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, ...);
    // Then proceed as if you were in the daughter clustering alg

    PRR_IF PandoraContentApi::SaveList<Cluster>(*this, m_chosenClusterListName);
    */
    /// -----------------------------

  
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayTrackConsolidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_minTrackLength = 7.5f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinTrackLength", m_minTrackLength));

    m_maxClusterLength = 15.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxClusterLength", m_maxClusterLength));

    m_maxTransverseDisplacement = 1.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTransverseDisplacement", m_maxTransverseDisplacement));

    m_minAssociatedSpan = 1.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinAssociatedSpan", m_minAssociatedSpan));

    m_minAssociatedFraction = 0.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinAssociatedFraction", m_minAssociatedFraction));

    m_halfWindowLayers = 20;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitHalfWindow", m_halfWindowLayers));


    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, xmlHandle, "Clustering",
        m_clusteringAlgorithmName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
