/**
 *  @file   LArContent/src/LArTwoDReco/LArClusterSplitting/TrackConsolidationAlgorithm.cc
 *
 *  @brief  Implementation of the track consolidation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArTwoDReco/LArClusterSplitting/TrackConsolidationAlgorithm.h"

using namespace pandora;

namespace lar
{

void TrackConsolidationAlgorithm::GetReclusteredHits(const TwoDSlidingFitResultList &slidingFitResultListI,
    const ClusterVector &showerClustersJ, ClusterToHitMap &caloHitsToAddI, ClusterToHitMap &caloHitsToRemoveJ) const
{
    for (TwoDSlidingFitResultList::const_iterator iterI = slidingFitResultListI.begin(), iterEndI = slidingFitResultListI.end(); iterI != iterEndI; ++iterI)
    {
        const TwoDSlidingFitResult &slidingFitResultI = *iterI;

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

            this->GetReclusteredHits(slidingFitResultI, pClusterJ, caloHitsToAddI, caloHitsToRemoveJ);
        }
    }

    // Finalize the re-clustering
    this->FinalizeReclusteredHits(caloHitsToAddI, caloHitsToRemoveJ);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackConsolidationAlgorithm::GetReclusteredHits(const TwoDSlidingFitResult& slidingFitResultI, const Cluster* pClusterJ,
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

void TrackConsolidationAlgorithm::FinalizeReclusteredHits(ClusterToHitMap &caloHitsToAdd, ClusterToHitMap &caloHitsToRemove) const
{
    for (ClusterToHitMap::const_iterator iter1 = caloHitsToRemove.begin(), iterEnd1 = caloHitsToRemove.end(); iter1 != iterEnd1; ++iter1)
    {
        const Cluster* pCluster = iter1->first;

        if (caloHitsToAdd[pCluster].size() > 0)
            continue;

        CaloHitList caloHitList;
        pCluster->GetOrderedCaloHitList().GetCaloHitList(caloHitList);
        for (CaloHitList::const_iterator iter2 = caloHitList.begin(), iterEnd2 = caloHitList.end(); iter2 != iterEnd2; ++iter2)
        {
            CaloHit* pCaloHit = *iter2;

            if (caloHitsToRemove[pCluster].count(pCaloHit) > 0)
                continue;

            caloHitsToRemove[pCluster].insert(pCaloHit);
        }
    }

// --- EVENT DISPLAY [BEGIN] ---
// CaloHitList tempCaloHitsToStart, tempCaloHitsToAdd, tempCaloHitsToRebuild;
// ClusterList tempList;

// for (ClusterToHitMap::const_iterator iterI = caloHitsToAdd.begin(), iterEndI = caloHitsToAdd.end(); iterI != iterEndI; ++iterI)
// {
// const CaloHitList &caloHitListI = iterI->second;
// tempCaloHitsToAdd.insert(caloHitListI.begin(),caloHitListI.end());
// }

// for (ClusterToHitMap::const_iterator iterI = caloHitsToAdd.begin(), iterEndI = caloHitsToAdd.end(); iterI != iterEndI; ++iterI)
// {
// const Cluster* pClusterI = iterI->first;
// const CaloHitList &caloHitListI = iterI->second;
// if(caloHitListI.empty())
// continue;

// CaloHitList clusterHitsI;
// pClusterI->GetOrderedCaloHitList().GetCaloHitList(clusterHitsI);
// for (CaloHitList::const_iterator iterK = clusterHitsI.begin(), iterEndK = clusterHitsI.end(); iterK != iterEndK; ++iterK)
// {
// CaloHit* pCaloHit = *iterK;
// if (0 == tempCaloHitsToAdd.count(pCaloHit))
// tempCaloHitsToStart.insert((CaloHit*)pCaloHit);
// }
// }

// for (ClusterToHitMap::const_iterator iterJ = caloHitsToRemove.begin(), iterEndJ = caloHitsToRemove.end(); iterJ != iterEndJ; ++iterJ)
// {
// const Cluster* pClusterJ = iterJ->first;
// CaloHitList clusterHitsJ;
// pClusterJ->GetOrderedCaloHitList().GetCaloHitList(clusterHitsJ);
// for (CaloHitList::const_iterator iterK = clusterHitsJ.begin(), iterEndK = clusterHitsJ.end(); iterK != iterEndK; ++iterK)
// {
// CaloHit* pCaloHit = *iterK;
// if (0 == tempCaloHitsToAdd.count(pCaloHit) && 0==tempCaloHitsToStart.count(pCaloHit))
// tempCaloHitsToRebuild.insert((CaloHit*)pCaloHit);
// }
// }

// if (!tempCaloHitsToAdd.empty())
// {
// PandoraMonitoringApi::VisualizeCaloHits(&tempCaloHitsToStart, "Initial Clusters", BLUE);
// PandoraMonitoringApi::VisualizeCaloHits(&tempCaloHitsToAdd, "Hits To Add", RED);
// PandoraMonitoringApi::VisualizeCaloHits(&tempCaloHitsToRebuild, "Hits To Rebuild", GREEN);
// PandoraMonitoringApi::ViewEvent();
// }
// --- EVENT DISPLAY [END] ---
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackConsolidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_maxTransverseDisplacement = 1.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTransverseDisplacement", m_maxTransverseDisplacement));

    m_minAssociatedSpan = 1.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinAssociatedSpan", m_minAssociatedSpan));

    m_minAssociatedFraction = 0.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinAssociatedFraction", m_minAssociatedFraction));

    return TwoDSlidingFitConsolidationAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
