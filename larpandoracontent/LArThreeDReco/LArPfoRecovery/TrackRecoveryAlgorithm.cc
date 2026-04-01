/**
 *  @file   larpandoracontent/LArThreeDReco/LArPfoRecovery/TrackRecoveryAlgorithm.cc
 *
 *  @brief  Implementation of the particle recovery algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArPlaneContextObject.h"

#include "larpandoracontent/LArThreeDReco/LArPfoRecovery/TrackRecoveryAlgorithm.h"

#include <algorithm>

using namespace pandora;

namespace lar_content
{

TrackRecoveryAlgorithm::TrackRecoveryAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackRecoveryAlgorithm::Run()
{
    // Get the list of track-like PFOs
    const PfoList *pPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputPfoListName, pPfoList));

    // Get maps between PFOs, views, clusters and hits
    PfoToViewClusterMap pfoToViewClusterMap;
    ClusterToPfoMap clusterToPfoMap;
    ClusterToHitMap clusterToHitMap;
    ClusterToFitMap clusterToFitMap;
    HitToClusterToMap hitToClusterMap;
    for (const Pfo *const pPfo : *pPfoList)
    {
        for (const HitType view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, TPC_3D})
        {
            ClusterList pfoClusters;
            LArPfoHelper::GetClusters(pPfo, view, pfoClusters);
            if (!pfoClusters.empty())
            {
                const Cluster *pCluster{pfoClusters.front()};
                pfoToViewClusterMap[pPfo][view] = pCluster;
                clusterToPfoMap[pCluster] = pPfo;
                if (view != TPC_3D)
                {
                    auto result = this->FitAndOrderCluster(pCluster, clusterToHitMap[pCluster]);
                    if (result)
                        clusterToFitMap.emplace(pCluster, *result);
                }
                else
                    LArClusterHelper::GetAllHits(pCluster, clusterToHitMap[pCluster]);
                for (const CaloHit *const pCaloHit : clusterToHitMap[pCluster])
                    hitToClusterMap[pCaloHit] = pCluster;
            }
        }
    }
    // Loop over all clusters to catch any that aren't associated to a PFO
    for (const std::string &clusterListName : m_inputClusterListNames)
    {
        const ClusterList *pClusterList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));
        for (const Cluster *const pCluster : *pClusterList)
        {
            if (clusterToPfoMap.find(pCluster) == clusterToPfoMap.end())
            {
                LArClusterHelper::GetAllHits(pCluster, clusterToHitMap[pCluster]);
                for (const CaloHit *const pCaloHit : clusterToHitMap[pCluster])
                    hitToClusterMap[pCaloHit] = pCluster;
            }
        }
    }

    // Loop over the PFOs and project 3D hits from two views into the third view to look for unassociated hits
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    std::unordered_map<const Pfo *, std::unordered_map<HitType, CaloHitList>> pfoToMergeHitsMap;
    for (const Pfo *const pPfo : *pPfoList)
    {
        const Cluster *const pClusterU{pfoToViewClusterMap[pPfo][TPC_VIEW_U]}, *const pClusterV{pfoToViewClusterMap[pPfo][TPC_VIEW_V]},
            *const pClusterW{pfoToViewClusterMap[pPfo][TPC_VIEW_W]};
        const CaloHitList &hitsU{clusterToHitMap[pClusterU]}, &hitsV{clusterToHitMap[pClusterV]}, &hitsW{clusterToHitMap[pClusterW]};

        CaloHitSet unmatchedHitsU, unmatchedHitsV, unmatchedHitsW;
        CaloHitList mergeHitsU, mergeHitsV, mergeHitsW;
        this->FindUnmatchedHits(hitsU, hitsV, hitsW, unmatchedHitsV, unmatchedHitsW);
        this->FindUnmatchedHits(hitsV, hitsU, hitsW, unmatchedHitsU, unmatchedHitsW);
        this->FindUnmatchedHits(hitsW, hitsU, hitsV, unmatchedHitsU, unmatchedHitsV);

        // Construct 3D hits from the combinatorics of each view pair and project into the third view
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &hitsU, "U", RED));
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &hitsV, "V", GREEN));
        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &hitsW, "W", BLUE));

        this->IdentifyHitsToMerge(pClusterU, hitsU, unmatchedHitsU, clusterToFitMap, mergeHitsU);
        this->IdentifyHitsToMerge(pClusterV, hitsV, unmatchedHitsV, clusterToFitMap, mergeHitsV);
        this->IdentifyHitsToMerge(pClusterW, hitsW, unmatchedHitsW, clusterToFitMap, mergeHitsW);

        CartesianPointVector newPosU;
        for (const CaloHit *const pCaloHit : hitsU)
            newPosU.emplace_back(pCaloHit->GetPositionVector());
        for (const CaloHit *const pCaloHit : mergeHitsU)
            newPosU.emplace_back(pCaloHit->GetPositionVector());
        TwoDSlidingFitResult sfrU(&newPosU, 3, LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_U));
        CartesianPointVector newPosV;
        for (const CaloHit *const pCaloHit : hitsV)
            newPosV.emplace_back(pCaloHit->GetPositionVector());
        for (const CaloHit *const pCaloHit : mergeHitsV)
            newPosV.emplace_back(pCaloHit->GetPositionVector());
        TwoDSlidingFitResult sfrV(&newPosV, 3, LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_V));
        CartesianPointVector newPosW;
        for (const CaloHit *const pCaloHit : hitsW)
            newPosW.emplace_back(pCaloHit->GetPositionVector());
        for (const CaloHit *const pCaloHit : mergeHitsW)
            newPosW.emplace_back(pCaloHit->GetPositionVector());
        TwoDSlidingFitResult sfrW(&newPosW, 3, LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_W));

        this->FilterHitsToMerge(sfrU, mergeHitsU);
        this->FilterHitsToMerge(sfrV, mergeHitsV);
        this->FilterHitsToMerge(sfrW, mergeHitsW);

        for (const CaloHit *const pCaloHit : mergeHitsU)
        {
            pfoToMergeHitsMap[pPfo][TPC_VIEW_U].emplace_back(pCaloHit);
            const CartesianVector position{pCaloHit->GetPositionVector()};
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &position, "Merge U", RED, 2.f));
        }
        for (const CaloHit *const pCaloHit : mergeHitsV)
        {
            pfoToMergeHitsMap[pPfo][TPC_VIEW_V].emplace_back(pCaloHit);
            const CartesianVector position{pCaloHit->GetPositionVector()};
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &position, "Merge V", GREEN, 2.f));
        }
        for (const CaloHit *const pCaloHit : mergeHitsW)
        {
            pfoToMergeHitsMap[pPfo][TPC_VIEW_W].emplace_back(pCaloHit);
            const CartesianVector position{pCaloHit->GetPositionVector()};
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &position, "Merge W", BLUE, 2.f));
        }
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    std::map<const CaloHit *, int> hitToNumMergesMap;
    for (auto &[pPfo, viewToMergeHitsMap] : pfoToMergeHitsMap)
    {
        for (auto &[view, mergeHits] : viewToMergeHitsMap)
        {
            const Cluster *pNewCluster{pfoToViewClusterMap[pPfo][view]};
            // Need to allow for the possibility that the view was entirely missing, in which case we need to create a new cluster to merge into
            for (const CaloHit *const pCaloHit : mergeHits)
            {
                if (hitToNumMergesMap.find(pCaloHit) == hitToNumMergesMap.end())
                    hitToNumMergesMap[pCaloHit] = 0;
                hitToNumMergesMap[pCaloHit]++;
                if (hitToClusterMap.find(pCaloHit) != hitToClusterMap.end())
                {
                    const Cluster *pOldCluster{hitToClusterMap[pCaloHit]};
                    PandoraContentApi::RemoveFromCluster(*this, pOldCluster, pCaloHit);
                    CaloHitList &oldClusterHits{clusterToHitMap[pOldCluster]};
                    std::remove(oldClusterHits.begin(), oldClusterHits.end(), pCaloHit);
                }
                CaloHitList preHitList;
                LArClusterHelper::GetAllHits(pNewCluster, preHitList);
                PandoraContentApi::AddToCluster(*this, pNewCluster, pCaloHit);
                CaloHitList postHitList;
                LArClusterHelper::GetAllHits(pNewCluster, postHitList);
                hitToClusterMap[pCaloHit] = pNewCluster;
            }
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::optional<TwoDSlidingFitResult> TrackRecoveryAlgorithm::FitAndOrderCluster(const Cluster *const pCluster, CaloHitList &orderedHits) const
{
    if (pCluster->GetNCaloHits() < 3)
    {
        LArClusterHelper::GetAllHits(pCluster, orderedHits);
        return std::nullopt;
    }
    const HitType view{LArClusterHelper::GetClusterHitType(pCluster)};
    TwoDSlidingFitResult sfr(pCluster, 3, LArGeometryHelper::GetWirePitch(this->GetPandora(), view));
    LArClusterHelper::OrderHitsAlongTrajectory(pCluster, sfr, orderedHits);

    return sfr;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackRecoveryAlgorithm::FindUnmatchedHits(const CaloHitList &hitsA, const CaloHitList &hitsB, const CaloHitList &hitsC, CaloHitSet &unmatchedHitsB,
    CaloHitSet &unmatchedHitsC) const
{
    if (hitsA.empty() || hitsB.empty() || hitsC.empty())
        return;
    const LArPlaneContextObject *pPlaneContextObject{dynamic_cast<const LArPlaneContextObject*>(PandoraContentApi::GetEventContextObject(
        *this, "PlaneContext"))};
    const HitType viewB{hitsB.front()->GetHitType()}, viewC{hitsC.front()->GetHitType()};
    for (const CaloHit *const pCaloHitA : hitsA)
    {
        const LArPlaneContextObject::HitTriplet *pTriplet{pPlaneContextObject->GetHitTriplet(pCaloHitA)};
        if (!pTriplet)
            continue;
        const CaloHit *const pHitB{viewB == TPC_VIEW_W ? pTriplet->m_wHit : viewB == TPC_VIEW_V ? pTriplet->m_vHit : pTriplet->m_uHit};
        const CaloHit *const pHitC{viewC == TPC_VIEW_W ? pTriplet->m_wHit : viewC == TPC_VIEW_V ? pTriplet->m_vHit : pTriplet->m_uHit};
        auto bIterator{std::find(hitsB.begin(), hitsB.end(), pHitB)};
        auto cIterator{std::find(hitsC.begin(), hitsC.end(), pHitC)};
        // Only flag hits as unmatched if the other hits in the triplet are present in the PFO
        if (bIterator == hitsB.end() && cIterator != hitsC.end())
            unmatchedHitsB.insert(pHitB);
        if (cIterator == hitsC.end() && bIterator != hitsB.end())
            unmatchedHitsC.insert(pHitC);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackRecoveryAlgorithm::IdentifyHitsToMerge(const Cluster *pCluster, const CaloHitList &clusterHits, const CaloHitSet &unmatchedHits,
    const ClusterToFitMap &clusterToFitMap, CaloHitList &mergeHits) const
{
    if (pCluster && (clusterHits.size() > 1) && (clusterToFitMap.find(pCluster) != clusterToFitMap.end()))
    {
        float rLFront{0.f}, rLBack{0.f}, dummy{0.f};
        clusterToFitMap.at(pCluster).GetLocalPosition(clusterHits.front()->GetPositionVector(), rLFront, dummy);
        clusterToFitMap.at(pCluster).GetLocalPosition(clusterHits.back()->GetPositionVector(), rLBack, dummy);
        for (const CaloHit *const pCaloHit : unmatchedHits)
        {
            float rL{0.f};
            clusterToFitMap.at(pCluster).GetLocalPosition(pCaloHit->GetPositionVector(), rL, dummy);
            if (rLFront <= rL && rL <= rLBack)
            {
                // We're along the track, look for this hit acting as a blocking path
                for (auto iter = clusterHits.begin(); iter != std::prev(clusterHits.end()); ++iter)
                {
                    const CaloHit *const pHit1{*iter}, *const pHit2{*(std::next(iter))};
                    if (LArClusterHelper::HasBlockedPath({pCaloHit}, pHit1, pHit2))
                    {
                        mergeHits.emplace_back(pCaloHit);
                        break;
                    }
                }
            }
            else
            {
                // We're outside the track, extend it
                mergeHits.emplace_back(pCaloHit);
            }
        }
    }
    else
    {
        // Entire view is missing, or has less than 3 hits, so we should merge all unmatched hits
        for (const CaloHit *const pCaloHit : unmatchedHits)
            mergeHits.emplace_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackRecoveryAlgorithm::FilterHitsToMerge(const TwoDSlidingFitResult &sfr, CaloHitList &mergeHits) const
{
    for (auto iter = mergeHits.begin(); iter != mergeHits.end();)
    {
        const CaloHit *const pCaloHit{*iter};
        const CartesianVector position{pCaloHit->GetPositionVector()};
        float rL(0.f), rT(0.f);
        sfr.GetLocalPosition(position, rL, rT);

        CartesianVector fitPosition(0, 0, 0);
        if (STATUS_CODE_SUCCESS != sfr.GetGlobalFitPosition(rL, fitPosition))
        {
            mergeHits.erase(iter++);
            continue;
        }

        float rTFit(0.f), dummy(0.f);
        sfr.GetLocalPosition(fitPosition, dummy, rTFit);

        CartesianVector fitDirection(0, 0, 0);
        if (STATUS_CODE_SUCCESS != sfr.GetGlobalFitDirection(rL, fitDirection))
        {
            mergeHits.erase(iter++);
            continue;
        }

        float dTdL(0.f);
        sfr.GetLocalDirection(fitDirection, dTdL);

        const float residual{(rT - rTFit) / std::sqrt(1.f + dTdL * dTdL)};
        if (std::abs(residual) > 0.5f)
            mergeHits.erase(iter++);
        else
            ++iter;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackRecoveryAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputClusterListNames", m_inputClusterListNames));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputPfoListName", m_inputPfoListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
