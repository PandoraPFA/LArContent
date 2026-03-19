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
                    this->FitAndOrderCluster(pCluster, clusterToHitMap[pCluster]);
                else
                    LArClusterHelper::GetAllHits(pCluster, clusterToHitMap[pCluster]);
                for (const CaloHit *const pCaloHit : clusterToHitMap[pCluster])
                    hitToClusterMap[pCaloHit] = pCluster;
            }
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackRecoveryAlgorithm::FitAndOrderCluster(const Cluster *const pCluster, CaloHitList &orderedHits) const
{
    if (pCluster->GetNCaloHits() < 3)
    {
        LArClusterHelper::GetAllHits(pCluster, orderedHits);
        return;
    }
    const HitType view{LArClusterHelper::GetClusterHitType(pCluster)};
    TwoDSlidingFitResult sfr(pCluster, 3, LArGeometryHelper::GetWirePitch(this->GetPandora(), view));
    LArClusterHelper::OrderHitsAlongTrajectory(pCluster, sfr, orderedHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackRecoveryAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputClusterListNames", m_inputClusterListNames));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputPfoListName", m_inputPfoListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
