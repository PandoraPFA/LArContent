/**
 *  @file   larpandoracontent/LArCheating/PfoSubDominantHitRemovalAlgorithm.cc
 *
 *  @brief  Implementation of the cheating cluster creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArCheating/CheatingPfoSubDominantHitRemovalAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CheatingPfoSubDominantHitRemovalAlgorithm::CheatingPfoSubDominantHitRemovalAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingPfoSubDominantHitRemovalAlgorithm::Run()
{
    //PandoraMonitoringApi::Create(this->GetPandora());
    //PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_DEFAULT, -1.f, 1.f, 1.f);

    const PfoList *pCheatedPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_cheatedPfoListName, pCheatedPfoList));

    if (!pCheatedPfoList || pCheatedPfoList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "CheatingPfoSubDominantHitRemovalAlgorithm: unable to find pfo list " << m_cheatedPfoListName << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    for (const ParticleFlowObject *const pCheatedPfo : *pCheatedPfoList)
    {
        const MCParticle *pMatchedMCParticle(nullptr);

        this->FindMatchedMCParticle(pCheatedPfo, pMatchedMCParticle);

        if (!pMatchedMCParticle)
        {
            std::cout << "ISOBEL WHAT HAS HAPPENED!?" << std::endl;
            return STATUS_CODE_SUCCESS;
        }

        this->RemoveSubDominantHits(pCheatedPfo, pMatchedMCParticle);

        this->RemoveSpacePoints(pCheatedPfo);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingPfoSubDominantHitRemovalAlgorithm::FindMatchedMCParticle(const ParticleFlowObject *const pCheatedPfo, const MCParticle *&pMatchedMCParticle) const
{
    CaloHitList caloHitList;

    LArPfoHelper::GetCaloHits(pCheatedPfo, TPC_VIEW_U, caloHitList);
    LArPfoHelper::GetCaloHits(pCheatedPfo, TPC_VIEW_V, caloHitList);
    LArPfoHelper::GetCaloHits(pCheatedPfo, TPC_VIEW_W, caloHitList);

    MCParticleWeightMap mcParticleWeightMap;

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        for (const MCParticleWeightMap::value_type &mapEntry : pCaloHit->GetMCParticleWeightMap())
            mcParticleWeightMap[mapEntry.first] += mapEntry.second;
    }

    MCParticleList mcParticleList;

    for (const auto &mapEntry : mcParticleWeightMap)
        mcParticleList.push_back(mapEntry.first);

    mcParticleList.sort(LArMCParticleHelper::SortByMomentum);

    float bestWeight(0.f);

    for (const MCParticle *const pMCParticle : mcParticleList)
    {
        const float weight(mcParticleWeightMap.at(pMCParticle));

        if (weight > bestWeight)
        {
            pMatchedMCParticle = pMCParticle;
            bestWeight = weight;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingPfoSubDominantHitRemovalAlgorithm::RemoveSubDominantHits(const ParticleFlowObject *const pCheatedPfo, const MCParticle *const pMatchedMCParticle) const
{   
    ClusterList twoDClusterList;

    LArPfoHelper::GetClusters(pCheatedPfo, TPC_VIEW_U, twoDClusterList);
    LArPfoHelper::GetClusters(pCheatedPfo, TPC_VIEW_V, twoDClusterList);
    LArPfoHelper::GetClusters(pCheatedPfo, TPC_VIEW_W, twoDClusterList);

    for (const Cluster *const pCluster : twoDClusterList)
    {
        CaloHitList caloHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

        for (const CaloHit *const pCaloHit : caloHitList)
        {

            const MCParticle *pDominantMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));

            if (pDominantMCParticle != pMatchedMCParticle)
            {
                //CartesianVector hitPosition(pCaloHit->GetPositionVector());
                //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "REMOVED", RED, 2);

                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromCluster(*this, pCluster, pCaloHit));
            }
        }
    }

    //PandoraMonitoringApi::ViewEvent(this->GetPandora());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingPfoSubDominantHitRemovalAlgorithm::RemoveSpacePoints(const ParticleFlowObject *const pCheatedPfo) const
{   
    ClusterList threeDClusterList;
    LArPfoHelper::GetClusters(pCheatedPfo, TPC_3D, threeDClusterList);

    for (const Cluster *const pCluster : threeDClusterList)
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromPfo(*this, pCheatedPfo, pCluster));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingPfoSubDominantHitRemovalAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CheatedPfoListName", m_cheatedPfoListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
