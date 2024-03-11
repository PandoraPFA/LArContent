/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterCreation/SimpleClusterCreationAlgorithm.cc
 *
 *  @brief  Implementation of the cluster creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterCreation/SimpleClusterCreationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

SimpleClusterCreationAlgorithm::SimpleClusterCreationAlgorithm() :
    m_clusteringWindowSquared(1.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SimpleClusterCreationAlgorithm::Run()
{
    const CaloHitList *pCaloHitList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));

    // Select available calo hits for clustering
    CaloHitList caloHitList;
    this->SelectCaloHits(pCaloHitList, caloHitList);

    if (caloHitList.empty())
        return STATUS_CODE_SUCCESS;

    // Build map of associations between selected calo hits
    HitAssociationMap hitAssociationMap;
    this->BuildAssociationMap(caloHitList, hitAssociationMap);

    // Create new clusters
    this->CreateClusters(caloHitList, hitAssociationMap);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SimpleClusterCreationAlgorithm::SelectCaloHits(const CaloHitList *const pInputList, CaloHitList &outputList) const
{
    for (const CaloHit *const pCaloHit : *pInputList)
    {
        if (PandoraContentApi::IsAvailable(*this, pCaloHit))
            outputList.push_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SimpleClusterCreationAlgorithm::BuildAssociationMap(const CaloHitList &caloHitList, HitAssociationMap &hitAssociationMap) const
{
    for (const CaloHit *const pCaloHitI : caloHitList)
    {
        for (const CaloHit *const pCaloHitJ : caloHitList)
        {
            if (pCaloHitI == pCaloHitJ)
                continue;

            if ((pCaloHitI->GetPositionVector() - pCaloHitJ->GetPositionVector()).GetMagnitudeSquared() < m_clusteringWindowSquared)
            {
                CaloHitList &caloHitListI(hitAssociationMap[pCaloHitI]);

                if (caloHitListI.end() == std::find(caloHitListI.begin(), caloHitListI.end(), pCaloHitJ))
                    caloHitListI.push_back(pCaloHitJ);

                CaloHitList &caloHitListJ(hitAssociationMap[pCaloHitI]);

                if (caloHitListJ.end() == std::find(caloHitListJ.begin(), caloHitListJ.end(), pCaloHitI))
                    caloHitListJ.push_back(pCaloHitI);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SimpleClusterCreationAlgorithm::CreateClusters(const CaloHitList &caloHitList, const HitAssociationMap &hitAssociationMap) const
{
    CaloHitSet vetoList;
    CaloHitVector caloHitVector(caloHitList.begin(), caloHitList.end());
    std::sort(caloHitVector.begin(), caloHitVector.end(), LArClusterHelper::SortHitsByPosition);

    for (const CaloHit *const pSeedCaloHit : caloHitVector)
    {
        if (vetoList.count(pSeedCaloHit))
            continue;

        CaloHitList mergeList;
        this->CollectAssociatedHits(pSeedCaloHit, pSeedCaloHit, hitAssociationMap, vetoList, mergeList);

        const Cluster *pCluster = NULL;
        PandoraContentApi::Cluster::Parameters parameters;
        parameters.m_caloHitList.push_back(pSeedCaloHit);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pCluster));
        vetoList.insert(pSeedCaloHit);

        for (const CaloHit *const pAssociatedCaloHit : mergeList)
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pCluster, pAssociatedCaloHit));
            vetoList.insert(pAssociatedCaloHit);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SimpleClusterCreationAlgorithm::CollectAssociatedHits(const CaloHit *const pSeedCaloHit, const CaloHit *const pCurrentCaloHit,
    const HitAssociationMap &hitAssociationMap, const CaloHitSet &vetoList, CaloHitList &mergeList) const
{
    if (vetoList.count(pCurrentCaloHit))
        return;

    HitAssociationMap::const_iterator iter1 = hitAssociationMap.find(pCurrentCaloHit);
    if (iter1 == hitAssociationMap.end())
        return;

    CaloHitVector caloHitVector(iter1->second.begin(), iter1->second.end());
    std::sort(caloHitVector.begin(), caloHitVector.end(), LArClusterHelper::SortHitsByPosition);

    for (const CaloHit *const pAssociatedCaloHit : caloHitVector)
    {
        if (pAssociatedCaloHit == pSeedCaloHit)
            continue;

        if (mergeList.end() != std::find(mergeList.begin(), mergeList.end(), pAssociatedCaloHit))
            continue;

        mergeList.push_back(pAssociatedCaloHit);

        this->CollectAssociatedHits(pSeedCaloHit, pAssociatedCaloHit, hitAssociationMap, vetoList, mergeList);
    }
}
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SimpleClusterCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    float clusteringWindow = std::sqrt(m_clusteringWindowSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ClusteringWindow", clusteringWindow));
    m_clusteringWindowSquared = clusteringWindow * clusteringWindow;

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
