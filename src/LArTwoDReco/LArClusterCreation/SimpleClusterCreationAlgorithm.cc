/**
 *  @file   LArContent/src/LArTwoDReco/LArClusterCreation/SimpleClusterCreationAlgorithm.cc
 *
 *  @brief  Implementation of the cluster creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArTwoDReco/LArClusterCreation/SimpleClusterCreationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

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

void SimpleClusterCreationAlgorithm::SelectCaloHits(const CaloHitList* pCaloHitList, CaloHitList &caloHitList) const
{
    for (CaloHitList::const_iterator iter = pCaloHitList->begin(), iterEnd = pCaloHitList->end(); iter != iterEnd; ++iter)
    {
        CaloHit* pCaloHit = *iter;
        if (PandoraContentApi::IsAvailable(*this, pCaloHit))
            caloHitList.insert(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SimpleClusterCreationAlgorithm::BuildAssociationMap(const CaloHitList &caloHitList, HitAssociationMap &hitAssociationMap) const
{
    for (CaloHitList::const_iterator iterI = caloHitList.begin(), iterEndI = caloHitList.end(); iterI != iterEndI; ++iterI)
    {
        CaloHit* pCaloHitI = *iterI;

        for (CaloHitList::const_iterator iterJ = iterI, iterEndJ = iterEndI; iterJ != iterEndJ; ++iterJ)
        {
            CaloHit* pCaloHitJ = *iterJ;

            if (pCaloHitI == pCaloHitJ)
                continue;

            if ((pCaloHitI->GetPositionVector() - pCaloHitJ->GetPositionVector()).GetMagnitudeSquared() < m_clusteringWindowSquared)
            {
                hitAssociationMap[pCaloHitI].insert(pCaloHitJ);
                hitAssociationMap[pCaloHitJ].insert(pCaloHitI);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SimpleClusterCreationAlgorithm::CreateClusters(const CaloHitList &caloHitList, const HitAssociationMap &hitAssociationMap) const
{
    CaloHitList vetoList;

    for (CaloHitList::const_iterator iterI = caloHitList.begin(), iterEndI = caloHitList.end(); iterI != iterEndI; ++iterI)
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
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SimpleClusterCreationAlgorithm::CollectAssociatedHits(CaloHit *pSeedCaloHit, CaloHit *pCurrentCaloHit,
    const HitAssociationMap &hitAssociationMap, const CaloHitList &vetoList, CaloHitList &mergeList) const
{
    if (vetoList.count(pCurrentCaloHit))
        return;

    HitAssociationMap::const_iterator iter1 = hitAssociationMap.find(pCurrentCaloHit);
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

StatusCode SimpleClusterCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_inputCaloHitListName = "";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "InputCaloHitListName", m_inputCaloHitListName));

    m_outputClusterListName = "PrimaryClusterList";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputClusterListName", m_outputClusterListName));

    float clusteringWindow = 1.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusteringWindow", clusteringWindow));
    m_clusteringWindowSquared = clusteringWindow * clusteringWindow;

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
