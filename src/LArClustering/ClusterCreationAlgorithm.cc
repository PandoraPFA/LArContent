/**
 *  @file   LArContent/src/LArClustering/ClusterCreationAlgorithm.cc
 * 
 *  @brief  Implementation of the cluster creation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArClustering/ClusterCreationAlgorithm.h"

#include <cassert>

using namespace pandora;

namespace lar
{

StatusCode ClusterCreationAlgorithm::Run()
{
    const CaloHitList *pCaloHitList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentCaloHitList(*this, pCaloHitList));

    OrderedCaloHitList orderedCaloHitList;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, orderedCaloHitList.Add(*pCaloHitList));

    HitAssociationMap forwardHitAssociationMap, backwardHitAssociationMap;
    this->MakePrimaryAssociations(orderedCaloHitList, forwardHitAssociationMap, backwardHitAssociationMap);
    this->MakeSecondaryAssociations(orderedCaloHitList, forwardHitAssociationMap, backwardHitAssociationMap);

    HitJoinMap hitJoinMap;
    this->IdentifyJoins(orderedCaloHitList, forwardHitAssociationMap, backwardHitAssociationMap, hitJoinMap);
    this->CreateClusters(orderedCaloHitList, hitJoinMap);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterCreationAlgorithm::MakePrimaryAssociations(const OrderedCaloHitList &orderedCaloHitList, HitAssociationMap &forwardHitAssociationMap,
    HitAssociationMap &backwardHitAssociationMap) const
{
    for (OrderedCaloHitList::const_iterator iterI = orderedCaloHitList.begin(), iterIEnd = orderedCaloHitList.end(); iterI != iterIEnd; ++iterI)
    {
        unsigned int nLayersConsidered(0);

        for (OrderedCaloHitList::const_iterator iterJ = iterI, iterJEnd = orderedCaloHitList.end(); (nLayersConsidered++ < 3) && (iterJ != iterJEnd); ++iterJ)
        {
            if (iterI->first == iterJ->first)
                continue;

            for (CaloHitList::const_iterator hitIterI = iterI->second->begin(), hitIterIEnd = iterI->second->end(); hitIterI != hitIterIEnd; ++hitIterI)
            {
                for (CaloHitList::const_iterator hitIterJ = iterJ->second->begin(), hitIterJEnd = iterJ->second->end(); hitIterJ != hitIterJEnd; ++hitIterJ)
                    this->CreatePrimaryAssociation(*hitIterI, *hitIterJ, forwardHitAssociationMap, backwardHitAssociationMap);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterCreationAlgorithm::MakeSecondaryAssociations(const OrderedCaloHitList &orderedCaloHitList, HitAssociationMap &forwardHitAssociationMap,
    HitAssociationMap &backwardHitAssociationMap) const
{
    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            CaloHit *pCaloHit = *hitIter;

            HitAssociationMap::const_iterator fwdIter = forwardHitAssociationMap.find(pCaloHit);
            CaloHit *pForwardHit((forwardHitAssociationMap.end() == fwdIter) ? NULL : fwdIter->second.GetPrimaryTarget());

            HitAssociationMap::const_iterator fwdCheckIter = backwardHitAssociationMap.find(pForwardHit);
            CaloHit *pForwardHitCheck((backwardHitAssociationMap.end() == fwdCheckIter) ? NULL : fwdCheckIter->second.GetPrimaryTarget());

            if ((NULL != pForwardHit) && (pForwardHitCheck != pCaloHit))
                this->CreateSecondaryAssociation(pCaloHit, pForwardHit, forwardHitAssociationMap, backwardHitAssociationMap);

            HitAssociationMap::const_iterator bwdIter = backwardHitAssociationMap.find(pCaloHit);
            CaloHit *pBackwardHit((backwardHitAssociationMap.end() == bwdIter) ? NULL : bwdIter->second.GetPrimaryTarget());

            HitAssociationMap::const_iterator bwdCheckIter = forwardHitAssociationMap.find(pBackwardHit);
            CaloHit *pBackwardHitCheck((forwardHitAssociationMap.end() == bwdCheckIter) ? NULL : bwdCheckIter->second.GetPrimaryTarget());

            if ((NULL != pBackwardHit) && (pBackwardHitCheck != pCaloHit))
                this->CreateSecondaryAssociation(pBackwardHit, pCaloHit, forwardHitAssociationMap, backwardHitAssociationMap);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterCreationAlgorithm::IdentifyJoins(const OrderedCaloHitList &orderedCaloHitList, const HitAssociationMap &forwardHitAssociationMap,
    const HitAssociationMap &backwardHitAssociationMap, HitJoinMap &hitJoinMap) const
{
    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            CaloHit *pCaloHit = *hitIter;
            CaloHit *pForwardJoinHit = this->GetJoinHit(pCaloHit, forwardHitAssociationMap, backwardHitAssociationMap);
            CaloHit *pBackwardJoinHit = this->GetJoinHit(pForwardJoinHit, backwardHitAssociationMap, forwardHitAssociationMap);

            if ((NULL == pForwardJoinHit) || (NULL == pBackwardJoinHit) || (pBackwardJoinHit != pCaloHit))
                continue;

            HitJoinMap::const_iterator joinIter = hitJoinMap.find(pCaloHit);

            if (hitJoinMap.end() == joinIter)
                hitJoinMap.insert(HitJoinMap::value_type(pCaloHit, pForwardJoinHit));

            if ((hitJoinMap.end() != joinIter) && (joinIter->second != pForwardJoinHit))
                throw StatusCodeException(STATUS_CODE_FAILURE);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterCreationAlgorithm::CreateClusters(const OrderedCaloHitList &orderedCaloHitList, const HitJoinMap &hitJoinMap) const
{
    HitToClusterMap hitToClusterMap;

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            CaloHit *pCaloHit = *hitIter;
            Cluster *pCluster = NULL;

            HitToClusterMap::const_iterator mapIter = hitToClusterMap.find(pCaloHit);

            if (hitToClusterMap.end() == mapIter)
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, pCaloHit, pCluster));
                hitToClusterMap.insert(HitToClusterMap::value_type(pCaloHit, pCluster));
            }
            else
            {
                pCluster = mapIter->second;
            }

            HitJoinMap::const_iterator joinIter = hitJoinMap.find(pCaloHit);

            if (hitJoinMap.end() == joinIter)
                continue;

            if (hitToClusterMap.end() != hitToClusterMap.find(joinIter->second))
                throw StatusCodeException(STATUS_CODE_FAILURE);

            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddCaloHitToCluster(*this, pCluster, joinIter->second));
            hitToClusterMap.insert(HitToClusterMap::value_type(joinIter->second, pCluster));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterCreationAlgorithm::CreatePrimaryAssociation(CaloHit *pCaloHitI, CaloHit *pCaloHitJ, HitAssociationMap &forwardHitAssociationMap,
    HitAssociationMap &backwardHitAssociationMap) const
{
    const float distanceSquared((pCaloHitJ->GetPositionVector() - pCaloHitI->GetPositionVector()).GetMagnitudeSquared());

    if (distanceSquared > HitAssociation::m_maxSeparationSquared)
        return;

    HitAssociationMap::iterator forwardIter = forwardHitAssociationMap.find(pCaloHitI);

    if (forwardHitAssociationMap.end() == forwardIter)
    {
        forwardHitAssociationMap.insert(HitAssociationMap::value_type(pCaloHitI, HitAssociation(pCaloHitJ, distanceSquared)));
    }
    else if (distanceSquared < forwardIter->second.GetPrimaryDistanceSquared())
    {
        forwardIter->second = HitAssociation(pCaloHitJ, distanceSquared);
    }

    HitAssociationMap::iterator backwardIter = backwardHitAssociationMap.find(pCaloHitJ);

    if (backwardHitAssociationMap.end() == backwardIter)
    {
        backwardHitAssociationMap.insert(HitAssociationMap::value_type(pCaloHitJ, HitAssociation(pCaloHitI, distanceSquared)));
    }
    else if (distanceSquared < backwardIter->second.GetPrimaryDistanceSquared())
    {
        backwardIter->second = HitAssociation(pCaloHitI, distanceSquared);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClusterCreationAlgorithm::CreateSecondaryAssociation(CaloHit *pCaloHitI, CaloHit *pCaloHitJ, HitAssociationMap &forwardHitAssociationMap,
    HitAssociationMap &backwardHitAssociationMap) const
{
    HitAssociationMap::iterator forwardIter = forwardHitAssociationMap.find(pCaloHitI);
    HitAssociationMap::iterator backwardIter = backwardHitAssociationMap.find(pCaloHitJ);

    if ((forwardHitAssociationMap.end() == forwardIter) || (backwardHitAssociationMap.end() == backwardIter))
        return;

    HitAssociation &forwardAssociation(forwardIter->second);
    HitAssociation &backwardAssociation(backwardIter->second);

    if ((forwardAssociation.GetPrimaryTarget() != pCaloHitJ) && (backwardAssociation.GetPrimaryTarget() == pCaloHitI))
    {
        if (backwardAssociation.GetPrimaryDistanceSquared() < forwardAssociation.GetSecondaryDistanceSquared())
            forwardAssociation.SetSecondaryTarget(pCaloHitJ, backwardAssociation.GetPrimaryDistanceSquared());
    }

    if ((backwardAssociation.GetPrimaryTarget() != pCaloHitI) && (forwardAssociation.GetPrimaryTarget() == pCaloHitJ))
    {
        if (forwardAssociation.GetPrimaryDistanceSquared() < backwardAssociation.GetSecondaryDistanceSquared())
            backwardAssociation.SetSecondaryTarget(pCaloHitI, forwardAssociation.GetPrimaryDistanceSquared());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

CaloHit *ClusterCreationAlgorithm::GetJoinHit(CaloHit *pCaloHit, const HitAssociationMap &hitAssociationMapI,
    const HitAssociationMap &hitAssociationMapJ) const
{
    HitAssociationMap::const_iterator iterI = hitAssociationMapI.find(pCaloHit);

    if (hitAssociationMapI.end() == iterI)
        return NULL;

    CaloHit *pPrimaryTarget = iterI->second.GetPrimaryTarget();
    CaloHit *pSecondaryTarget = iterI->second.GetSecondaryTarget();

    if (NULL == pSecondaryTarget)
        return pPrimaryTarget;

    unsigned int primaryNSteps(0), secondaryNSteps(0);
    CaloHit *pPrimaryTrace = this->TraceHitAssociation(pPrimaryTarget, hitAssociationMapI, hitAssociationMapJ, primaryNSteps);
    CaloHit *pSecondaryTrace = this->TraceHitAssociation(pSecondaryTarget, hitAssociationMapI, hitAssociationMapJ, secondaryNSteps);

    if ((pPrimaryTrace == pSecondaryTrace) || (secondaryNSteps < 5))
        return pPrimaryTarget;

    return NULL;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CaloHit *ClusterCreationAlgorithm::TraceHitAssociation(CaloHit *pCaloHit, const HitAssociationMap &hitAssociationMapI,
    const HitAssociationMap &hitAssociationMapJ, unsigned int &nSteps) const
{
    nSteps = 0;
    CaloHit *pThisHit = pCaloHit;
    CaloHit *pLastHit = pCaloHit;

    while (true)
    {
        ++nSteps;
        pThisHit = pLastHit;
        HitAssociationMap::const_iterator iterI = hitAssociationMapI.find(pThisHit);

        if (hitAssociationMapI.end() == iterI)
            break;

        pLastHit = iterI->second.GetPrimaryTarget();
        HitAssociationMap::const_iterator iterJ = hitAssociationMapJ.find(pLastHit);

        if (hitAssociationMapJ.end() == iterJ)
            break;

        if (iterJ->second.GetPrimaryTarget() != pThisHit)
            break;
    }

    return pThisHit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ClusterCreationAlgorithm::HitAssociation::m_maxSeparationSquared = std::numeric_limits<float>::max();
float ClusterCreationAlgorithm::HitAssociation::m_closeSeparationSquared = std::numeric_limits<float>::max();

StatusCode ClusterCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_inputCaloHitListName = "";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "InputCaloHitListName", m_inputCaloHitListName));

    m_outputClusterListName = "PrimaryClusterList";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputClusterListName", m_outputClusterListName));

    float maxSeparation = 2.7f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxSeparation", maxSeparation));
    HitAssociation::m_maxSeparationSquared = maxSeparation * maxSeparation;

    float closeSeparation = 0.9f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CloseSeparation", closeSeparation));
    HitAssociation::m_closeSeparationSquared = closeSeparation * closeSeparation;

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
