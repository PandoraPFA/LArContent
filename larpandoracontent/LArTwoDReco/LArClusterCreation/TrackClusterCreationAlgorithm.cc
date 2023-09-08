/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterCreation/TrackClusterCreationAlgorithm.cc
 *
 *  @brief  Implementation of the cluster creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterCreation/TrackClusterCreationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

TrackClusterCreationAlgorithm::TrackClusterCreationAlgorithm() :
    m_mergeBackFilteredHits(true),
    m_maxGapLayers(2),
    m_maxCaloHitSeparationSquared(1.3f * 1.3f),
    m_minCaloHitSeparationSquared(0.4f * 0.4f),
    m_closeSeparationSquared(0.9f * 0.9f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackClusterCreationAlgorithm::Run()
{
    const CaloHitList *pCaloHitList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));

    OrderedCaloHitList selectedCaloHitList, rejectedCaloHitList;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->FilterCaloHits(pCaloHitList, selectedCaloHitList, rejectedCaloHitList));

    HitAssociationMap forwardHitAssociationMap, backwardHitAssociationMap;
    this->MakePrimaryAssociations(selectedCaloHitList, forwardHitAssociationMap, backwardHitAssociationMap);
    this->MakeSecondaryAssociations(selectedCaloHitList, forwardHitAssociationMap, backwardHitAssociationMap);

    HitJoinMap hitJoinMap;
    HitToClusterMap hitToClusterMap;
    this->IdentifyJoins(selectedCaloHitList, forwardHitAssociationMap, backwardHitAssociationMap, hitJoinMap);
    this->CreateClusters(selectedCaloHitList, hitJoinMap, hitToClusterMap);

    if (!m_mergeBackFilteredHits)
        this->CreateClusters(rejectedCaloHitList, hitJoinMap, hitToClusterMap);

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->AddFilteredCaloHits(selectedCaloHitList, rejectedCaloHitList, hitToClusterMap));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackClusterCreationAlgorithm::FilterCaloHits(
    const CaloHitList *const pCaloHitList, OrderedCaloHitList &selectedCaloHitList, OrderedCaloHitList &rejectedCaloHitList) const
{
    CaloHitList availableHitList;

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        if (PandoraContentApi::IsAvailable(*this, pCaloHit))
            availableHitList.push_back(pCaloHit);
    }

    if (availableHitList.empty())
        return STATUS_CODE_SUCCESS;

    HitType view{availableHitList.front()->GetHitType()};
    const float ratio{LArGeometryHelper::GetWirePitchRatio(this->GetPandora(), view)};
    const float minSeparationSquaredAdjusted{ratio * ratio * m_minCaloHitSeparationSquared};

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, selectedCaloHitList.Add(availableHitList));

    for (OrderedCaloHitList::const_iterator iter = selectedCaloHitList.begin(), iterEnd = selectedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        CaloHitVector caloHits(iter->second->begin(), iter->second->end());
        std::sort(caloHits.begin(), caloHits.end(), LArClusterHelper::SortHitsByPosition);

        for (const CaloHit *const pCaloHitI : caloHits)
        {
            bool useCaloHit(true);

            for (const CaloHit *const pCaloHitJ : caloHits)
            {
                if (pCaloHitI == pCaloHitJ)
                    continue;

                if ((pCaloHitI->GetMipEquivalentEnergy() < pCaloHitJ->GetMipEquivalentEnergy()) &&
                    ((pCaloHitI->GetPositionVector() - pCaloHitJ->GetPositionVector()).GetMagnitudeSquared() < minSeparationSquaredAdjusted))
                {
                    useCaloHit = false;
                    break;
                }
            }

            if (!useCaloHit)
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, rejectedCaloHitList.Add(pCaloHitI));
        }
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, selectedCaloHitList.Remove(rejectedCaloHitList));
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackClusterCreationAlgorithm::AddFilteredCaloHits(
    const OrderedCaloHitList &selectedCaloHitList, const OrderedCaloHitList &rejectedCaloHitList, HitToClusterMap &hitToClusterMap) const
{
    for (OrderedCaloHitList::const_iterator iter = rejectedCaloHitList.begin(), iterEnd = rejectedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        CaloHitList *pCaloHitList = NULL;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, selectedCaloHitList.GetCaloHitsInPseudoLayer(iter->first, pCaloHitList));

        if (!pCaloHitList || pCaloHitList->empty())
            continue;
        HitType view{pCaloHitList->front()->GetHitType()};
        const float ratio{LArGeometryHelper::GetWirePitchRatio(this->GetPandora(), view)};
        const float minSeparationSquaredAdjusted{ratio * ratio * m_minCaloHitSeparationSquared};

        CaloHitSet unavailableHits;

        CaloHitVector inputAvailableHits(iter->second->begin(), iter->second->end());
        std::sort(inputAvailableHits.begin(), inputAvailableHits.end(), LArClusterHelper::SortHitsByPosition);

        CaloHitVector clusteredHits(pCaloHitList->begin(), pCaloHitList->end());
        std::sort(clusteredHits.begin(), clusteredHits.end(), LArClusterHelper::SortHitsByPosition);

        bool carryOn(true);

        while (carryOn)
        {
            carryOn = false;
            CaloHitVector newClusteredHits;

            for (const CaloHit *const pCaloHitI : inputAvailableHits)
            {
                if (unavailableHits.count(pCaloHitI))
                    continue;

                if (hitToClusterMap.end() != hitToClusterMap.find(pCaloHitI))
                    continue;

                const CaloHit *pClosestHit = NULL;
                float closestSeparationSquared(minSeparationSquaredAdjusted);

                for (const CaloHit *const pCaloHitJ : clusteredHits)
                {
                    if (pCaloHitI->GetMipEquivalentEnergy() > pCaloHitJ->GetMipEquivalentEnergy())
                        continue;

                    const float separationSquared((pCaloHitI->GetPositionVector() - pCaloHitJ->GetPositionVector()).GetMagnitudeSquared());

                    if (separationSquared < closestSeparationSquared)
                    {
                        closestSeparationSquared = separationSquared;
                        pClosestHit = pCaloHitJ;
                    }
                }

                if (!pClosestHit)
                    continue;

                HitToClusterMap::const_iterator mapIter = hitToClusterMap.find(pClosestHit);

                if (hitToClusterMap.end() == mapIter)
                    throw StatusCodeException(STATUS_CODE_FAILURE);

                const Cluster *const pCluster = mapIter->second;
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pCluster, pCaloHitI));
                (void)hitToClusterMap.insert(HitToClusterMap::value_type(pCaloHitI, pCluster));

                newClusteredHits.push_back(pCaloHitI);
                carryOn = true;
            }

            for (const CaloHit *const pCaloHit : newClusteredHits)
            {
                clusteredHits.push_back(pCaloHit);
                unavailableHits.insert(pCaloHit);
            }
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackClusterCreationAlgorithm::MakePrimaryAssociations(const OrderedCaloHitList &orderedCaloHitList,
    HitAssociationMap &forwardHitAssociationMap, HitAssociationMap &backwardHitAssociationMap) const
{
    for (OrderedCaloHitList::const_iterator iterI = orderedCaloHitList.begin(), iterIEnd = orderedCaloHitList.end(); iterI != iterIEnd; ++iterI)
    {
        unsigned int nLayersConsidered(0);

        CaloHitVector caloHitsI(iterI->second->begin(), iterI->second->end());
        std::sort(caloHitsI.begin(), caloHitsI.end(), LArClusterHelper::SortHitsByPosition);

        for (OrderedCaloHitList::const_iterator iterJ = iterI, iterJEnd = orderedCaloHitList.end();
             (nLayersConsidered++ <= m_maxGapLayers + 1) && (iterJ != iterJEnd); ++iterJ)
        {
            if (iterJ->first == iterI->first || iterJ->first > iterI->first + m_maxGapLayers + 1)
                continue;

            CaloHitVector caloHitsJ(iterJ->second->begin(), iterJ->second->end());
            std::sort(caloHitsJ.begin(), caloHitsJ.end(), LArClusterHelper::SortHitsByPosition);

            for (const CaloHit *const pCaloHitI : caloHitsI)
            {
                for (const CaloHit *const pCaloHitJ : caloHitsJ)
                    this->CreatePrimaryAssociation(pCaloHitI, pCaloHitJ, forwardHitAssociationMap, backwardHitAssociationMap);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackClusterCreationAlgorithm::MakeSecondaryAssociations(const OrderedCaloHitList &orderedCaloHitList,
    HitAssociationMap &forwardHitAssociationMap, HitAssociationMap &backwardHitAssociationMap) const
{
    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        CaloHitVector caloHits(iter->second->begin(), iter->second->end());
        std::sort(caloHits.begin(), caloHits.end(), LArClusterHelper::SortHitsByPosition);

        for (const CaloHit *const pCaloHit : caloHits)
        {
            HitAssociationMap::const_iterator fwdIter = forwardHitAssociationMap.find(pCaloHit);
            const CaloHit *const pForwardHit((forwardHitAssociationMap.end() == fwdIter) ? NULL : fwdIter->second.GetPrimaryTarget());

            HitAssociationMap::const_iterator fwdCheckIter = backwardHitAssociationMap.find(pForwardHit);
            const CaloHit *const pForwardHitCheck((backwardHitAssociationMap.end() == fwdCheckIter) ? NULL : fwdCheckIter->second.GetPrimaryTarget());

            if ((NULL != pForwardHit) && (pForwardHitCheck != pCaloHit))
                this->CreateSecondaryAssociation(pCaloHit, pForwardHit, forwardHitAssociationMap, backwardHitAssociationMap);

            HitAssociationMap::const_iterator bwdIter = backwardHitAssociationMap.find(pCaloHit);
            const CaloHit *const pBackwardHit((backwardHitAssociationMap.end() == bwdIter) ? NULL : bwdIter->second.GetPrimaryTarget());

            HitAssociationMap::const_iterator bwdCheckIter = forwardHitAssociationMap.find(pBackwardHit);
            const CaloHit *const pBackwardHitCheck((forwardHitAssociationMap.end() == bwdCheckIter) ? NULL : bwdCheckIter->second.GetPrimaryTarget());

            if ((NULL != pBackwardHit) && (pBackwardHitCheck != pCaloHit))
                this->CreateSecondaryAssociation(pBackwardHit, pCaloHit, forwardHitAssociationMap, backwardHitAssociationMap);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackClusterCreationAlgorithm::IdentifyJoins(const OrderedCaloHitList &orderedCaloHitList,
    const HitAssociationMap &forwardHitAssociationMap, const HitAssociationMap &backwardHitAssociationMap, HitJoinMap &hitJoinMap) const
{
    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        CaloHitVector caloHits(iter->second->begin(), iter->second->end());
        std::sort(caloHits.begin(), caloHits.end(), LArClusterHelper::SortHitsByPosition);

        for (const CaloHit *const pCaloHit : caloHits)
        {
            const CaloHit *const pForwardJoinHit = this->GetJoinHit(pCaloHit, forwardHitAssociationMap, backwardHitAssociationMap);
            const CaloHit *const pBackwardJoinHit = this->GetJoinHit(pForwardJoinHit, backwardHitAssociationMap, forwardHitAssociationMap);

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

void TrackClusterCreationAlgorithm::CreateClusters(
    const OrderedCaloHitList &orderedCaloHitList, const HitJoinMap &hitJoinMap, HitToClusterMap &hitToClusterMap) const
{
    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        CaloHitVector caloHits(iter->second->begin(), iter->second->end());
        std::sort(caloHits.begin(), caloHits.end(), LArClusterHelper::SortHitsByPosition);

        for (const CaloHit *const pCaloHit : caloHits)
        {
            const Cluster *pCluster = NULL;

            HitToClusterMap::const_iterator mapIter = hitToClusterMap.find(pCaloHit);

            if (hitToClusterMap.end() == mapIter)
            {
                PandoraContentApi::Cluster::Parameters parameters;
                parameters.m_caloHitList.push_back(pCaloHit);
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pCluster));
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

            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pCluster, joinIter->second));
            hitToClusterMap.insert(HitToClusterMap::value_type(joinIter->second, pCluster));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackClusterCreationAlgorithm::CreatePrimaryAssociation(const CaloHit *const pCaloHitI, const CaloHit *const pCaloHitJ,
    HitAssociationMap &forwardHitAssociationMap, HitAssociationMap &backwardHitAssociationMap) const
{
    HitType view{pCaloHitI->GetHitType()};
    const float ratio{LArGeometryHelper::GetWirePitchRatio(this->GetPandora(), view)};
    const float maxSeparationSquaredAdjusted{ratio * ratio * m_maxCaloHitSeparationSquared};

    const float distanceSquared((pCaloHitJ->GetPositionVector() - pCaloHitI->GetPositionVector()).GetMagnitudeSquared());

    if (distanceSquared > maxSeparationSquaredAdjusted)
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

void TrackClusterCreationAlgorithm::CreateSecondaryAssociation(const CaloHit *const pCaloHitI, const CaloHit *const pCaloHitJ,
    HitAssociationMap &forwardHitAssociationMap, HitAssociationMap &backwardHitAssociationMap) const
{
    HitType view{pCaloHitI->GetHitType()};
    const float ratio{LArGeometryHelper::GetWirePitchRatio(this->GetPandora(), view)};
    const float closeSeparationSquaredAdjusted{ratio * ratio * m_closeSeparationSquared};

    HitAssociationMap::iterator forwardIter = forwardHitAssociationMap.find(pCaloHitI);
    HitAssociationMap::iterator backwardIter = backwardHitAssociationMap.find(pCaloHitJ);

    if ((forwardHitAssociationMap.end() == forwardIter) || (backwardHitAssociationMap.end() == backwardIter))
        return;

    HitAssociation &forwardAssociation(forwardIter->second);
    HitAssociation &backwardAssociation(backwardIter->second);

    if ((forwardAssociation.GetPrimaryTarget() != pCaloHitJ) && (backwardAssociation.GetPrimaryTarget() == pCaloHitI))
    {
        if ((backwardAssociation.GetPrimaryDistanceSquared() < forwardAssociation.GetSecondaryDistanceSquared()) &&
            (backwardAssociation.GetPrimaryDistanceSquared() < closeSeparationSquaredAdjusted))
        {
            forwardAssociation.SetSecondaryTarget(pCaloHitJ, backwardAssociation.GetPrimaryDistanceSquared());
        }
    }

    if ((backwardAssociation.GetPrimaryTarget() != pCaloHitI) && (forwardAssociation.GetPrimaryTarget() == pCaloHitJ))
    {
        if ((forwardAssociation.GetPrimaryDistanceSquared() < backwardAssociation.GetSecondaryDistanceSquared()) &&
            (forwardAssociation.GetPrimaryDistanceSquared() < closeSeparationSquaredAdjusted))
        {
            backwardAssociation.SetSecondaryTarget(pCaloHitI, forwardAssociation.GetPrimaryDistanceSquared());
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CaloHit *TrackClusterCreationAlgorithm::GetJoinHit(
    const CaloHit *const pCaloHit, const HitAssociationMap &hitAssociationMapI, const HitAssociationMap &hitAssociationMapJ) const
{
    HitAssociationMap::const_iterator iterI = hitAssociationMapI.find(pCaloHit);

    if (hitAssociationMapI.end() == iterI)
        return NULL;

    const CaloHit *const pPrimaryTarget = iterI->second.GetPrimaryTarget();
    const CaloHit *const pSecondaryTarget = iterI->second.GetSecondaryTarget();

    if (NULL == pSecondaryTarget)
        return pPrimaryTarget;

    unsigned int primaryNSteps(0), secondaryNSteps(0);
    const CaloHit *const pPrimaryTrace = this->TraceHitAssociation(pPrimaryTarget, hitAssociationMapI, hitAssociationMapJ, primaryNSteps);
    const CaloHit *const pSecondaryTrace = this->TraceHitAssociation(pSecondaryTarget, hitAssociationMapI, hitAssociationMapJ, secondaryNSteps);

    if ((pPrimaryTrace == pSecondaryTrace) || (secondaryNSteps < 5))
        return pPrimaryTarget;

    return NULL;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CaloHit *TrackClusterCreationAlgorithm::TraceHitAssociation(const CaloHit *const pCaloHit,
    const HitAssociationMap &hitAssociationMapI, const HitAssociationMap &hitAssociationMapJ, unsigned int &nSteps) const
{
    nSteps = 0;
    const CaloHit *pThisHit = pCaloHit;
    const CaloHit *pLastHit = pCaloHit;

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

StatusCode TrackClusterCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MergeBackFilteredHits", m_mergeBackFilteredHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxGapLayers", m_maxGapLayers));

    float maxCaloHitSeparation = std::sqrt(m_maxCaloHitSeparationSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxCaloHitSeparation", maxCaloHitSeparation));
    m_maxCaloHitSeparationSquared = maxCaloHitSeparation * maxCaloHitSeparation;

    float minCaloHitSeparation = std::sqrt(m_minCaloHitSeparationSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinCaloHitSeparation", minCaloHitSeparation));
    m_minCaloHitSeparationSquared = minCaloHitSeparation * minCaloHitSeparation;

    float closeSeparation = std::sqrt(m_closeSeparationSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CloseSeparation", closeSeparation));
    m_closeSeparationSquared = closeSeparation * closeSeparation;

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
