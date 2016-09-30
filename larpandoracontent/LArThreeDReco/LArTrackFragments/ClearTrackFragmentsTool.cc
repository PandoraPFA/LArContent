/**
 *  @file   larpandoracontent/LArThreeDReco/LArTrackFragments/ClearTrackFragmentsTool.cc
 *
 *  @brief  Implementation of the clear track fragments tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArThreeDReco/LArTrackFragments/ClearTrackFragmentsTool.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

using namespace pandora;

namespace lar_content
{

ClearTrackFragmentsTool::ClearTrackFragmentsTool() :
    m_minMatchedSamplingPointFraction(0.5f),
    m_minMatchedHits(5)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ClearTrackFragmentsTool::Run(ThreeDTrackFragmentsAlgorithm *const pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << this->GetType() << std::endl;

    return this->FindTrackFragments(pAlgorithm, overlapTensor);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ClearTrackFragmentsTool::FindTrackFragments(ThreeDTrackFragmentsAlgorithm *const pAlgorithm, const TensorType &overlapTensor) const
{
    ClusterVector sortedKeyClusters;
    overlapTensor.GetSortedKeyClusters(sortedKeyClusters);

    for (const Cluster *const pKeyCluster : sortedKeyClusters)
    {
        if (!pKeyCluster->IsAvailable())
            continue;

        TensorType::ElementList elementList;

        if (!this->GetAndCheckElementList(overlapTensor, pKeyCluster, elementList))
            continue;

        IteratorList iteratorList;
        this->SelectClearElements(elementList, iteratorList);

        bool particlesMade(false);
        ClusterSet deletedClusters, unavailableClusters;
        ClusterList modifiedClusters, newAvailableClusters;

        for (IteratorList::const_iterator iIter = iteratorList.begin(), iIterEnd = iteratorList.end(); iIter != iIterEnd; ++iIter)
        {
            const TensorType::OverlapResult &overlapResult((*iIter)->GetOverlapResult());

            if (!this->CheckOverlapResult(overlapResult))
                continue;

            const Cluster *pFragmentCluster(NULL);
            this->ProcessTensorElement(pAlgorithm, overlapResult, modifiedClusters, deletedClusters, pFragmentCluster);

            if (NULL == pFragmentCluster)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            const HitType fragmentHitType(overlapResult.GetFragmentHitType());
            const Cluster *const pClusterU((TPC_VIEW_U == fragmentHitType) ? pFragmentCluster : (*iIter)->GetClusterU());
            const Cluster *const pClusterV((TPC_VIEW_V == fragmentHitType) ? pFragmentCluster : (*iIter)->GetClusterV());
            const Cluster *const pClusterW((TPC_VIEW_W == fragmentHitType) ? pFragmentCluster : (*iIter)->GetClusterW());

            unavailableClusters.insert(pClusterU);
            unavailableClusters.insert(pClusterV);
            unavailableClusters.insert(pClusterW);

            if (!(pClusterU->IsAvailable() && pClusterV->IsAvailable() && pClusterW->IsAvailable()))
                throw StatusCodeException(STATUS_CODE_FAILURE);

            ProtoParticle protoParticle;
            ProtoParticleVector protoParticleVector;
            protoParticle.m_clusterListU.push_back(pClusterU);
            protoParticle.m_clusterListV.push_back(pClusterV);
            protoParticle.m_clusterListW.push_back(pClusterW);
            protoParticleVector.push_back(protoParticle);
            particlesMade |= pAlgorithm->CreateThreeDParticles(protoParticleVector);
        }

        unavailableClusters.insert(modifiedClusters.begin(), modifiedClusters.end());

        this->RebuildClusters(pAlgorithm, modifiedClusters, deletedClusters, newAvailableClusters);
        this->UpdateTensor(pAlgorithm, overlapTensor, unavailableClusters, newAvailableClusters);

        if (particlesMade)
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ClearTrackFragmentsTool::GetAndCheckElementList(const TensorType &overlapTensor, const Cluster *const pCluster, TensorType::ElementList &elementList) const
{
    // Get list of connected elements from tensor
    unsigned int nU(0), nV(0), nW(0);
    overlapTensor.GetConnectedElements(pCluster, true, elementList, nU, nV, nW);

    // Only allow one fragment hit type
    HitType fragmentHitType(HIT_CUSTOM);

    for (TensorType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
    {
        const HitType thisHitType(eIter->GetOverlapResult().GetFragmentHitType());

        if (!((TPC_VIEW_U == thisHitType) || (TPC_VIEW_V == thisHitType) || (TPC_VIEW_W == thisHitType)))
            throw StatusCodeException(STATUS_CODE_FAILURE);

        if (thisHitType != fragmentHitType && HIT_CUSTOM != fragmentHitType)
            return false;

        fragmentHitType = thisHitType;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ClearTrackFragmentsTool::CheckForHitAmbiguities(const TensorType::ElementList &elementList) const
{
    CaloHitList allMatchedHits;

    for (TensorType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
    {
        const CaloHitList &fragmentHits(eIter->GetOverlapResult().GetFragmentCaloHitList());

        for (CaloHitList::const_iterator hIter = fragmentHits.begin(), hIterEnd = fragmentHits.end(); hIter != hIterEnd; ++hIter)
        {
            if (allMatchedHits.end() != std::find(allMatchedHits.begin(), allMatchedHits.end(), *hIter))
                return false;

            allMatchedHits.push_back(*hIter);
        }
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ClearTrackFragmentsTool::CheckOverlapResult(const TensorType::OverlapResult &overlapResult) const
{
    // ATTN This method is currently mirrored in ThreeDTrackFragmentsAlgorithm algorithm

    if (overlapResult.GetMatchedFraction() < m_minMatchedSamplingPointFraction)
        return false;

    if (overlapResult.GetFragmentCaloHitList().size() < m_minMatchedHits)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClearTrackFragmentsTool::SelectClearElements(const TensorType::ElementList &elementList, IteratorList &iteratorList) const
{
    for (TensorType::ElementList::const_iterator eIter1 = elementList.begin(), eIterEnd1 = elementList.end(); eIter1 != eIterEnd1; ++eIter1)
    {
        const CaloHitList &fragmentHits1(eIter1->GetOverlapResult().GetFragmentCaloHitList());
        const float nCaloHits1(static_cast<float>(eIter1->GetClusterU()->GetNCaloHits() + eIter1->GetClusterV()->GetNCaloHits() +
            eIter1->GetClusterW()->GetNCaloHits()));

        bool isClearElement(true);

        for (TensorType::ElementList::const_iterator eIter2 = elementList.begin(), eIterEnd2 = elementList.end(); eIter2 != eIterEnd2; ++eIter2)
        {
            const CaloHitList &fragmentHits2(eIter2->GetOverlapResult().GetFragmentCaloHitList());
            const float nCaloHits2(static_cast<float>(eIter2->GetClusterU()->GetNCaloHits() + eIter2->GetClusterV()->GetNCaloHits() +
                    eIter2->GetClusterW()->GetNCaloHits()));

            const bool commonClusterU(eIter1->GetClusterU() == eIter2->GetClusterU());
            const bool commonClusterV(eIter1->GetClusterV() == eIter2->GetClusterV());
            const bool commonClusterW(eIter1->GetClusterW() == eIter2->GetClusterW());

            if (commonClusterU && commonClusterV && commonClusterW)
                continue;

            if (eIter1->GetOverlapResult().GetFragmentHitType() != eIter2->GetOverlapResult().GetFragmentHitType())
                throw StatusCodeException(STATUS_CODE_FAILURE);

            bool isAmbiguousElement(commonClusterU || commonClusterV || commonClusterW);

            if (!isAmbiguousElement)
            {
                for (CaloHitList::const_iterator hIter2 = fragmentHits2.begin(), hIterEnd2 = fragmentHits2.end(); hIter2 != hIterEnd2; ++hIter2)
                {
                    if (fragmentHits1.end() != std::find(fragmentHits1.begin(), fragmentHits1.end(), *hIter2))
                    {
                        isAmbiguousElement = true;
                        break;
                    }
                }
            }

            if (isAmbiguousElement && nCaloHits2 > 0.25f * nCaloHits1)
            {
                isClearElement = false;
                break;
            }
        }

// --- DEBUG BEGIN ---
// if(isClearElement)
// {
// std::cout << " pU=" << eIter1->GetClusterU() << " pV=" << eIter1->GetClusterV() << " pW=" << eIter1->GetClusterW() << "  (" <<  eIter1->GetClusterU()->IsAvailable() << "," << eIter1->GetClusterV()->IsAvailable() << "," << eIter1->GetClusterW()->IsAvailable() << ")" << std::endl;
// const ClusterList &clusterList(eIter1->GetOverlapResult().GetFragmentClusterList());
// for (ClusterList::const_iterator cIter = clusterList.begin(), cIterEnd = clusterList.end(); cIter != cIterEnd; ++cIter)
// {
// Cluster *const pCluster = *cIter;
// std::cout << "   Fragment: " << pCluster << " " << pCluster->IsAvailable() << std::endl;
// }
// }
// --- DEBUG END ---

        if (isClearElement)
            iteratorList.push_back(eIter1);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClearTrackFragmentsTool::ProcessTensorElement(ThreeDTrackFragmentsAlgorithm *const pAlgorithm, const TensorType::OverlapResult &overlapResult,
    ClusterList &modifiedClusters, ClusterSet &deletedClusters, const Cluster *&pFragmentCluster) const
{
    pFragmentCluster = NULL;

    const HitType fragmentHitType(overlapResult.GetFragmentHitType());
    const CaloHitList &caloHitList(overlapResult.GetFragmentCaloHitList());
    ClusterVector clusterVector(overlapResult.GetFragmentClusterList().begin(), overlapResult.GetFragmentClusterList().end());
    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);

    const std::string currentListName((TPC_VIEW_U == fragmentHitType) ? pAlgorithm->GetClusterListNameU() :
        (TPC_VIEW_V == fragmentHitType) ? pAlgorithm->GetClusterListNameV() : pAlgorithm->GetClusterListNameW());

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*pAlgorithm, currentListName));

    for (const Cluster *const pCluster : clusterVector)
    {
        if (!pCluster->IsAvailable())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        CaloHitList clusterHitList;
        pCluster->GetOrderedCaloHitList().GetCaloHitList(clusterHitList);

        CaloHitList daughterHits, separateHits;
        for (const CaloHit *const pCaloHit : clusterHitList)
        {
            if (caloHitList.end() != std::find(caloHitList.begin(), caloHitList.end(), pCaloHit))
            {
                daughterHits.push_back(pCaloHit);
            }
            else
            {
                separateHits.push_back(pCaloHit);
            }
        }

        if (daughterHits.empty())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        modifiedClusters.push_back(pCluster);
        this->Recluster(pAlgorithm, pCluster, daughterHits, separateHits, deletedClusters, pFragmentCluster);
    }

    if (NULL == pFragmentCluster)
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClearTrackFragmentsTool::Recluster(ThreeDTrackFragmentsAlgorithm *const pAlgorithm, const Cluster *const pCluster, const CaloHitList &daughterHits,
    const CaloHitList &separateHits, ClusterSet &deletedClusters, const Cluster *&pFragmentCluster) const
{
    const Cluster *pDaughterCluster(NULL);

    if (separateHits.empty())
    {
        pDaughterCluster = pCluster;
    }
    else
    {
        // ATTN Can't delete these clusters yet
        for (CaloHitList::const_iterator hIter = daughterHits.begin(), hIterEnd = daughterHits.end(); hIter != hIterEnd; ++hIter)
        {
            const CaloHit *const pCaloHit = *hIter;
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromCluster(*pAlgorithm, pCluster, pCaloHit));
        }

        const ClusterList *pTemporaryList = NULL;
        std::string temporaryListName, currentListName;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*pAlgorithm, currentListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent<ClusterList>(*pAlgorithm,
            pTemporaryList, temporaryListName));

        PandoraContentApi::Cluster::Parameters hitParameters;
        hitParameters.m_caloHitList = daughterHits;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*pAlgorithm, hitParameters, pDaughterCluster));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*pAlgorithm, temporaryListName, currentListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*pAlgorithm, currentListName));
    }

    if (NULL == pDaughterCluster)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    if (NULL == pFragmentCluster)
    {
        pFragmentCluster = pDaughterCluster;
    }
    else
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*pAlgorithm, pFragmentCluster, pDaughterCluster));
        deletedClusters.insert(pDaughterCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClearTrackFragmentsTool::RebuildClusters(ThreeDTrackFragmentsAlgorithm *const pAlgorithm, const ClusterList &modifiedClusters,
    const ClusterSet &deletedClusters, ClusterList &newClusters) const
{
    ClusterList rebuildList;

    for (ClusterList::const_iterator cIter = modifiedClusters.begin(), cIterEnd = modifiedClusters.end(); cIter != cIterEnd; ++cIter)
    {
        const Cluster *const pCluster = *cIter;

        if (!pCluster->IsAvailable() || (deletedClusters.count(pCluster)))
            continue;

        rebuildList.push_back(pCluster);
    }

    if (rebuildList.empty())
        return;

    pAlgorithm->RebuildClusters(rebuildList, newClusters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClearTrackFragmentsTool::UpdateTensor(ThreeDTrackFragmentsAlgorithm *const pAlgorithm, const TensorType &overlapTensor,
    const ClusterSet &unavailableClusters, const ClusterList &newAvailableClusters) const
{
    // TODO might need these to be handled in a well-defined order (but may contain dangling pointers)
    for (ClusterSet::const_iterator iter = unavailableClusters.begin(), iterEnd = unavailableClusters.end(); iter != iterEnd; ++iter)
        pAlgorithm->UpdateUponDeletion(*iter);

    ClusterList newAvailableKeyClusters;
    pAlgorithm->SelectInputClusters(&newAvailableClusters, newAvailableKeyClusters);

    for (ClusterList::const_iterator iter = newAvailableKeyClusters.begin(), iterEnd = newAvailableKeyClusters.end(); iter != iterEnd; ++iter)
        pAlgorithm->UpdateForNewCluster(*iter);

    ClusterList affectedKeyClusters;
    this->GetAffectedKeyClusters(overlapTensor, unavailableClusters, affectedKeyClusters);

    for (ClusterList::const_iterator iter = affectedKeyClusters.begin(), iterEnd = affectedKeyClusters.end(); iter != iterEnd; ++iter)
        pAlgorithm->UpdateUponDeletion(*iter);

    for (ClusterList::const_iterator iter = affectedKeyClusters.begin(), iterEnd = affectedKeyClusters.end(); iter != iterEnd; ++iter)
        pAlgorithm->UpdateForNewCluster(*iter);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClearTrackFragmentsTool::GetAffectedKeyClusters(const TensorType &overlapTensor, const ClusterSet &unavailableClusters,
    ClusterList &affectedKeyClusters) const
{
    for (TensorType::const_iterator tIterU = overlapTensor.begin(), tIterUEnd = overlapTensor.end(); tIterU != tIterUEnd; ++tIterU)
    {
        for (TensorType::OverlapMatrix::const_iterator tIterV = tIterU->second.begin(), tIterVEnd = tIterU->second.end(); tIterV != tIterVEnd; ++tIterV)
        {
            for (TensorType::OverlapList::const_iterator tIterW = tIterV->second.begin(), tIterWEnd = tIterV->second.end(); tIterW != tIterWEnd; ++tIterW)
            {
                const TensorType::OverlapResult &overlapResult(tIterW->second);
                const HitType fragmentHitType(overlapResult.GetFragmentHitType());
                const ClusterList &fragmentClusters(overlapResult.GetFragmentClusterList());

                for (ClusterList::const_iterator fIter = fragmentClusters.begin(), fIterEnd = fragmentClusters.end(); fIter != fIterEnd; ++fIter)
                {
                    if (!unavailableClusters.count(*fIter))
                        continue;

                    if (TPC_VIEW_U != fragmentHitType)
                        affectedKeyClusters.push_back(tIterU->first);

                    if (TPC_VIEW_V != fragmentHitType)
                        affectedKeyClusters.push_back(tIterV->first);

                    if (TPC_VIEW_W != fragmentHitType)
                        affectedKeyClusters.push_back(tIterW->first);

                    break;
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClearTrackFragmentsTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedSamplingPointFraction", m_minMatchedSamplingPointFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedHits", m_minMatchedHits));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
