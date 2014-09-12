/**
 *  @file   LArContent/src/LArThreeDReco/LArTrackFragments/ClearTrackFragmentsTool.cc
 *
 *  @brief  Implementation of the clear track fragments tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArThreeDReco/LArTrackFragments/ClearTrackFragmentsTool.h"

#include "LArHelpers/LArClusterHelper.h"

using namespace pandora;

namespace lar_content
{

bool ClearTrackFragmentsTool::Run(ThreeDTrackFragmentsAlgorithm *pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << this->GetType() << std::endl;

    return this->FindTrackFragments(pAlgorithm, overlapTensor);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ClearTrackFragmentsTool::FindTrackFragments(ThreeDTrackFragmentsAlgorithm *pAlgorithm, const TensorType &overlapTensor) const
{
    for (TensorType::const_iterator iterU = overlapTensor.begin(), iterUEnd = overlapTensor.end(); iterU != iterUEnd; ++iterU)
    {
        if (!iterU->first->IsAvailable())
            continue;

        TensorType::ElementList elementList;

        if (!this->GetAndCheckElementList(overlapTensor, iterU->first, elementList))
            continue;

        IteratorList iteratorList;
        this->SelectClearElements(elementList, iteratorList);

        bool particlesMade(false);
        ClusterList modifiedClusters, deletedClusters, unavailableClusters, newAvailableClusters;

        for (IteratorList::const_iterator iIter = iteratorList.begin(), iIterEnd = iteratorList.end(); iIter != iIterEnd; ++iIter)
        {
            const TensorType::OverlapResult &overlapResult((*iIter)->GetOverlapResult());

            if (!this->CheckOverlapResult(overlapResult))
                continue;

            Cluster *pFragmentCluster(NULL);
            this->ProcessTensorElement(pAlgorithm, overlapResult, modifiedClusters, deletedClusters, pFragmentCluster);

            if (NULL == pFragmentCluster)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            const HitType fragmentHitType(overlapResult.GetFragmentHitType());
            Cluster *pClusterU((TPC_VIEW_U == fragmentHitType) ? pFragmentCluster : (*iIter)->GetClusterU());
            Cluster *pClusterV((TPC_VIEW_V == fragmentHitType) ? pFragmentCluster : (*iIter)->GetClusterV());
            Cluster *pClusterW((TPC_VIEW_W == fragmentHitType) ? pFragmentCluster : (*iIter)->GetClusterW());

            unavailableClusters.insert(pClusterU);
            unavailableClusters.insert(pClusterV);
            unavailableClusters.insert(pClusterW);

            if (!(pClusterU->IsAvailable() && pClusterV->IsAvailable() && pClusterW->IsAvailable()))
                throw StatusCodeException(STATUS_CODE_FAILURE);

            ProtoParticle protoParticle;
            ProtoParticleVector protoParticleVector;
            protoParticle.m_clusterListU.insert(pClusterU);
            protoParticle.m_clusterListV.insert(pClusterV);
            protoParticle.m_clusterListW.insert(pClusterW);
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

bool ClearTrackFragmentsTool::GetAndCheckElementList(const TensorType &overlapTensor, Cluster *pCluster, TensorType::ElementList &elementList) const
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
            if (!allMatchedHits.insert(*hIter).second)
                return false;
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
                    if (fragmentHits1.count(*hIter2))
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
// Cluster *pCluster = *cIter;
// std::cout << "   Fragment: " << pCluster << " " << pCluster->IsAvailable() << std::endl;
// }
// }
// --- DEBUG END ---

        if (isClearElement)
            iteratorList.push_back(eIter1);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClearTrackFragmentsTool::ProcessTensorElement(ThreeDTrackFragmentsAlgorithm *pAlgorithm, const TensorType::OverlapResult &overlapResult,
    ClusterList &modifiedClusters, ClusterList &deletedClusters, Cluster *&pFragmentCluster) const
{
    pFragmentCluster = NULL;

    const CaloHitList &caloHitList(overlapResult.GetFragmentCaloHitList());
    const ClusterList &clusterList(overlapResult.GetFragmentClusterList());
    const HitType fragmentHitType(overlapResult.GetFragmentHitType());

    const std::string currentListName((TPC_VIEW_U == fragmentHitType) ? pAlgorithm->GetClusterListNameU() :
        (TPC_VIEW_V == fragmentHitType) ? pAlgorithm->GetClusterListNameV() : pAlgorithm->GetClusterListNameW());

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*pAlgorithm, currentListName));

    for (ClusterList::const_iterator cIter = clusterList.begin(), cIterEnd = clusterList.end(); cIter != cIterEnd; ++cIter)
    {
        Cluster *pCluster = *cIter;

        if (!pCluster->IsAvailable())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        CaloHitList clusterHitList;
        pCluster->GetOrderedCaloHitList().GetCaloHitList(clusterHitList);

        CaloHitList daughterHits, separateHits;
        for (CaloHitList::const_iterator hIter = clusterHitList.begin(), hIterEnd = clusterHitList.end(); hIter != hIterEnd; ++hIter)
        {
            CaloHit *pCaloHit = *hIter;

            if (caloHitList.count(pCaloHit))
            {
                daughterHits.insert(pCaloHit);
            }
            else
            {
                separateHits.insert(pCaloHit);
            }
        }

        if (daughterHits.empty())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        modifiedClusters.insert(pCluster);
        this->Recluster(pAlgorithm, pCluster, daughterHits, separateHits, deletedClusters, pFragmentCluster);
    }

    if (NULL == pFragmentCluster)
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClearTrackFragmentsTool::Recluster(ThreeDTrackFragmentsAlgorithm *pAlgorithm, Cluster *pCluster, const CaloHitList &daughterHits,
    const CaloHitList &separateHits, ClusterList &deletedClusters, Cluster *&pFragmentCluster) const
{
    Cluster *pDaughterCluster(NULL);

    if (separateHits.empty())
    {
        pDaughterCluster = pCluster;
    }
    else
    {
        // ATTN Can't delete these clusters yet
        for (CaloHitList::const_iterator hIter = daughterHits.begin(), hIterEnd = daughterHits.end(); hIter != hIterEnd; ++hIter)
        {
            CaloHit *pCaloHit = *hIter;
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

void ClearTrackFragmentsTool::RebuildClusters(ThreeDTrackFragmentsAlgorithm *pAlgorithm, const ClusterList &modifiedClusters,
    const ClusterList &deletedClusters, ClusterList &newClusters) const
{
    ClusterList rebuildList;

    for (ClusterList::const_iterator cIter = modifiedClusters.begin(), cIterEnd = modifiedClusters.end(); cIter != cIterEnd; ++cIter)
    {
        Cluster *pCluster = *cIter;

        if (deletedClusters.count(pCluster) || !pCluster->IsAvailable())
            continue;

        rebuildList.insert(pCluster);
    }

    if (rebuildList.empty())
        return;
        
    pAlgorithm->RebuildClusters(rebuildList, newClusters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClearTrackFragmentsTool::UpdateTensor(ThreeDTrackFragmentsAlgorithm *pAlgorithm, const TensorType &overlapTensor,
    const ClusterList &unavailableClusters, const ClusterList &newAvailableClusters) const
{
    for (ClusterList::const_iterator iter = unavailableClusters.begin(), iterEnd = unavailableClusters.end(); iter != iterEnd; ++iter)
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

void ClearTrackFragmentsTool::GetAffectedKeyClusters(const TensorType &overlapTensor, const ClusterList &unavailableClusters,
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
                        affectedKeyClusters.insert(tIterU->first);

                    if (TPC_VIEW_V != fragmentHitType)
                        affectedKeyClusters.insert(tIterV->first);

                    if (TPC_VIEW_W != fragmentHitType)
                        affectedKeyClusters.insert(tIterW->first);

                    break;
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClearTrackFragmentsTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_minMatchedSamplingPointFraction = 0.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedSamplingPointFraction", m_minMatchedSamplingPointFraction));

    m_minMatchedHits = 5;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedHits", m_minMatchedHits));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
