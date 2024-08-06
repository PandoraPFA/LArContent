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

bool ClearTrackFragmentsTool::Run(ThreeViewTrackFragmentsAlgorithm *const pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    return this->FindTrackFragments(pAlgorithm, overlapTensor);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ClearTrackFragmentsTool::FindTrackFragments(ThreeViewTrackFragmentsAlgorithm *const pAlgorithm, const TensorType &overlapTensor) const
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

        if (iteratorList.empty())
            return false;

        // ATTN Cache information as will later modify tensor during reclustering. Approach relies on updating tensor
        // prior to reclustering operations and then exiting from tool (which will be scheduled to run again by algorithm)
        TensorType::ElementList::const_iterator iter(iteratorList.front());

        const TensorType::OverlapResult overlapResult(iter->GetOverlapResult());
        const HitType fragmentHitType(overlapResult.GetFragmentHitType());
        const Cluster *pClusterU(iter->GetClusterU()), *pClusterV(iter->GetClusterV()), *pClusterW(iter->GetClusterW());

        if (!this->CheckOverlapResult(overlapResult))
            continue;

        // ATTN No longer guaranteed safe to use iterator after processing tensor element, hence caching results above
        const Cluster *pFragmentCluster(nullptr);
        this->ProcessTensorElement(pAlgorithm, overlapTensor, overlapResult, pFragmentCluster);

        if (!pFragmentCluster)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        if (TPC_VIEW_U == fragmentHitType)
            pClusterU = pFragmentCluster;
        if (TPC_VIEW_V == fragmentHitType)
            pClusterV = pFragmentCluster;
        if (TPC_VIEW_W == fragmentHitType)
            pClusterW = pFragmentCluster;

        if (!(pClusterU->IsAvailable() && pClusterV->IsAvailable() && pClusterW->IsAvailable()))
            throw StatusCodeException(STATUS_CODE_FAILURE);

        // ATTN For safety, remove all clusters associated with this fragment particle from the tensor
        ClusterList fragmentClusterList, affectedKeyClusters;
        fragmentClusterList.push_back(pClusterU);
        fragmentClusterList.push_back(pClusterV);
        fragmentClusterList.push_back(pClusterW);
        this->GetAffectedKeyClusters(overlapTensor, fragmentClusterList, affectedKeyClusters);

        for (const Cluster *const pCluster : affectedKeyClusters)
            pAlgorithm->UpdateUponDeletion(pCluster);

        // Now make the particle
        ProtoParticle protoParticle;
        ProtoParticleVector protoParticleVector;
        protoParticle.m_clusterList.push_back(pClusterU);
        protoParticle.m_clusterList.push_back(pClusterV);
        protoParticle.m_clusterList.push_back(pClusterW);
        protoParticleVector.push_back(protoParticle);
        return pAlgorithm->CreateThreeDParticles(protoParticleVector);
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

bool ClearTrackFragmentsTool::CheckOverlapResult(const TensorType::OverlapResult &overlapResult) const
{
    // ATTN This method is currently mirrored in ThreeViewTrackFragmentsAlgorithm algorithm
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
        const float nCaloHits1(static_cast<float>(
            eIter1->GetClusterU()->GetNCaloHits() + eIter1->GetClusterV()->GetNCaloHits() + eIter1->GetClusterW()->GetNCaloHits()));

        bool isClearElement(true);

        for (TensorType::ElementList::const_iterator eIter2 = elementList.begin(), eIterEnd2 = elementList.end(); eIter2 != eIterEnd2; ++eIter2)
        {
            const CaloHitList &fragmentHits2(eIter2->GetOverlapResult().GetFragmentCaloHitList());
            const float nCaloHits2(static_cast<float>(
                eIter2->GetClusterU()->GetNCaloHits() + eIter2->GetClusterV()->GetNCaloHits() + eIter2->GetClusterW()->GetNCaloHits()));

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

        if (isClearElement)
            iteratorList.push_back(eIter1);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClearTrackFragmentsTool::ProcessTensorElement(ThreeViewTrackFragmentsAlgorithm *const pAlgorithm, const TensorType &overlapTensor,
    const TensorType::OverlapResult &overlapResult, const Cluster *&pFragmentCluster) const
{
    pFragmentCluster = nullptr;

    const HitType fragmentHitType(overlapResult.GetFragmentHitType());
    const std::string &currentListName(pAlgorithm->GetClusterListName(fragmentHitType));

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*pAlgorithm, currentListName));

    ClusterList fragmentClusterList(overlapResult.GetFragmentClusterList());
    fragmentClusterList.sort(LArClusterHelper::SortByNHits);
    const CaloHitSet caloHitSet(overlapResult.GetFragmentCaloHitList().begin(), overlapResult.GetFragmentCaloHitList().end());

    // Remove any clusters to be modified (or affected by modifications) from tensor
    ClusterList affectedKeyClusters;
    this->GetAffectedKeyClusters(overlapTensor, fragmentClusterList, affectedKeyClusters);

    for (const Cluster *const pCluster : affectedKeyClusters)
        pAlgorithm->UpdateUponDeletion(pCluster);

    for (const Cluster *const pCluster : fragmentClusterList)
        pAlgorithm->UpdateUponDeletion(pCluster);

    ClusterList clustersToRebuild;
    ClusterSet badClusters, deletedClusters;

    for (const Cluster *const pCluster : fragmentClusterList)
    {
        if (deletedClusters.count(pCluster))
            continue;

        if (!pCluster->IsAvailable())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        CaloHitList clusterHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(clusterHitList);

        CaloHitList daughterHits, separateHits;
        for (const CaloHit *const pCaloHit : clusterHitList)
        {
            if (caloHitSet.count(pCaloHit))
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

        this->Recluster(pAlgorithm, pCluster, daughterHits, separateHits, deletedClusters, badClusters, pFragmentCluster);

        // ATTN Fragment cluster will be used to build particle, so it shouldn't ever be bad, or deleted
        if (badClusters.count(pFragmentCluster) || deletedClusters.count(pFragmentCluster))
            throw StatusCodeException(STATUS_CODE_FAILURE);

        // ATTN Keep track of clusters to be rebuilt; does not include those for which address has been deleted at any time.
        // Note distinction between list of all deletions and the up-to-date list of bad clusters.
        // Fragment cluster will be automatically added to the output particle and never rebuilt.
        ClusterList::iterator rebuildIter(std::find(clustersToRebuild.begin(), clustersToRebuild.end(), pCluster));
        if (deletedClusters.count(pCluster))
        {
            if (clustersToRebuild.end() != rebuildIter)
                clustersToRebuild.erase(rebuildIter);
        }
        else if ((clustersToRebuild.end() == rebuildIter) && (pCluster != pFragmentCluster))
        {
            clustersToRebuild.push_back(pCluster);
        }
    }

    if (!pFragmentCluster)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    // Rebuild fragmented clusters into something better defined
    ClusterList clustersToAddToTensor;
    this->RebuildClusters(pAlgorithm, clustersToRebuild, clustersToAddToTensor);

    // ATTN Repopulate tensor according to modifications performed above
    ClusterList newKeyClusters;
    pAlgorithm->SelectInputClusters(&clustersToAddToTensor, newKeyClusters);

    for (const Cluster *const pCluster : newKeyClusters)
        pAlgorithm->UpdateForNewCluster(pCluster);

    for (const Cluster *const pCluster : affectedKeyClusters)
        pAlgorithm->UpdateForNewCluster(pCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClearTrackFragmentsTool::Recluster(ThreeViewTrackFragmentsAlgorithm *const pAlgorithm, const Cluster *const pCluster, const CaloHitList &daughterHits,
    const CaloHitList &separateHits, ClusterSet &deletedClusters, ClusterSet &badClusters, const Cluster *&pFragmentCluster) const
{
    if (separateHits.empty())
    {
        if (!pFragmentCluster)
        {
            pFragmentCluster = pCluster;
        }
        else
        {
            // ATTN Addresses can be re-used by new allocations, so can't just use list of all cluster deletions; keep track of "bad" clusters
            if (!badClusters.insert(pCluster).second)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            (void)deletedClusters.insert(pCluster);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*pAlgorithm, pFragmentCluster, pCluster));
        }
    }
    else
    {
        for (const CaloHit *const pCaloHit : daughterHits)
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromCluster(*pAlgorithm, pCluster, pCaloHit));

        if (!pFragmentCluster)
        {
            const ClusterList *pTemporaryList(nullptr);
            std::string temporaryListName, currentListName;
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*pAlgorithm, currentListName));
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                PandoraContentApi::CreateTemporaryListAndSetCurrent<ClusterList>(*pAlgorithm, pTemporaryList, temporaryListName));

            PandoraContentApi::Cluster::Parameters hitParameters;
            hitParameters.m_caloHitList = daughterHits;
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*pAlgorithm, hitParameters, pFragmentCluster));
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*pAlgorithm, temporaryListName, currentListName));
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*pAlgorithm, currentListName));

            (void)badClusters.erase(pFragmentCluster);
        }
        else
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*pAlgorithm, pFragmentCluster, &daughterHits));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClearTrackFragmentsTool::RebuildClusters(
    ThreeViewTrackFragmentsAlgorithm *const pAlgorithm, const ClusterList &modifiedClusters, ClusterList &newClusters) const
{
    ClusterList rebuildList;

    for (const Cluster *const pCluster : modifiedClusters)
    {
        if (pCluster->IsAvailable())
            rebuildList.push_back(pCluster);
    }

    if (!rebuildList.empty())
        pAlgorithm->RebuildClusters(rebuildList, newClusters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClearTrackFragmentsTool::GetAffectedKeyClusters(
    const TensorType &overlapTensor, const ClusterList &clustersToRemoveFromTensor, ClusterList &affectedKeyClusters) const
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
                    if (clustersToRemoveFromTensor.end() == std::find(clustersToRemoveFromTensor.begin(), clustersToRemoveFromTensor.end(), *fIter))
                        continue;

                    if ((TPC_VIEW_U != fragmentHitType) &&
                        (affectedKeyClusters.end() == std::find(affectedKeyClusters.begin(), affectedKeyClusters.end(), tIterU->first)))
                        affectedKeyClusters.push_back(tIterU->first);

                    if ((TPC_VIEW_V != fragmentHitType) &&
                        (affectedKeyClusters.end() == std::find(affectedKeyClusters.begin(), affectedKeyClusters.end(), tIterV->first)))
                        affectedKeyClusters.push_back(tIterV->first);

                    if ((TPC_VIEW_W != fragmentHitType) &&
                        (affectedKeyClusters.end() == std::find(affectedKeyClusters.begin(), affectedKeyClusters.end(), tIterW->first)))
                        affectedKeyClusters.push_back(tIterW->first);

                    break;
                }
            }
        }
    }

    affectedKeyClusters.sort(LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClearTrackFragmentsTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinMatchedSamplingPointFraction", m_minMatchedSamplingPointFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinMatchedHits", m_minMatchedHits));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
