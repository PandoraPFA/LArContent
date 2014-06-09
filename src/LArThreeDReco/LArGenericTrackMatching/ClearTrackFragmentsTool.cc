/**
 *  @file   LArContent/src/LArThreeDReco/LArGenericTrackMatching/ClearTrackFragmentsTool.cc
 *
 *  @brief  Implementation of the clear track fragments tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArThreeDReco/LArGenericTrackMatching/ClearTrackFragmentsTool.h"

using namespace pandora;

namespace lar
{

bool ClearTrackFragmentsTool::Run(ThreeDFragmentsBaseAlgorithm *pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraSettings::ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << m_algorithmToolType << std::endl;

    return this->FindTrackFragments(pAlgorithm, overlapTensor);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ClearTrackFragmentsTool::FindTrackFragments(ThreeDFragmentsBaseAlgorithm *pAlgorithm, const TensorType &overlapTensor) const
{
    for (TensorType::const_iterator iterU = overlapTensor.begin(), iterUEnd = overlapTensor.end(); iterU != iterUEnd; ++iterU)
    {
        if (!iterU->first->IsAvailable())
            continue;

        TensorType::ElementList elementList;

        if (!this->GetAndCheckElementList(overlapTensor, iterU->first, elementList))
            continue;

        if (!this->CheckForHitAmbiguities(elementList))
            continue;

        bool particlesMade(false);
        ClusterList unavailableClusters, newlyAvailableClusters;

        for (TensorType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
        {
            const TensorType::OverlapResult &overlapResult(eIter->GetOverlapResult());

            if (!this->CheckOverlapResult(overlapResult))
                continue;

            Cluster *pFragmentCluster(NULL);
            this->ProcessTensorElement(pAlgorithm,overlapResult, unavailableClusters, newlyAvailableClusters, pFragmentCluster);

            if (NULL == pFragmentCluster)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            ProtoParticle protoParticle;
            const HitType fragmentHitType(overlapResult.GetFragmentHitType());
            protoParticle.m_clusterListU.insert((TPC_VIEW_U == fragmentHitType) ? pFragmentCluster : eIter->GetClusterU());
            protoParticle.m_clusterListV.insert((TPC_VIEW_V == fragmentHitType) ? pFragmentCluster : eIter->GetClusterV());
            protoParticle.m_clusterListW.insert((TPC_VIEW_W == fragmentHitType) ? pFragmentCluster : eIter->GetClusterW());

            ProtoParticleVector protoParticleVector;
            protoParticleVector.push_back(protoParticle);
            particlesMade |= pAlgorithm->CreateThreeDParticles(protoParticleVector);
        }

        this->UpdateTensor(pAlgorithm, overlapTensor, unavailableClusters, newlyAvailableClusters);

        if (particlesMade)
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ClearTrackFragmentsTool::GetAndCheckElementList(const TensorType &overlapTensor, Cluster *pCluster, TensorType::ElementList &elementList) const
{
    unsigned int nU(0), nV(0), nW(0);
    overlapTensor.GetConnectedElements(pCluster, true, elementList, nU, nV, nW);

    const HitType fragmentHitType = ((nU * nV * nW == 1) ? (elementList.begin())->GetOverlapResult().GetFragmentHitType() :
        ((nV * nW == 1) && (nU > 1)) ? TPC_VIEW_U :
        ((nU * nW == 1) && (nV > 1)) ? TPC_VIEW_V :
        ((nU * nV == 1) && (nW > 1)) ? TPC_VIEW_W : CUSTOM);

    if (!((TPC_VIEW_U == fragmentHitType) || (TPC_VIEW_V == fragmentHitType) || (TPC_VIEW_W == fragmentHitType)))
        return false;

    for (TensorType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
    {
        if (eIter->GetOverlapResult().GetFragmentHitType() != fragmentHitType)
            return false;
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
    if (overlapResult.GetNMatchedSamplingPoints() < m_minMatchedSamplingPoints)
        return false;

    if (overlapResult.GetMatchedFraction() < m_minMatchedSamplingPointFraction)
        return false;

    if (overlapResult.GetFragmentCaloHitList().size() < m_minMatchedHits)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClearTrackFragmentsTool::ProcessTensorElement(ThreeDFragmentsBaseAlgorithm *pAlgorithm, const TensorType::OverlapResult &overlapResult,
    ClusterList &unavailableClusters, ClusterList &newlyAvailableClusters, Cluster *&pFragmentCluster) const
{
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

        unavailableClusters.insert(pCluster);
        this->Recluster(pAlgorithm, pCluster, daughterHits, separateHits, newlyAvailableClusters, pFragmentCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClearTrackFragmentsTool::Recluster(ThreeDFragmentsBaseAlgorithm *pAlgorithm, Cluster *pCluster, const CaloHitList &daughterHits,
    const CaloHitList &separateHits, ClusterList &newlyAvailableClusters, Cluster *&pFragmentCluster) const
{
    Cluster *pDaughterCluster(NULL);

    if (separateHits.empty())
    {
        pDaughterCluster = pCluster;
    }
    else
    {
        ClusterList clusterToFragment;
        clusterToFragment.insert(pCluster);
        std::string originalListName, fragmentListName;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*pAlgorithm, clusterToFragment, originalListName, fragmentListName));

        PandoraContentApi::Cluster::Parameters daughterParameters;
        daughterParameters.m_caloHitList = daughterHits;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*pAlgorithm, daughterParameters, pDaughterCluster));

        Cluster *pSeparateCluster(NULL);
        PandoraContentApi::Cluster::Parameters separateParameters;
        separateParameters.m_caloHitList = separateHits;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*pAlgorithm, separateParameters, pSeparateCluster));

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*pAlgorithm, fragmentListName, originalListName));
        newlyAvailableClusters.insert(pSeparateCluster);
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
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClearTrackFragmentsTool::UpdateTensor(ThreeDFragmentsBaseAlgorithm *pAlgorithm, const TensorType &overlapTensor,
    const ClusterList &unavailableClusters, const ClusterList &newlyAvailableClusters) const
{
    for (ClusterList::const_iterator iter = unavailableClusters.begin(), iterEnd = unavailableClusters.end(); iter != iterEnd; ++iter)
        pAlgorithm->UpdateUponDeletion(*iter);

    ClusterList newlyAvailableKeyClusters;
    pAlgorithm->SelectInputClusters(&newlyAvailableClusters, newlyAvailableKeyClusters);

    for (ClusterList::const_iterator iter = newlyAvailableKeyClusters.begin(), iterEnd = newlyAvailableKeyClusters.end(); iter != iterEnd; ++iter)
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
    m_minMatchedSamplingPoints = 10;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedSamplingPoints", m_minMatchedSamplingPoints));

    m_minMatchedSamplingPointFraction = 0.6f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedSamplingPointFraction", m_minMatchedSamplingPointFraction));

    m_minMatchedHits = 5;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedHits", m_minMatchedHits));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
