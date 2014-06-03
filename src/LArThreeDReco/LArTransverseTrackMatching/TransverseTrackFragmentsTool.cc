/**
 *  @file   LArContent/src/LArThreeDReco/LArTransverseTrackMatching/TransverseTrackFragmentsTool.cc
 *
 *  @brief  Implementation of the transverse track fragments tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArThreeDReco/LArTransverseTrackMatching/TransverseTrackFragmentsTool.h"

using namespace pandora;

namespace lar
{

bool TransverseTrackFragmentsTool::Run(ThreeDTransverseTrackFragmentsAlg *pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraSettings::ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << m_algorithmToolType << std::endl;

    ProtoParticleVector protoParticleVector;
    this->FindTrackFragments(pAlgorithm, overlapTensor, protoParticleVector);

    const bool particlesMade(pAlgorithm->CreateThreeDParticles(protoParticleVector));
    return particlesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseTrackFragmentsTool::FindTrackFragments(ThreeDTransverseTrackFragmentsAlg *pAlgorithm, const TensorType &overlapTensor,
    ProtoParticleVector &protoParticleVector) const
{
    for (TensorType::const_iterator iterU = overlapTensor.begin(), iterUEnd = overlapTensor.end(); iterU != iterUEnd; ++iterU)
    {
        if (!iterU->first->IsAvailable())
            continue;

        HitType fragmentHitType(CUSTOM);
        TensorType::ElementList elementList;

        if (!this->GetElementList(overlapTensor, iterU->first, elementList, fragmentHitType))
            continue;

        if (!this->CheckForHitAmbiguities(elementList))
            continue;

        // Now reassign the hits as requested
        const std::string currentListName((TPC_VIEW_U == fragmentHitType) ? pAlgorithm->GetClusterListNameU() :
            (TPC_VIEW_V == fragmentHitType) ? pAlgorithm->GetClusterListNameV() : pAlgorithm->GetClusterListNameW());

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*pAlgorithm, currentListName));

        ClusterList unavailableClusters, newlyAvailableClusters;

        for (TensorType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
        {
std::cout << " msp " << eIter->GetOverlapResult().GetNMatchedSamplingPoints() << " nsp " << eIter->GetOverlapResult().GetNSamplingPoints() << " chi2Sum " << eIter->GetOverlapResult().GetChi2() << std::endl;
std::cout << " FragmentHitType " << eIter->GetOverlapResult().GetFragmentHitType() << std::endl;
ClusterList tempList1, tempList2, tempList3;
tempList1.insert(eIter->GetClusterU());
tempList2.insert(eIter->GetClusterV());
tempList3.insert(eIter->GetClusterW());
PANDORA_MONITORING_API(VisualizeClusters(&tempList1, "ClusterU", RED));
PANDORA_MONITORING_API(VisualizeClusters(&tempList2, "ClusterV", GREEN));
PANDORA_MONITORING_API(VisualizeClusters(&tempList3, "ClusterW", BLUE));
PANDORA_MONITORING_API(VisualizeClusters(&(eIter->GetOverlapResult().GetFragmentClusterList()), "matchedClusters", MAGENTA));
PANDORA_MONITORING_API(VisualizeCaloHits(&(eIter->GetOverlapResult().GetFragmentCaloHitList()), "matchedCaloHits", LIGHTGREEN));
PANDORA_MONITORING_API(ViewEvent());

            const CaloHitList &caloHitList(eIter->GetOverlapResult().GetFragmentCaloHitList());
            const ClusterList &clusterList(eIter->GetOverlapResult().GetFragmentClusterList());

            Cluster *pParentCluster(NULL);

            for (ClusterList::const_iterator cIter = clusterList.begin(), cIterEnd = clusterList.end(); cIter != cIterEnd; ++cIter)
            {
                Cluster *pCluster = *cIter;
                CaloHitList clusterHitList;
                pCluster->GetOrderedCaloHitList().GetCaloHitList(clusterHitList);

                CaloHitList daughterHits, separateHits;
                for (CaloHitList::const_iterator hIter = clusterHitList.begin(), hIterEnd = clusterHitList.end(); hIter != hIterEnd; ++hIter)
                {
                    CaloHit *pCaloHit = *hIter;
                    caloHitList.count(pCaloHit) ? daughterHits.insert(pCaloHit) : separateHits.insert(pCaloHit);
                }

                if (daughterHits.empty())
                    throw StatusCodeException(STATUS_CODE_FAILURE);

                if (separateHits.empty())
                {
                    unavailableClusters.insert(pCluster);

                    if (NULL == pParentCluster)
                    {
                        pParentCluster = pCluster;
                    }
                    else
                    {
                        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*pAlgorithm, pParentCluster, pCluster));
                    }
                }
                else
                {
                    ClusterList clusterToFragment;
                    clusterToFragment.insert(pCluster);
                    std::string originalListName, fragmentListName;
                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*pAlgorithm, clusterToFragment, originalListName, fragmentListName));

                    Cluster *pDaughterCluster(NULL);
                    PandoraContentApi::Cluster::Parameters daughterParameters;
                    daughterParameters.m_caloHitList = daughterHits;
                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*pAlgorithm, daughterParameters, pDaughterCluster));

                    Cluster *pSeparateCluster(NULL);
                    PandoraContentApi::Cluster::Parameters separateParameters;
                    separateParameters.m_caloHitList = separateHits;
                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*pAlgorithm, separateParameters, pSeparateCluster));

                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*pAlgorithm, fragmentListName, originalListName));

                    unavailableClusters.insert(pCluster);
                    newlyAvailableClusters.insert(pSeparateCluster);

                    if (NULL == pParentCluster)
                    {
                        pParentCluster = pDaughterCluster;
                    }
                    else
                    {
                        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*pAlgorithm, pParentCluster, pDaughterCluster));
                    }
                }
            }

            if (NULL == pParentCluster)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            ProtoParticle protoParticle;
            protoParticle.m_clusterListU.insert((TPC_VIEW_U == fragmentHitType) ? pParentCluster : eIter->GetClusterU());
            protoParticle.m_clusterListV.insert((TPC_VIEW_V == fragmentHitType) ? pParentCluster : eIter->GetClusterV());
            protoParticle.m_clusterListW.insert((TPC_VIEW_W == fragmentHitType) ? pParentCluster : eIter->GetClusterW());
            protoParticleVector.push_back(protoParticle);
std::cout << " Will make pfo " << std::endl;
PANDORA_MONITORING_API(VisualizeClusters(&protoParticle.m_clusterListU, "ClusterU", RED));
PANDORA_MONITORING_API(VisualizeClusters(&protoParticle.m_clusterListV, "ClusterV", GREEN));
PANDORA_MONITORING_API(VisualizeClusters(&protoParticle.m_clusterListW, "ClusterW", BLUE));
PANDORA_MONITORING_API(ViewEvent());
        }

        for (ClusterList::const_iterator kIter = unavailableClusters.begin(), kIterEnd = unavailableClusters.end(); kIter != kIterEnd; ++kIter)
            pAlgorithm->UpdateUponDeletion(*kIter);

        ClusterList newlyAvailableKeyClusters;
        pAlgorithm->SelectInputClusters(&newlyAvailableClusters, newlyAvailableKeyClusters);

        for (ClusterList::const_iterator kIter = newlyAvailableKeyClusters.begin(), kIterEnd = newlyAvailableKeyClusters.end(); kIter != kIterEnd; ++kIter)
            pAlgorithm->UpdateForNewCluster(*kIter);

        ClusterList keyClustersAffected;

        for (TensorType::const_iterator tIterU = overlapTensor.begin(), tIterUEnd = overlapTensor.end(); tIterU != tIterUEnd; ++tIterU)
        {
            for (TensorType::OverlapMatrix::const_iterator tIterV = tIterU->second.begin(), tIterVEnd = tIterU->second.end(); tIterV != tIterVEnd; ++tIterV)
            {
                for (TensorType::OverlapList::const_iterator tIterW = tIterV->second.begin(), tIterWEnd = tIterV->second.end(); tIterW != tIterWEnd; ++tIterW)
                {
                    const TensorType::OverlapResult &overlapResult(tIterW->second);
                    const ClusterList &fragmentClusters(overlapResult.GetFragmentClusterList());

                    for (ClusterList::const_iterator fIter = fragmentClusters.begin(), fIterEnd = fragmentClusters.end(); fIter != fIterEnd; ++fIter)
                    {
                        if (unavailableClusters.count(*fIter))
                        {
                            if (TPC_VIEW_U != fragmentHitType)
                                keyClustersAffected.insert(tIterU->first);

                            if (TPC_VIEW_V != fragmentHitType)
                                keyClustersAffected.insert(tIterV->first);

                            if (TPC_VIEW_W != fragmentHitType)
                                keyClustersAffected.insert(tIterW->first);

                            break;
                        }
                    }
                }
            }
        }

        for (ClusterList::const_iterator kIter = keyClustersAffected.begin(), kIterEnd = keyClustersAffected.end(); kIter != kIterEnd; ++kIter)
            pAlgorithm->UpdateUponDeletion(*kIter);

        for (ClusterList::const_iterator kIter = keyClustersAffected.begin(), kIterEnd = keyClustersAffected.end(); kIter != kIterEnd; ++kIter)
            pAlgorithm->UpdateForNewCluster(*kIter);

        return;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TransverseTrackFragmentsTool::GetElementList(const TensorType &overlapTensor, Cluster *pCluster, TensorType::ElementList &elementList,
    HitType &fragmentHitType) const
{
    unsigned int nU(0), nV(0), nW(0);
    overlapTensor.GetConnectedElements(pCluster, true, elementList, nU, nV, nW);

    fragmentHitType = ((nU * nV * nW == 1) ? (elementList.begin())->GetOverlapResult().GetFragmentHitType() :
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

bool TransverseTrackFragmentsTool::CheckForHitAmbiguities(const TensorType::ElementList &elementList) const
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

StatusCode TransverseTrackFragmentsTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar
