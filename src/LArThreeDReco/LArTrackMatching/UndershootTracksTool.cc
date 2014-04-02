/**
 *  @file   LArContent/src/LArThreeDReco/LArTrackMatching/UndershootTracksTool.cc
 * 
 *  @brief  Implementation of the undershoot tracks tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArPointingClusterHelper.h"
#include "LArHelpers/LArThreeDHelper.h"

#include "LArObjects/LArPointingCluster.h"

#include "LArThreeDReco/LArTrackMatching/UndershootTracksTool.h"

using namespace pandora;

namespace lar
{

bool UndershootTracksTool::Run(ThreeDTransverseTracksAlgorithm *pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraSettings::ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << m_algorithmToolType << std::endl;

    ProtoParticleVector protoParticleVector;
    this->FindUndershootTracks(overlapTensor, protoParticleVector);
    const bool changesMade(this->ApplyChanges(pAlgorithm, protoParticleVector));

    return changesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void UndershootTracksTool::FindUndershootTracks(const TensorType &overlapTensor, ProtoParticleVector &protoParticleVector) const
{
    ClusterList usedClusters;

    for (TensorType::const_iterator iterU = overlapTensor.begin(), iterUEnd = overlapTensor.end(); iterU != iterUEnd; ++iterU)
    {
        if (!iterU->first->IsAvailable())
            continue;

        unsigned int nU(0), nV(0), nW(0);
        TensorType::ElementList elementList;
        overlapTensor.GetConnectedElements(iterU->first, true, elementList, nU, nV, nW);

        if (nU * nV * nW < 2)
            continue;

        std::sort(elementList.begin(), elementList.end(), ThreeDTransverseTracksAlgorithm::SortByNMatchedSamplingPoints);

        for (TensorType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
        {
            if (!this->PassesElementCuts(eIter, usedClusters))
                continue;

            IteratorList iteratorList;
            this->SelectUndershootElements(eIter, elementList, usedClusters, iteratorList);

            if (iteratorList.size() < 2)
                continue;

            ProtoParticle protoParticle;
            this->BuildProtoParticle(iteratorList, protoParticle);

            if (protoParticle.m_clusterListU.empty() && protoParticle.m_clusterListV.empty() && protoParticle.m_clusterListW.empty())
                continue;

            protoParticleVector.push_back(protoParticle);
            usedClusters.insert(protoParticle.m_clusterListU.begin(), protoParticle.m_clusterListU.end());
            usedClusters.insert(protoParticle.m_clusterListV.begin(), protoParticle.m_clusterListV.end());
            usedClusters.insert(protoParticle.m_clusterListW.begin(), protoParticle.m_clusterListW.end());
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool UndershootTracksTool::ApplyChanges(ThreeDTransverseTracksAlgorithm *pAlgorithm, const ProtoParticleVector &protoParticleVector) const
{
    bool changesMade(false);

    for (ProtoParticleVector::const_iterator iter = protoParticleVector.begin(), iterEnd = protoParticleVector.end(); iter != iterEnd; ++iter)
    {
        changesMade |= this->MakeClusterMerges(pAlgorithm, iter->m_clusterListU);
        changesMade |= this->MakeClusterMerges(pAlgorithm, iter->m_clusterListV);
        changesMade |= this->MakeClusterMerges(pAlgorithm, iter->m_clusterListW);
    }

    return changesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool UndershootTracksTool::MakeClusterMerges(ThreeDTransverseTracksAlgorithm *pAlgorithm, const ClusterList &clusterList) const
{
    if (clusterList.size() < 2)
        return false;

    bool changesMade(false);
    Cluster *pParentCluster = *(clusterList.begin());

    const HitType hitType(LArThreeDHelper::GetClusterHitType(pParentCluster));
    const std::string clusterListName((VIEW_U == hitType) ? pAlgorithm->GetClusterListNameU() : (VIEW_V == hitType) ? pAlgorithm->GetClusterListNameV() :
        (VIEW_W == hitType) ? pAlgorithm->GetClusterListNameW() : throw StatusCodeException(STATUS_CODE_FAILURE));

    for (ClusterList::const_iterator iter = clusterList.begin(), iterEnd = clusterList.end(); iter != iterEnd; ++iter)
    {
        Cluster *pDaughterCluster = *iter;

        if (pParentCluster == pDaughterCluster)
            continue;

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*pAlgorithm, pParentCluster, pDaughterCluster, clusterListName, clusterListName));
        pAlgorithm->UpdateTensorUponMerge(pParentCluster, pDaughterCluster);
    }

    return changesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void UndershootTracksTool::SelectUndershootElements(TensorType::ElementList::const_iterator eIter, const TensorType::ElementList &elementList,
    const ClusterList &usedClusters, IteratorList &iteratorList) const
{
    iteratorList.push_back(eIter);

    for (TensorType::ElementList::const_iterator eIter2 = elementList.begin(); eIter2 != elementList.end(); ++eIter2)
    {
        if (eIter == eIter2)
            continue;

        if (!this->PassesElementCuts(eIter2, usedClusters))
            continue;

        for (IteratorList::const_iterator iIter = iteratorList.begin(); iIter != iteratorList.end(); ++iIter)
        {
            if ((*iIter) == eIter2)
                continue;

            unsigned int nMatchedClusters(0);

            if ((*iIter)->GetClusterU() == eIter2->GetClusterU())
                ++nMatchedClusters;

            if ((*iIter)->GetClusterV() == eIter2->GetClusterV())
                ++nMatchedClusters;

            if ((*iIter)->GetClusterW() == eIter2->GetClusterW())
                ++nMatchedClusters;

            if (2 == nMatchedClusters)
            {
                iteratorList.push_back(eIter2);
                return;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void UndershootTracksTool::BuildProtoParticle(const IteratorList &iteratorList, ProtoParticle &protoParticle) const
{
    for (IteratorList::const_iterator iIter1 = iteratorList.begin(), iIter1End = iteratorList.end(); iIter1 != iIter1End; ++iIter1)
    {
        for (IteratorList::const_iterator iIter2 = iIter1; iIter2 != iIter1End; ++iIter2) // TODO check this
        {
            if (iIter1 == iIter2)
                continue;

            Cluster *pCluster1(((*iIter1)->GetClusterU() != (*iIter2)->GetClusterU()) ? (*iIter1)->GetClusterU() :
                ((*iIter1)->GetClusterV() != (*iIter2)->GetClusterV()) ? (*iIter1)->GetClusterV() : (*iIter1)->GetClusterW());

            Cluster *pCluster2(((*iIter1)->GetClusterU() != (*iIter2)->GetClusterU()) ? (*iIter2)->GetClusterU() :
                ((*iIter1)->GetClusterV() != (*iIter2)->GetClusterV()) ? (*iIter2)->GetClusterV() : (*iIter2)->GetClusterW());

            if (pCluster1 == pCluster2)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            try
            {
                const LArPointingCluster pointingCluster1(pCluster1);
                const LArPointingCluster pointingCluster2(pCluster2);

                LArPointingCluster::Vertex vertex1, vertex2;
                LArPointingClusterHelper::GetClosestVertices(pointingCluster1, pointingCluster2, vertex1, vertex2);

                float transverse12(std::numeric_limits<float>::max()), transverse21(std::numeric_limits<float>::max());
                float longitudinal12(-std::numeric_limits<float>::max()), longitudinal21(-std::numeric_limits<float>::max());

                LArPointingClusterHelper::GetImpactParameters(vertex1, vertex2, longitudinal12, transverse12);
                LArPointingClusterHelper::GetImpactParameters(vertex2, vertex1, longitudinal21, transverse21);

                if (std::min(longitudinal12, longitudinal21) < m_minLongitudinalImpactParameter)
                    continue;

                if (std::min(transverse12, transverse21) > m_maxTransverseImpactParameter)
                    continue;

                const float cosTheta(-vertex1.GetDirection().GetCosOpeningAngle(vertex2.GetDirection()));

                if (cosTheta < m_minImpactParameterCosTheta)
                    continue;
            }
            catch (StatusCodeException &)
            {
                continue;
            }

            protoParticle.m_clusterListU.insert((*iIter1)->GetClusterU());
            protoParticle.m_clusterListV.insert((*iIter1)->GetClusterV());
            protoParticle.m_clusterListW.insert((*iIter1)->GetClusterW());

            protoParticle.m_clusterListU.insert((*iIter2)->GetClusterU());
            protoParticle.m_clusterListV.insert((*iIter2)->GetClusterV());
            protoParticle.m_clusterListW.insert((*iIter2)->GetClusterW());
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool UndershootTracksTool::PassesElementCuts(TensorType::ElementList::const_iterator eIter, const ClusterList &usedClusters) const
{
    if (usedClusters.count(eIter->GetClusterU()) || usedClusters.count(eIter->GetClusterV()) || usedClusters.count(eIter->GetClusterW()))
        return false;

    if (eIter->GetOverlapResult().GetMatchedFraction() < m_minMatchedFraction)
        return false;

    if (eIter->GetOverlapResult().GetNMatchedSamplingPoints() < m_minMatchedSamplingPoints)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode UndershootTracksTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_minMatchedFraction = 0.75f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedFraction", m_minMatchedFraction));

    m_minMatchedSamplingPoints = 10;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedSamplingPoints", m_minMatchedSamplingPoints));

    m_minLongitudinalImpactParameter = -1.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinLongitudinalImpactParameter", m_minLongitudinalImpactParameter));

    m_maxTransverseImpactParameter = 5.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTransverseImpactParameter", m_maxTransverseImpactParameter));

    m_minImpactParameterCosTheta = 0.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinImpactParameterCosTheta", m_minImpactParameterCosTheta));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
