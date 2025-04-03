/**
 *  @file   larpandoracontent/LArThreeDReco/LArShowerMatching/ClearShowersTool.cc
 *
 *  @brief  Implementation of the clear showers tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArThreeDReco/LArShowerMatching/ClearShowersTool.h"

using namespace pandora;

namespace lar_content
{

ClearShowersTool::ClearShowersTool() :
    m_minMatchedFraction(0.2f),
    m_minMatchedSamplingPoints(40),
    m_minXOverlapFraction(0.5f),
    m_minMatchedSamplingPointRatio(3),
    m_minXOverlapSpanRatio(3.f),
    m_visualize(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ClearShowersTool::HasLargeDirectConnections(IteratorList::const_iterator iIter, const IteratorList &iteratorList)
{
    for (IteratorList::const_iterator iIter2 = iteratorList.begin(), iIter2End = iteratorList.end(); iIter2 != iIter2End; ++iIter2)
    {
        if (iIter == iIter2)
            continue;

        if (((*iIter)->GetClusterU() == (*iIter2)->GetClusterU()) || ((*iIter)->GetClusterV() == (*iIter2)->GetClusterV()) ||
            ((*iIter)->GetClusterW() == (*iIter2)->GetClusterW()))
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ClearShowersTool::IsLargerThanDirectConnections(IteratorList::const_iterator iIter, const TensorType::ElementList &elementList,
    const unsigned int minMatchedSamplingPointRatio, const float minXOverlapSpanRatio, const ClusterSet &usedClusters)
{
    const unsigned int nMatchedSamplingPoints((*iIter)->GetOverlapResult().GetNMatchedSamplingPoints());
    const unsigned int xOverlapSpan((*iIter)->GetOverlapResult().GetXOverlap().GetXOverlapSpan());

    for (TensorType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
    {
        if ((*iIter) == eIter)
            continue;

        if (usedClusters.count(eIter->GetClusterU()) || usedClusters.count(eIter->GetClusterV()) || usedClusters.count(eIter->GetClusterW()))
            continue;

        if (((*iIter)->GetClusterU() != eIter->GetClusterU()) && ((*iIter)->GetClusterV() != eIter->GetClusterV()) &&
            ((*iIter)->GetClusterW() != eIter->GetClusterW()))
            continue;

        if (nMatchedSamplingPoints < minMatchedSamplingPointRatio * eIter->GetOverlapResult().GetNMatchedSamplingPoints())
            return false;

        if (xOverlapSpan < minXOverlapSpanRatio * eIter->GetOverlapResult().GetXOverlap().GetXOverlapSpan())
            return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ClearShowersTool::Run(ThreeViewShowersAlgorithm *const pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    ProtoParticleVector protoParticleVector;
    this->FindClearShowers(overlapTensor, protoParticleVector);

    const bool particlesMade(pAlgorithm->CreateThreeDParticles(protoParticleVector));
    return particlesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClearShowersTool::FindClearShowers(const TensorType &overlapTensor, ProtoParticleVector &protoParticleVector) const
{
    ClusterSet usedClusters;
    ClusterVector sortedKeyClusters;
    overlapTensor.GetSortedKeyClusters(sortedKeyClusters);

    for (const Cluster *const pKeyCluster : sortedKeyClusters)
    {
        if (!pKeyCluster->IsAvailable())
            continue;

        unsigned int nU(0), nV(0), nW(0);
        TensorType::ElementList elementList;
        overlapTensor.GetConnectedElements(pKeyCluster, true, elementList, nU, nV, nW);

        IteratorList iteratorList;
        this->SelectLargeShowerElements(elementList, usedClusters, iteratorList);

        // Check that elements are not directly connected and are significantly longer than any other directly connected elements
        for (IteratorList::const_iterator iIter = iteratorList.begin(), iIterEnd = iteratorList.end(); iIter != iIterEnd; ++iIter)
        {
            if (ClearShowersTool::HasLargeDirectConnections(iIter, iteratorList))
                continue;

            if (!ClearShowersTool::IsLargerThanDirectConnections(iIter, elementList, m_minMatchedSamplingPointRatio, m_minXOverlapSpanRatio, usedClusters))
                continue;

            if (m_visualize)
            {
                ClusterList clusterList{(*iIter)->GetClusterU(), (*iIter)->GetClusterV(), (*iIter)->GetClusterW()};
                PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clusterList, "Selected", BLUE));
                PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
            }

            ProtoParticle protoParticle;
            protoParticle.m_clusterList.emplace_back((*iIter)->GetClusterU());
            protoParticle.m_clusterList.emplace_back((*iIter)->GetClusterV());
            protoParticle.m_clusterList.emplace_back((*iIter)->GetClusterW());
            protoParticleVector.emplace_back(protoParticle);

            usedClusters.insert((*iIter)->GetClusterU());
            usedClusters.insert((*iIter)->GetClusterV());
            usedClusters.insert((*iIter)->GetClusterW());
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClearShowersTool::SelectLargeShowerElements(
    const TensorType::ElementList &elementList, const pandora::ClusterSet &usedClusters, IteratorList &iteratorList) const
{
    for (TensorType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
    {
        if (m_visualize)
        {
            ClusterList clusterList{eIter->GetClusterU(), eIter->GetClusterV(), eIter->GetClusterW()};
            PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clusterList, "Considered", RED));
            PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        }

        if (usedClusters.count(eIter->GetClusterU()) || usedClusters.count(eIter->GetClusterV()) || usedClusters.count(eIter->GetClusterW()))
            continue;

        if (eIter->GetOverlapResult().GetMatchedFraction() < m_minMatchedFraction)
            continue;

        if (eIter->GetOverlapResult().GetNMatchedSamplingPoints() < m_minMatchedSamplingPoints)
            continue;

        const XOverlap &xOverlap(eIter->GetOverlapResult().GetXOverlap());

        if ((xOverlap.GetXSpanU() > std::numeric_limits<float>::epsilon()) && (xOverlap.GetXOverlapSpan() / xOverlap.GetXSpanU() > m_minXOverlapFraction) &&
            (xOverlap.GetXSpanV() > std::numeric_limits<float>::epsilon()) && (xOverlap.GetXOverlapSpan() / xOverlap.GetXSpanV() > m_minXOverlapFraction) &&
            (xOverlap.GetXSpanW() > std::numeric_limits<float>::epsilon()) && (xOverlap.GetXOverlapSpan() / xOverlap.GetXSpanW() > m_minXOverlapFraction))
        {
            iteratorList.push_back(eIter);
            if (m_visualize)
            {
                ClusterList clusterList{eIter->GetClusterU(), eIter->GetClusterV(), eIter->GetClusterW()};
                PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clusterList, "Large", GREEN));
                PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClearShowersTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinMatchedFraction", m_minMatchedFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinMatchedSamplingPoints", m_minMatchedSamplingPoints));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinXOverlapFraction", m_minXOverlapFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinMatchedSamplingPointRatio", m_minMatchedSamplingPointRatio));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinXOverlapSpanRatio", m_minXOverlapSpanRatio));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualize", m_visualize));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
