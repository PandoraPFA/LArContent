/**
 *  @file   LArContent/src/LArThreeDReco/LArTransverseTrackMatching/LongTracksTool.cc
 * 
 *  @brief  Implementation of the long tracks tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArThreeDReco/LArTransverseTrackMatching/LongTracksTool.h"

using namespace pandora;

namespace lar_content
{

bool LongTracksTool::HasLongDirectConnections(IteratorList::const_iterator iIter, const IteratorList &iteratorList)
{
    for (IteratorList::const_iterator iIter2 = iteratorList.begin(), iIter2End = iteratorList.end(); iIter2 != iIter2End; ++iIter2)
    {
        if (iIter == iIter2)
            continue;

        if (((*iIter)->GetClusterU() == (*iIter2)->GetClusterU()) || ((*iIter)->GetClusterV() == (*iIter2)->GetClusterV()) || ((*iIter)->GetClusterW() == (*iIter2)->GetClusterW()) )
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LongTracksTool::IsLongerThanDirectConnections(IteratorList::const_iterator iIter, const TensorType::ElementList &elementList,
    const unsigned int minMatchedSamplingPointRatio, const pandora::ClusterList &usedClusters)
{
    const unsigned int nMatchedSamplingPoints((*iIter)->GetOverlapResult().GetNMatchedSamplingPoints());

    for (TensorType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
    {
        if ((*iIter) == eIter)
            continue;

        if (usedClusters.count(eIter->GetClusterU()) || usedClusters.count(eIter->GetClusterV()) || usedClusters.count(eIter->GetClusterW()))
            continue;

        if (((*iIter)->GetClusterU() != eIter->GetClusterU()) && ((*iIter)->GetClusterV() != eIter->GetClusterV()) && ((*iIter)->GetClusterW() != eIter->GetClusterW()))
            continue;

        if (nMatchedSamplingPoints < minMatchedSamplingPointRatio * eIter->GetOverlapResult().GetNMatchedSamplingPoints())
            return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LongTracksTool::Run(ThreeDTransverseTracksAlgorithm *pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << this->GetType() << std::endl;

    ProtoParticleVector protoParticleVector;
    this->FindLongTracks(overlapTensor, protoParticleVector);

    const bool particlesMade(pAlgorithm->CreateThreeDParticles(protoParticleVector));
    return particlesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LongTracksTool::FindLongTracks(const TensorType &overlapTensor, ProtoParticleVector &protoParticleVector) const
{
    ClusterList usedClusters;

    for (TensorType::const_iterator iterU = overlapTensor.begin(), iterUEnd = overlapTensor.end(); iterU != iterUEnd; ++iterU)
    {
        if (!iterU->first->IsAvailable())
            continue;

        unsigned int nU(0), nV(0), nW(0);
        TensorType::ElementList elementList;
        overlapTensor.GetConnectedElements(iterU->first, true, elementList, nU, nV, nW);

        IteratorList iteratorList;
        this->SelectLongElements(elementList, usedClusters, iteratorList);

        // Check that elements are not directly connected and are significantly longer than any other directly connected elements
        for (IteratorList::const_iterator iIter = iteratorList.begin(), iIterEnd = iteratorList.end(); iIter != iIterEnd; ++iIter)
        {
            if (LongTracksTool::HasLongDirectConnections(iIter, iteratorList))
                continue;

            if (!LongTracksTool::IsLongerThanDirectConnections(iIter, elementList, m_minMatchedSamplingPointRatio, usedClusters))
                continue;

            ProtoParticle protoParticle;
            protoParticle.m_clusterListU.insert((*iIter)->GetClusterU());
            protoParticle.m_clusterListV.insert((*iIter)->GetClusterV());
            protoParticle.m_clusterListW.insert((*iIter)->GetClusterW());
            protoParticleVector.push_back(protoParticle);

            usedClusters.insert((*iIter)->GetClusterU());
            usedClusters.insert((*iIter)->GetClusterV());
            usedClusters.insert((*iIter)->GetClusterW());
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LongTracksTool::SelectLongElements(const TensorType::ElementList &elementList, const pandora::ClusterList &usedClusters,
    IteratorList &iteratorList) const
{
    for (TensorType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
    {
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
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LongTracksTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_minMatchedFraction = 0.9f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedFraction", m_minMatchedFraction));

    m_minMatchedSamplingPoints = 20;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedSamplingPoints", m_minMatchedSamplingPoints));

    m_minXOverlapFraction = 0.9f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinXOverlapFraction", m_minXOverlapFraction));

    m_minMatchedSamplingPointRatio = 2;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedSamplingPointRatio", m_minMatchedSamplingPointRatio));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
