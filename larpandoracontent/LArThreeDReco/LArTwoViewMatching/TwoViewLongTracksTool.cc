/**
 *  @file   larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewLongTracksTool.cc
 *
 *  @brief  Implementation of the long tracks tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewLongTracksTool.h"

using namespace pandora;

namespace lar_content
{

TwoViewLongTracksTool::TwoViewLongTracksTool() :
    m_minMatchedFraction(0.3f), m_minMatchingScore(0.95f), m_minMatchedSamplingPoints(15), m_minXOverlapFraction(0.9f), m_minMatchedSamplingPointRatio(2)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewLongTracksTool::HasLongDirectConnections(IteratorList::const_iterator iIter, const IteratorList &iteratorList)
{
    for (IteratorList::const_iterator iIter2 = iteratorList.begin(), iIter2End = iteratorList.end(); iIter2 != iIter2End; ++iIter2)
    {
        if (iIter == iIter2)
            continue;

        if (((*iIter)->GetCluster1() == (*iIter2)->GetCluster1()) || ((*iIter)->GetCluster2() == (*iIter2)->GetCluster2()))
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewLongTracksTool::IsLongerThanDirectConnections(IteratorList::const_iterator iIter, const MatrixType::ElementList &elementList,
    const unsigned int minMatchedSamplingPointRatio, const pandora::ClusterSet &usedClusters)
{
    const unsigned int nMatchedReUpsampledSamplingPoints((*iIter)->GetOverlapResult().GetNMatchedReUpsampledSamplingPoints());

    for (MatrixType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
    {
        if ((*iIter) == eIter)
            continue;

        if (usedClusters.count(eIter->GetCluster1()) || usedClusters.count(eIter->GetCluster2()))
            continue;

        if (((*iIter)->GetCluster1() != eIter->GetCluster1()) && ((*iIter)->GetCluster2() != eIter->GetCluster2()))
            continue;

        if (nMatchedReUpsampledSamplingPoints < minMatchedSamplingPointRatio * eIter->GetOverlapResult().GetNMatchedReUpsampledSamplingPoints())
            return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewLongTracksTool::Run(TwoViewTransverseTracksAlgorithm *const pAlgorithm, MatrixType &overlapMatrix)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    ProtoParticleVector protoParticleVector;
    this->FindLongTracks(overlapMatrix, protoParticleVector);

    const bool particlesMade(pAlgorithm->CreateThreeDParticles(protoParticleVector));
    return particlesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewLongTracksTool::FindLongTracks(const MatrixType &overlapMatrix, ProtoParticleVector &protoParticleVector) const
{
    ClusterSet usedClusters;
    ClusterVector sortedKeyClusters;
    overlapMatrix.GetSortedKeyClusters(sortedKeyClusters);

    for (const Cluster *const pKeyCluster : sortedKeyClusters)
    {
        if (!pKeyCluster->IsAvailable())
            continue;

        unsigned int n0(0), n1(0);
        MatrixType::ElementList elementList;
        overlapMatrix.GetConnectedElements(pKeyCluster, true, elementList, n0, n1);

        IteratorList iteratorList;
        this->SelectLongElements(elementList, usedClusters, iteratorList);

        // Check that elements are not directly connected and are significantly longer than any other directly connected elements
        for (IteratorList::const_iterator iIter = iteratorList.begin(), iIterEnd = iteratorList.end(); iIter != iIterEnd; ++iIter)
        {
            if (TwoViewLongTracksTool::HasLongDirectConnections(iIter, iteratorList))
                continue;

            if (!TwoViewLongTracksTool::IsLongerThanDirectConnections(iIter, elementList, m_minMatchedSamplingPointRatio, usedClusters))
                continue;

            ProtoParticle protoParticle;
            protoParticle.m_clusterList.push_back((*iIter)->GetCluster1());
            protoParticle.m_clusterList.push_back((*iIter)->GetCluster2());
            protoParticleVector.push_back(protoParticle);

            usedClusters.insert((*iIter)->GetCluster1());
            usedClusters.insert((*iIter)->GetCluster2());
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewLongTracksTool::SelectLongElements(
    const MatrixType::ElementList &elementList, const pandora::ClusterSet &usedClusters, IteratorList &iteratorList) const
{
    for (MatrixType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
    {
        if (usedClusters.count(eIter->GetCluster1()) || usedClusters.count(eIter->GetCluster2()))
            continue;

        if (eIter->GetOverlapResult().GetLocallyMatchedFraction() < m_minMatchedFraction)
            continue;

        if (eIter->GetOverlapResult().GetMatchingScore() < m_minMatchingScore)
            continue;

        if (eIter->GetOverlapResult().GetNMatchedReUpsampledSamplingPoints() < m_minMatchedSamplingPoints)
            continue;

        const TwoViewXOverlap &xOverlap(eIter->GetOverlapResult().GetTwoViewXOverlap());

        if ((xOverlap.GetXOverlapFraction0() > m_minXOverlapFraction) && (xOverlap.GetXOverlapFraction1() > m_minXOverlapFraction))
        {
            iteratorList.push_back(eIter);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoViewLongTracksTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinMatchedFraction", m_minMatchedFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinMatchingScore", m_minMatchingScore));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinMatchedSamplingPoints", m_minMatchedSamplingPoints));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinXOverlapFraction", m_minXOverlapFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinMatchedSamplingPointRatio", m_minMatchedSamplingPointRatio));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
