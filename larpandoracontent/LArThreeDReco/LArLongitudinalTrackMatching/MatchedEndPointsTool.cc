/**
 *  @file   larpandoracontent/LArThreeDReco/LArLongitudinalTrackMatching/MatchedEndPointsTool.cc
 *
 *  @brief  Implementation of the matched end point tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArThreeDReco/LArLongitudinalTrackMatching/MatchedEndPointsTool.h"

using namespace pandora;

namespace lar_content
{

MatchedEndPointsTool::MatchedEndPointsTool() :
    m_minMatchedFraction(0.8f),
    m_maxEndPointChi2(3.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MatchedEndPointsTool::Run(ThreeViewLongitudinalTracksAlgorithm *const pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    ProtoParticleVector protoParticleVector;
    this->FindMatchedTracks(overlapTensor, protoParticleVector);

    const bool particlesMade(pAlgorithm->CreateThreeDParticles(protoParticleVector));
    return particlesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MatchedEndPointsTool::FindMatchedTracks(const TensorType &overlapTensor, ProtoParticleVector &protoParticleVector) const
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

        if (nU * nV * nW == 0)
            continue;

        std::sort(elementList.begin(), elementList.end(), MatchedEndPointsTool::SortByChiSquared);

        for (TensorType::ElementList::const_iterator iter = elementList.begin(); iter != elementList.end(); ++iter)
        {
            if (usedClusters.count(iter->GetClusterU()) || usedClusters.count(iter->GetClusterV()) || usedClusters.count(iter->GetClusterW()))
                continue;

            if (iter->GetOverlapResult().GetMatchedFraction() < m_minMatchedFraction)
                continue;

            if (std::max(iter->GetOverlapResult().GetInnerChi2(), iter->GetOverlapResult().GetOuterChi2()) > m_maxEndPointChi2)
                continue;

            ProtoParticle protoParticle;
            protoParticle.m_clusterList.push_back(iter->GetClusterU());
            protoParticle.m_clusterList.push_back(iter->GetClusterV());
            protoParticle.m_clusterList.push_back(iter->GetClusterW());
            protoParticleVector.push_back(protoParticle);

            usedClusters.insert(iter->GetClusterU());
            usedClusters.insert(iter->GetClusterV());
            usedClusters.insert(iter->GetClusterW());
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MatchedEndPointsTool::SortByChiSquared(const TensorType::Element &lhs, const TensorType::Element &rhs)
{
    return (lhs.GetOverlapResult().GetInnerChi2() + lhs.GetOverlapResult().GetOuterChi2() <
        rhs.GetOverlapResult().GetInnerChi2() + rhs.GetOverlapResult().GetOuterChi2());
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MatchedEndPointsTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinMatchedFraction", m_minMatchedFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxEndPointChi2", m_maxEndPointChi2));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
