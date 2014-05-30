/**
 *  @file   LArContent/src/LArThreeDReco/LArTransverseTrackMatching/MissingTrackTool.cc
 * 
 *  @brief  Implementation of the missing track tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "LArThreeDReco/LArTransverseTrackMatching/MissingTrackTool.h"

using namespace pandora;

namespace lar
{

bool MissingTrackTool::Run(ThreeDTransverseTracksAlgorithm *pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraSettings::ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << m_algorithmToolType << std::endl;

    ProtoParticleVector protoParticleVector;
    this->FindMissingTracks(overlapTensor, protoParticleVector);

    const bool particlesMade(pAlgorithm->CreateThreeDParticles(protoParticleVector));
    return particlesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MissingTrackTool::FindMissingTracks(const TensorType &overlapTensor, ProtoParticleVector &protoParticleVector) const
{
    ClusterList usedClusters;

    for (TensorType::const_iterator iterU = overlapTensor.begin(), iterUEnd = overlapTensor.end(); iterU != iterUEnd; ++iterU)
    {
        unsigned int nU(0), nV(0), nW(0);
        TensorType::ElementList elementList;
        overlapTensor.GetConnectedElements(iterU->first, false, elementList, nU, nV, nW);

        for (TensorType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
        {
            unsigned int nAvailable(0);

            if (eIter->GetClusterU()->IsAvailable() && !usedClusters.count(eIter->GetClusterU()))
                ++nAvailable;

            if (eIter->GetClusterV()->IsAvailable() && !usedClusters.count(eIter->GetClusterV()))
                ++nAvailable;

            if (eIter->GetClusterW()->IsAvailable() && !usedClusters.count(eIter->GetClusterW()))
                ++nAvailable;

            if (2 != nAvailable)
                continue;

            const TransverseOverlapResult &overlapResult(eIter->GetOverlapResult());

            if (overlapResult.GetNMatchedSamplingPoints() < m_minMatchedSamplingPoints)
                continue;

            if (overlapResult.GetMatchedFraction() < m_minMatchedFraction)
                continue;

            if (overlapResult.GetReducedChi2() > m_maxReducedChiSquared)
                continue;

            if ((overlapResult.GetXOverlap().GetXSpanU() < std::numeric_limits<float>::epsilon()) ||
                (overlapResult.GetXOverlap().GetXSpanV() < std::numeric_limits<float>::epsilon()) ||
                (overlapResult.GetXOverlap().GetXSpanW() < std::numeric_limits<float>::epsilon()))
            {
                continue;
            }

            const float xOverlapSpan(overlapResult.GetXOverlap().GetXOverlapSpan());

            if (eIter->GetClusterU()->IsAvailable() && (xOverlapSpan / overlapResult.GetXOverlap().GetXSpanU() < m_minXOverlapFraction))
                continue;

            if (eIter->GetClusterV()->IsAvailable() && (xOverlapSpan / overlapResult.GetXOverlap().GetXSpanV() < m_minXOverlapFraction))
                continue;

            if (eIter->GetClusterW()->IsAvailable() && (xOverlapSpan / overlapResult.GetXOverlap().GetXSpanW() < m_minXOverlapFraction))
                continue;

            ProtoParticle protoParticle;

            if (eIter->GetClusterU()->IsAvailable())
                protoParticle.m_clusterListU.insert(eIter->GetClusterU());

            if (eIter->GetClusterV()->IsAvailable())
                protoParticle.m_clusterListV.insert(eIter->GetClusterV());

            if (eIter->GetClusterW()->IsAvailable())
                protoParticle.m_clusterListW.insert(eIter->GetClusterW());

            protoParticleVector.push_back(protoParticle);
            usedClusters.insert(eIter->GetClusterU());
            usedClusters.insert(eIter->GetClusterV());
            usedClusters.insert(eIter->GetClusterW());
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MissingTrackTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_minMatchedSamplingPoints = 15;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedSamplingPoints", m_minMatchedSamplingPoints));

    m_minMatchedFraction = 0.95f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedFraction", m_minMatchedFraction));

    m_maxReducedChiSquared = 0.707f; // Based upon an agressive sigmaUVW of 0.5cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxReducedChiSquared", m_maxReducedChiSquared));

    m_minXOverlapFraction = 0.75f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinXOverlapFraction", m_minXOverlapFraction));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
