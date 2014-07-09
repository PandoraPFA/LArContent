/**
 *  @file   LArContent/src/LArThreeDReco/LArTransverseTrackMatching/ClearTracksTool.cc
 * 
 *  @brief  Implementation of the clear tracks tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "LArThreeDReco/LArTransverseTrackMatching/ClearTracksTool.h"

using namespace pandora;

namespace lar
{

bool ClearTracksTool::Run(ThreeDTransverseTracksAlgorithm *pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraSettings::ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << m_algorithmToolType << std::endl;

    bool particlesMade(false);

    TensorType::ElementList elementList;
    overlapTensor.GetUnambiguousElements(true, elementList);
    this->CreateThreeDParticles(pAlgorithm, elementList, particlesMade);

    elementList.clear();
    overlapTensor.GetUnambiguousElements(true, &TrackTrackTrackAmbiguity, elementList);
    this->CreateThreeDParticles(pAlgorithm, elementList, particlesMade);

    elementList.clear();
    overlapTensor.GetUnambiguousElements(true, &TrackTrackShowerAmbiguity, elementList);
    this->CreateThreeDParticles(pAlgorithm, elementList, particlesMade);

    elementList.clear();
    overlapTensor.GetUnambiguousElements(true, &TrackShowerShowerAmbiguity, elementList);
    this->CreateThreeDParticles(pAlgorithm, elementList, particlesMade);

    elementList.clear();
    overlapTensor.GetUnambiguousElements(true, &ShowerShowerShowerAmbiguity, elementList);
    this->CreateThreeDParticles(pAlgorithm, elementList, particlesMade);

    return particlesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClearTracksTool::CreateThreeDParticles(ThreeDTransverseTracksAlgorithm *pAlgorithm, const TensorType::ElementList &elementList,
    bool &particlesMade) const
{
    ProtoParticleVector protoParticleVector;

    for (TensorType::ElementList::const_iterator iter = elementList.begin(), iterEnd = elementList.end(); iter != iterEnd; ++iter)
    {
        if (iter->GetOverlapResult().GetMatchedFraction() < m_minMatchedFraction)
            continue;

        const TransverseOverlapResult::XOverlap &xOverlap(iter->GetOverlapResult().GetXOverlap());

        if ((xOverlap.GetXSpanU() < std::numeric_limits<float>::epsilon()) || (xOverlap.GetXOverlapSpan() / xOverlap.GetXSpanU() < m_minXOverlapFraction))
            continue;

        if ((xOverlap.GetXSpanV() < std::numeric_limits<float>::epsilon()) || (xOverlap.GetXOverlapSpan() / xOverlap.GetXSpanV() < m_minXOverlapFraction))
            continue;

        if ((xOverlap.GetXSpanW() < std::numeric_limits<float>::epsilon()) || (xOverlap.GetXOverlapSpan() / xOverlap.GetXSpanW() < m_minXOverlapFraction))
            continue;

        ProtoParticle protoParticle;
        protoParticle.m_clusterListU.insert(iter->GetClusterU());
        protoParticle.m_clusterListV.insert(iter->GetClusterV());
        protoParticle.m_clusterListW.insert(iter->GetClusterW());
        protoParticleVector.push_back(protoParticle);
    }

    particlesMade |= pAlgorithm->CreateThreeDParticles(protoParticleVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClearTracksTool::ClassifyClusters(const ClusterList &clusterList, ClusterList &trackClusterList, ClusterList &showerClusterList)
{
    for (ClusterList::const_iterator iter = clusterList.begin(), iterEnd = clusterList.end(); iter != iterEnd; ++iter)
    {
        if ((*iter)->IsMipTrack())
        {
            trackClusterList.insert(*iter);
        }
        else
        {
            showerClusterList.insert(*iter);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ClearTracksTool::TrackTrackTrackAmbiguity(const ClusterList &clusterListU, const ClusterList &clusterListV, const ClusterList &clusterListW,
    Cluster *&pClusterU, Cluster *&pClusterV, Cluster *&pClusterW)
{
    ClusterList trackClustersU, showerClustersU;
    ClearTracksTool::ClassifyClusters(clusterListU, trackClustersU, showerClustersU);

    if (1 != trackClustersU.size())
        return false;

    ClusterList trackClustersV, showerClustersV;
    ClearTracksTool::ClassifyClusters(clusterListV, trackClustersV, showerClustersV);

    if (1 != trackClustersV.size())
        return false;

    ClusterList trackClustersW, showerClustersW;
    ClearTracksTool::ClassifyClusters(clusterListW, trackClustersW, showerClustersW);

    if (1 != trackClustersW.size())
        return false;

    pClusterU = *(trackClustersU.begin());
    pClusterV = *(trackClustersV.begin());
    pClusterW = *(trackClustersW.begin());

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ClearTracksTool::TrackTrackShowerAmbiguity(const ClusterList &clusterListU, const ClusterList &clusterListV, const ClusterList &clusterListW,
    Cluster *&pClusterU, Cluster *&pClusterV, Cluster *&pClusterW)
{
    ClusterList trackClustersU, showerClustersU;
    ClearTracksTool::ClassifyClusters(clusterListU, trackClustersU, showerClustersU);

    ClusterList trackClustersV, showerClustersV;
    ClearTracksTool::ClassifyClusters(clusterListV, trackClustersV, showerClustersV);

    ClusterList trackClustersW, showerClustersW;
    ClearTracksTool::ClassifyClusters(clusterListW, trackClustersW, showerClustersW);

    const unsigned int nTracksU(trackClustersU.size()), nTracksV(trackClustersV.size()), nTracksW(trackClustersW.size());
    const unsigned int nShowersU(showerClustersU.size()), nShowersV(showerClustersV.size()), nShowersW(showerClustersW.size());

    if ((1 == nTracksU * nTracksV * nShowersW) && (1 != nTracksU * nShowersV * nTracksW) && (1 != nShowersU * nTracksV * nTracksW))
    {
        pClusterU = *(trackClustersU.begin());
        pClusterV = *(trackClustersV.begin());
        pClusterW = *(showerClustersW.begin());
        return true;
    }

    if ((1 != nTracksU * nTracksV * nShowersW) && (1 == nTracksU * nShowersV * nTracksW) && (1 != nShowersU * nTracksV * nTracksW))
    {
        pClusterU = *(trackClustersU.begin());
        pClusterV = *(showerClustersV.begin());
        pClusterW = *(trackClustersW.begin());
        return true;
    }

    if ((1 != nTracksU * nTracksV * nShowersW) && (1 != nTracksU * nShowersV * nTracksW) && (1 == nShowersU * nTracksV * nTracksW))
    {
        pClusterU = *(showerClustersU.begin());
        pClusterV = *(trackClustersV.begin());
        pClusterW = *(trackClustersW.begin());
        return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ClearTracksTool::TrackShowerShowerAmbiguity(const ClusterList &clusterListU, const ClusterList &clusterListV, const ClusterList &clusterListW,
    Cluster *&pClusterU, Cluster *&pClusterV, Cluster *&pClusterW)
{
    ClusterList trackClustersU, showerClustersU;
    ClearTracksTool::ClassifyClusters(clusterListU, trackClustersU, showerClustersU);

    ClusterList trackClustersV, showerClustersV;
    ClearTracksTool::ClassifyClusters(clusterListV, trackClustersV, showerClustersV);

    ClusterList trackClustersW, showerClustersW;
    ClearTracksTool::ClassifyClusters(clusterListW, trackClustersW, showerClustersW);

    const unsigned int nTracksU(trackClustersU.size()), nTracksV(trackClustersV.size()), nTracksW(trackClustersW.size());
    const unsigned int nShowersU(showerClustersU.size()), nShowersV(showerClustersV.size()), nShowersW(showerClustersW.size());

    if ((1 == nTracksU * nShowersV * nShowersW) && (1 != nShowersU * nTracksV * nShowersW) && (1 != nShowersU * nShowersV * nTracksW))
    {
        pClusterU = *(trackClustersU.begin());
        pClusterV = *(showerClustersV.begin());
        pClusterW = *(showerClustersW.begin());
        return true;
    }

    if ((1 != nTracksU * nShowersV * nShowersW) && (1 == nShowersU * nTracksV * nShowersW) && (1 != nShowersU * nShowersV * nTracksW))
    {
        pClusterU = *(showerClustersU.begin());
        pClusterV = *(trackClustersV.begin());
        pClusterW = *(showerClustersW.begin());
        return true;
    }

    if ((1 != nTracksU * nShowersV * nShowersW) && (1 != nShowersU * nTracksV * nShowersW) && (1 == nShowersU * nShowersV * nTracksW))
    {
        pClusterU = *(showerClustersU.begin());
        pClusterV = *(showerClustersV.begin());
        pClusterW = *(trackClustersW.begin());
        return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ClearTracksTool::ShowerShowerShowerAmbiguity(const ClusterList &clusterListU, const ClusterList &clusterListV, const ClusterList &clusterListW,
    Cluster *&pClusterU, Cluster *&pClusterV, Cluster *&pClusterW)
{
    ClusterList trackClustersU, showerClustersU;
    ClearTracksTool::ClassifyClusters(clusterListU, trackClustersU, showerClustersU);

    if (1 != showerClustersU.size())
        return false;

    ClusterList trackClustersV, showerClustersV;
    ClearTracksTool::ClassifyClusters(clusterListV, trackClustersV, showerClustersV);

    if (1 != showerClustersV.size())
        return false;

    ClusterList trackClustersW, showerClustersW;
    ClearTracksTool::ClassifyClusters(clusterListW, trackClustersW, showerClustersW);

    if (1 != showerClustersW.size())
        return false;

    pClusterU = *(showerClustersU.begin());
    pClusterV = *(showerClustersV.begin());
    pClusterW = *(showerClustersW.begin());

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClearTracksTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_minMatchedFraction = 0.9f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedFraction", m_minMatchedFraction));

    m_minXOverlapFraction = 0.9f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinXOverlapFraction", m_minXOverlapFraction));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
