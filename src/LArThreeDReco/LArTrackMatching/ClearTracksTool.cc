/**
 *  @file   LArContent/src/LArThreeDReco/LArTrackMatching/ClearTracksTool.cc
 * 
 *  @brief  Implementation of the clear tracks tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "LArThreeDReco/LArTrackMatching/ClearTracksTool.h"

using namespace pandora;

namespace lar
{

StatusCode ClearTracksTool::Run(ThreeDTransverseTracksAlgorithm *pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraSettings::ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << m_algorithmToolType << std::endl;

    TensorType::ElementList elementList;
    overlapTensor.GetUnambiguousElements(true, elementList);
    this->CreateThreeDParticles(pAlgorithm, elementList);

    elementList.clear();
    overlapTensor.GetUnambiguousElements(true, &TrackTrackTrackAmbiguity, elementList);
    this->CreateThreeDParticles(pAlgorithm, elementList);

    elementList.clear();
    overlapTensor.GetUnambiguousElements(true, &TrackTrackShowerAmbiguity, elementList);
    this->CreateThreeDParticles(pAlgorithm, elementList);

    elementList.clear();
    overlapTensor.GetUnambiguousElements(true, &TrackShowerShowerAmbiguity, elementList);
    this->CreateThreeDParticles(pAlgorithm, elementList);

    elementList.clear();
    overlapTensor.GetUnambiguousElements(true, &ShowerShowerShowerAmbiguity, elementList);
    this->CreateThreeDParticles(pAlgorithm, elementList);

    elementList.clear();
    this->ResolveSimpleAmbiguities(overlapTensor, elementList);
    this->CreateThreeDParticles(pAlgorithm, elementList);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ClearTracksTool::CreateThreeDParticles(ThreeDTransverseTracksAlgorithm *pAlgorithm, const TensorType::ElementList &elementList) const
{
    ProtoParticleVector protoParticleVector;

    for (TensorType::ElementList::const_iterator iter = elementList.begin(), iterEnd = elementList.end(); iter != iterEnd; ++iter)
    {
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

    pAlgorithm->CreateThreeDParticles(protoParticleVector);
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

void ClearTracksTool::ResolveSimpleAmbiguities(const TensorType &overlapTensor, TensorType::ElementList &elementList) const
{
    ClusterList usedClusters;

    for (TensorType::const_iterator iterU = overlapTensor.begin(), iterUEnd = overlapTensor.end(); iterU != iterUEnd; ++iterU)
    {
        unsigned int nU(0), nV(0), nW(0);
        TensorType::ElementList localElementList;
        overlapTensor.GetConnectedElements(iterU->first, true, localElementList, nU, nV, nW);

        if (nU * nV * nW < 2)
            continue;

        unsigned int nElementsPassingTrackOverlapCut(0), highestMSP(0), secondHighestMSP(0);
        TensorType::ElementList::const_iterator passedTrackOverlapElement(localElementList.end());
        TensorType::ElementList::const_iterator highestMSPElement(localElementList.end());

        for (TensorType::ElementList::const_iterator eIter = localElementList.begin(); eIter != localElementList.end(); ++eIter)
        {
            const TransverseOverlapResult::XOverlap &xOverlap(eIter->GetOverlapResult().GetXOverlap());

            if ((xOverlap.GetXSpanU() > std::numeric_limits<float>::epsilon()) && (xOverlap.GetXOverlapSpan() / xOverlap.GetXSpanU() > m_minXOverlapFraction) &&
                (xOverlap.GetXSpanV() > std::numeric_limits<float>::epsilon()) && (xOverlap.GetXOverlapSpan() / xOverlap.GetXSpanV() > m_minXOverlapFraction) &&
                (xOverlap.GetXSpanW() > std::numeric_limits<float>::epsilon()) && (xOverlap.GetXOverlapSpan() / xOverlap.GetXSpanW() > m_minXOverlapFraction))
            {
                ++nElementsPassingTrackOverlapCut;
                passedTrackOverlapElement = eIter;
            }

            const unsigned int nMatchedSamplingPoints(eIter->GetOverlapResult().GetNMatchedSamplingPoints());

            if (nMatchedSamplingPoints > highestMSP)
            {
                highestMSP = nMatchedSamplingPoints;
                highestMSPElement = eIter;
            }
            else if (nMatchedSamplingPoints > secondHighestMSP)
            {
                secondHighestMSP = nMatchedSamplingPoints;
            }
        }

        if ((1 != nElementsPassingTrackOverlapCut) || (highestMSPElement != passedTrackOverlapElement) || (highestMSP < m_minMatchedSamplingPointRatio * secondHighestMSP))
            continue;

        if (usedClusters.count(passedTrackOverlapElement->GetClusterU()) || usedClusters.count(passedTrackOverlapElement->GetClusterV()) || usedClusters.count(passedTrackOverlapElement->GetClusterW()))
            continue;

        elementList.push_back(*passedTrackOverlapElement);
        usedClusters.insert(passedTrackOverlapElement->GetClusterU());
        usedClusters.insert(passedTrackOverlapElement->GetClusterV());
        usedClusters.insert(passedTrackOverlapElement->GetClusterW());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClearTracksTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_minXOverlapFraction = 0.9f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinXOverlapFraction", m_minXOverlapFraction));

    m_minMatchedSamplingPointRatio = 2;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedSamplingPointRatio", m_minMatchedSamplingPointRatio));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
