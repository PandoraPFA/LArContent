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

StatusCode ClearTracksTool::Run(const SlidingFitResultMap &slidingFitResultMap, TensorType &overlapTensor, ProtoParticleVector &protoParticleVector)
{
    std::cout << "ClearTracksTool::Run() " << std::endl;

    // TODO x-overlap and definitions of unambiguous using x-overlap
    TensorType::ElementList elementList;
    overlapTensor.GetUnambiguousElements(true, elementList);
    overlapTensor.GetUnambiguousElements(true, &TrackTrackTrackAmbiguity, elementList);

    for (TensorType::ElementList::const_iterator iter = elementList.begin(), iterEnd = elementList.end(); iter != iterEnd; ++iter)
    {
        ProtoParticle protoParticle;
        protoParticle.m_clusterListU.insert(iter->GetClusterU());
        protoParticle.m_clusterListV.insert(iter->GetClusterV());
        protoParticle.m_clusterListW.insert(iter->GetClusterW());
        protoParticleVector.push_back(protoParticle);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ClearTracksTool::TrackTrackTrackAmbiguity(const ClusterList &clusterListU, const ClusterList &clusterListV, const ClusterList &clusterListW,
    Cluster *&pClusterU, Cluster *&pClusterV, Cluster *&pClusterW)
{
    if ((1 == clusterListU.size()) || (1 == clusterListV.size()) || (1 == clusterListW.size()))
        return false;

    unsigned int nTracksU(0), nTracksV(0), nTracksW(0);
    Cluster *pTrackClusterU(NULL), *pTrackClusterV(NULL), *pTrackClusterW(NULL);

    // ATTN For later development
    unsigned int nShowersU(0), nShowersV(0), nShowersW(0);
    Cluster *pShowerClusterU(NULL), *pShowerClusterV(NULL), *pShowerClusterW(NULL);

    for (ClusterList::const_iterator iter = clusterListU.begin(), iterEnd = clusterListU.end(); iter != iterEnd; ++iter)
    {
        if ((*iter)->IsMipTrack())
        {
            ++nTracksU;
            pTrackClusterU = *iter;
        }
        else
        {
            ++nShowersU;
            pShowerClusterU = *iter;
        }
    }

    for (ClusterList::const_iterator iter = clusterListV.begin(), iterEnd = clusterListV.end(); iter != iterEnd; ++iter)
    {
        if ((*iter)->IsMipTrack())
        {
            ++nTracksV;
            pTrackClusterV = *iter;
        }
        else
        {
            ++nShowersU;
            pShowerClusterU = *iter;
        }
    }

    for (ClusterList::const_iterator iter = clusterListW.begin(), iterEnd = clusterListW.end(); iter != iterEnd; ++iter)
    {
        if ((*iter)->IsMipTrack())
        {
            ++nTracksW;
            pTrackClusterW = *iter;
        }
        else
        {
            ++nShowersU;
            pShowerClusterU = *iter;
        }
    }

    if ((1 != nTracksU) || (1 != nTracksV) || (1 != nTracksW))
        return false;

    if ((NULL == pTrackClusterU) || (NULL == pTrackClusterV) || (NULL == pTrackClusterW))
        return false;

    pClusterU = pTrackClusterU;
    pClusterV = pTrackClusterV;
    pClusterW = pTrackClusterW;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClearTracksTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    // Read settings from xml file here

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
