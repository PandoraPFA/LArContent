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

bool ClearTracksTool::Run(const SlidingFitResultMap &slidingFitResultMap, TrackOverlapTensor &overlapTensor, ProtoParticleVector &protoParticleVector)
{
    std::cout << "ClearTracksTool::Run() " << std::endl;
    TrackOverlapResult bestTrackOverlapResult;
    Cluster *pBestClusterU(NULL), *pBestClusterV(NULL), *pBestClusterW(NULL);

    const ClusterList &clusterListU(overlapTensor.GetClusterListU());
    const ClusterList &clusterListV(overlapTensor.GetClusterListV());
    const ClusterList &clusterListW(overlapTensor.GetClusterListW());

    for (ClusterList::const_iterator iterU = clusterListU.begin(), iterUEnd = clusterListU.end(); iterU != iterUEnd; ++iterU)
    {
        for (ClusterList::const_iterator iterV = clusterListV.begin(), iterVEnd = clusterListV.end(); iterV != iterVEnd; ++iterV)
        {
            for (ClusterList::const_iterator iterW = clusterListW.begin(), iterWEnd = clusterListW.end(); iterW != iterWEnd; ++iterW)
            {
                try
                {
                    const TrackOverlapResult &trackOverlapResult(overlapTensor.GetOverlapResult(*iterU, *iterV, *iterW));

                    if (trackOverlapResult < bestTrackOverlapResult)
                        continue;

                    bestTrackOverlapResult = trackOverlapResult;
                    pBestClusterU = *iterU;
                    pBestClusterV = *iterV;
                    pBestClusterW = *iterW;
                }
                catch (StatusCodeException &)
                {
                }
            }
        }
    }

    if (!pBestClusterU || !pBestClusterV || !pBestClusterW)
        return false;

    ProtoParticle protoParticle;
    protoParticle.m_clusterListU.insert(pBestClusterU);
    protoParticle.m_clusterListV.insert(pBestClusterV);
    protoParticle.m_clusterListW.insert(pBestClusterW);

    protoParticleVector.push_back(protoParticle);
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClearTracksTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    // Read settings from xml file here

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
