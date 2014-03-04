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

    ClusterList usedClusters;

    while (true)
    {
        TrackOverlapResult bestTrackOverlapResult;
        Cluster *pBestClusterU(NULL), *pBestClusterV(NULL), *pBestClusterW(NULL);

        for (TensorType::const_iterator iterU = overlapTensor.begin(), iterUEnd = overlapTensor.end(); iterU != iterUEnd; ++iterU)
        {
            Cluster *pClusterU = iterU->first;
            const TensorType::OverlapMatrix &overlapMatrix(iterU->second);

            if (usedClusters.count(pClusterU))
                continue;

            if (overlapMatrix.size() != 1)
                continue;

            for (TensorType::OverlapMatrix::const_iterator iterV = overlapMatrix.begin(), iterVEnd = overlapMatrix.end(); iterV != iterVEnd; ++iterV)
            {
                Cluster *pClusterV = iterV->first;
                const TensorType::OverlapList &overlapList(iterV->second);

                if (usedClusters.count(pClusterV))
                    continue;

                if (overlapList.size() != 1)
                    continue;

                Cluster *pClusterW = overlapList.begin()->first;
                const TrackOverlapResult &trackOverlapResult(overlapList.begin()->second);

                if (usedClusters.count(pClusterW))
                    continue;

                if (trackOverlapResult < bestTrackOverlapResult)
                    continue;

                bestTrackOverlapResult = trackOverlapResult;
                pBestClusterU = pClusterU;
                pBestClusterV = iterV->first;
                pBestClusterW = overlapList.begin()->first;
            }
        }

        if (!pBestClusterU || !pBestClusterV || !pBestClusterW)
            break;

        ProtoParticle protoParticle;
        protoParticle.m_clusterListU.insert(pBestClusterU);
        protoParticle.m_clusterListV.insert(pBestClusterV);
        protoParticle.m_clusterListW.insert(pBestClusterW);
        protoParticleVector.push_back(protoParticle);

        usedClusters.insert(pBestClusterU); usedClusters.insert(pBestClusterV); usedClusters.insert(pBestClusterW);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClearTracksTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    // Read settings from xml file here

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
