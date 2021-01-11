/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/DeltaRayRecoveryTool.cc
 *
 *  @brief  Implementation of the delta ray recovery tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArThreeDReco/LArCosmicRay/DeltaRayRecoveryTool.h"

using namespace pandora;

namespace lar_content
{

DeltaRayRecoveryTool::DeltaRayRecoveryTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayRecoveryTool::Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    bool mergesMade(false);

    this->MakeMerges(pAlgorithm, overlapTensor, mergesMade);

    return mergesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayRecoveryTool::MakeMerges(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor, bool &mergesMade) const
{
    bool mergeMade(true);

    while (mergeMade)
    {
        mergeMade = false;

        ClusterVector sortedKeyClusters;
        overlapTensor.GetSortedKeyClusters(sortedKeyClusters);

        ClusterSet usedKeyClusters;
        for (const Cluster *const pKeyCluster : sortedKeyClusters)
        {
            if (usedKeyClusters.count(pKeyCluster))
                continue;

            ClusterSet checkedClusters;
            TensorType::ElementList elementList;
            pAlgorithm->GetConnectedElements(pKeyCluster, false, elementList, checkedClusters);

            if (elementList.empty())
                continue;

            for (const TensorType::Element &element : elementList)
	        {
                if (usedKeyClusters.count(element.GetClusterU()))
                    continue;

                usedKeyClusters.insert(element.GetClusterU());
            }

            if (this->PickOutGoodMatches(pAlgorithm, elementList))
            {
                mergeMade = true; mergesMade = true;
                break;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayRecoveryTool::PickOutGoodMatches(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::ElementList &elementList) const
{
    ProtoParticleVector protoParticleVector;
    
    bool found(true);

    ClusterSet usedClusters;
    
    while (found)
    {
        found = false;

        float highestHitCount(-std::numeric_limits<float>::max()), bestChiSquared(0.f);
        const Cluster *pBestClusterU(nullptr), *pBestClusterV(nullptr), *pBestClusterW(nullptr);

        for (const TensorType::Element &element : elementList)
        {
            if (element.GetOverlapResult().GetReducedChi2() > 1.f)
                continue;
            
            const Cluster *const pClusterU(element.GetCluster(TPC_VIEW_U)), *const pClusterV(element.GetCluster(TPC_VIEW_V)), *const pClusterW(element.GetCluster(TPC_VIEW_W));
            
            if (usedClusters.count(pClusterU) || usedClusters.count(pClusterV) || usedClusters.count(pClusterW))
                continue;

            const float chiSquared = element.GetOverlapResult().GetReducedChi2();            
            const unsigned int hitSum(pClusterU->GetNCaloHits() + pClusterV->GetNCaloHits() + pClusterW->GetNCaloHits());

            if ((hitSum == highestHitCount) && (chiSquared < bestChiSquared))
            {
                bestChiSquared = chiSquared;
                highestHitCount = hitSum;
                pBestClusterU = pClusterU; pBestClusterV = pClusterV; pBestClusterW = pClusterW;
                
                continue;
            }
            
            if (hitSum > highestHitCount)
            {
                bestChiSquared = chiSquared;
                highestHitCount = hitSum;
                pBestClusterU = pClusterU; pBestClusterV = pClusterV; pBestClusterW = pClusterW;
            }
        }

        if (pBestClusterU && pBestClusterV && pBestClusterW)
        {
            found = true;
            usedClusters.insert(pBestClusterU); usedClusters.insert(pBestClusterV); usedClusters.insert(pBestClusterW);
            
            ProtoParticle protoParticle;
            protoParticle.m_clusterList.push_back(pBestClusterU);
            protoParticle.m_clusterList.push_back(pBestClusterV);
            protoParticle.m_clusterList.push_back(pBestClusterW);
            protoParticleVector.push_back(protoParticle);
        }
    }        

    if (!protoParticleVector.empty())
    {
        pAlgorithm->CreatePfos(protoParticleVector);
        return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DeltaRayRecoveryTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{          
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content
