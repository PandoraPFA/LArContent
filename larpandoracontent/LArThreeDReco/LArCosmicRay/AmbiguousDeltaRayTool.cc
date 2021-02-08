/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/AmbiguousDeltaRayTool.cc
 *
 *  @brief  Implementation of the delta ray merge tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/AmbiguousDeltaRayTool.h"

using namespace pandora;

namespace lar_content
{

AmbiguousDeltaRayTool::AmbiguousDeltaRayTool() :
    m_maxGoodMatchReducedChiSquared(1.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool AmbiguousDeltaRayTool::Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    this->ExamineConnectedElements(pAlgorithm, overlapTensor);

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AmbiguousDeltaRayTool::ExamineConnectedElements(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor) const
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

            TensorType::ElementList elementList;
            overlapTensor.GetConnectedElements(pKeyCluster, true, elementList);            

            for (const TensorType::Element &element : elementList)
            {
                if (usedKeyClusters.count(element.GetClusterU()))
                    continue;

                usedKeyClusters.insert(element.GetClusterU());
            }

            if (elementList.size() < 2)
                continue;

            if (this->PickOutGoodMatches(pAlgorithm, elementList))
            {
                mergeMade = true;
                break;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool AmbiguousDeltaRayTool::PickOutGoodMatches(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, const TensorType::ElementList &elementList) const
{
    ProtoParticleVector protoParticleVector;
    
    bool found(true);

    ClusterSet usedClusters;
    
    while (found)
    {
        found = false;

        unsigned int highestHitCount(0);
        float bestChiSquared(std::numeric_limits<float>::max());
        const Cluster *pBestClusterU(nullptr), *pBestClusterV(nullptr), *pBestClusterW(nullptr);

        for (const TensorType::Element &element : elementList)
        {            
            const Cluster *const pClusterU(element.GetCluster(TPC_VIEW_U)), *const pClusterV(element.GetCluster(TPC_VIEW_V)), *const pClusterW(element.GetCluster(TPC_VIEW_W));
            
            if (usedClusters.count(pClusterU) || usedClusters.count(pClusterV) || usedClusters.count(pClusterW))
                continue;

            const float chiSquared = element.GetOverlapResult().GetReducedChi2();

            if (chiSquared > m_maxGoodMatchReducedChiSquared)
                continue;
            
            const unsigned int hitSum(pClusterU->GetNCaloHits() + pClusterV->GetNCaloHits() + pClusterW->GetNCaloHits());

            if ((hitSum > highestHitCount) || ((hitSum == highestHitCount) && (chiSquared < bestChiSquared)))
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

StatusCode AmbiguousDeltaRayTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxGoodMatchReducedChiSquared", m_maxGoodMatchReducedChiSquared)); 
    
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content
