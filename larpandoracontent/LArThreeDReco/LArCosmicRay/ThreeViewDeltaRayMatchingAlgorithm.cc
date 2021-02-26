/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/ThreeViewDeltaRayMatchingAlgorithm.cc
 *
 *  @brief  Implementation of the three view delta ray matching class
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/ThreeViewDeltaRayMatchingAlgorithm.h"

#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

using namespace pandora;

namespace lar_content
{

ThreeViewDeltaRayMatchingAlgorithm::ThreeViewDeltaRayMatchingAlgorithm() :
    m_minClusterCaloHits(5),
    m_nMaxTensorToolRepeats(10)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ThreeViewDeltaRayMatchingAlgorithm::DoesClusterPassTesorThreshold(const Cluster *const pCluster) const
{
    return (pCluster->GetNCaloHits() >= m_minClusterCaloHits);
}    

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewDeltaRayMatchingAlgorithm::CalculateOverlapResult(const Cluster *const pClusterU, const Cluster *const pClusterV, const Cluster *const pClusterW)
{
    DeltaRayOverlapResult overlapResult;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, this->CalculateOverlapResult(pClusterU, pClusterV, pClusterW, overlapResult));

    if (overlapResult.IsInitialized())
        this->GetMatchingControl().GetOverlapTensor().SetOverlapResult(pClusterU, pClusterV, pClusterW, overlapResult);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeViewDeltaRayMatchingAlgorithm::CalculateOverlapResult(const Cluster *const pClusterU, const Cluster *const pClusterV, const Cluster *const pClusterW,
    DeltaRayOverlapResult &overlapResult) const
{
    float chiSquaredSum(0.f);
    unsigned int nSamplingPoints(0), nMatchedSamplingPoints(0);
    XOverlap xOverlapObject(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f);
    
    StatusCode statusCode(this->PerformThreeViewMatching(pClusterU, pClusterV, pClusterW, chiSquaredSum, nSamplingPoints, nMatchedSamplingPoints, xOverlapObject));

    if (statusCode != STATUS_CODE_SUCCESS)
        return statusCode;

    PfoList commonMuonPfoList;
    this->FindCommonMuonParents(pClusterU, pClusterV, pClusterW, commonMuonPfoList);

    if (commonMuonPfoList.empty())
        return STATUS_CODE_NOT_FOUND;

    overlapResult = DeltaRayOverlapResult(nMatchedSamplingPoints, nSamplingPoints, chiSquaredSum, xOverlapObject, commonMuonPfoList);
    
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewDeltaRayMatchingAlgorithm::FindCommonMuonParents(const Cluster *const pClusterU, const Cluster *const pClusterV, const Cluster *const pClusterW,
    PfoList &commonMuonPfoList) const
{
    ClusterList consideredClustersU, consideredClustersV, consideredClustersW;
    PfoList nearbyMuonPfosU, nearbyMuonPfosV, nearbyMuonPfosW;
    
    this->GetNearbyMuonPfos(pClusterU, consideredClustersU, nearbyMuonPfosU);

    if (nearbyMuonPfosU.empty())
        return;
    
    this->GetNearbyMuonPfos(pClusterV, consideredClustersV, nearbyMuonPfosV);

    if (nearbyMuonPfosV.empty())
        return;
    
    this->GetNearbyMuonPfos(pClusterW, consideredClustersW, nearbyMuonPfosW);

    if (nearbyMuonPfosW.empty())
        return;    

    for (const ParticleFlowObject *const pNearbyMuonU : nearbyMuonPfosU)
    {
        for (const ParticleFlowObject *const pNearbyMuonV : nearbyMuonPfosV)
        {
            for (const ParticleFlowObject *const pNearbyMuonW : nearbyMuonPfosW)
            {
                if ((pNearbyMuonU == pNearbyMuonV) && (pNearbyMuonV == pNearbyMuonW))
                    commonMuonPfoList.push_back(pNearbyMuonU);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewDeltaRayMatchingAlgorithm::ExamineOverlapContainer()
{
    unsigned int repeatCounter(0);

    for (TensorToolVector::const_iterator toolIter = m_algorithmToolVector.begin(); toolIter != m_algorithmToolVector.end(); )
    {
        if ((*toolIter)->Run(this, this->GetMatchingControl().GetOverlapTensor()))
        {
            toolIter = m_algorithmToolVector.begin();

            if (++repeatCounter > m_nMaxTensorToolRepeats)
                break;
        }
        else
        {
            ++toolIter;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeViewDeltaRayMatchingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{   
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MuonPfoListName", m_muonPfoListName));
    
    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle,
        "DeltaRayTools", algorithmToolVector));

    for (AlgorithmToolVector::const_iterator iter = algorithmToolVector.begin(), iterEnd = algorithmToolVector.end(); iter != iterEnd; ++iter)
    {
        DeltaRayTensorTool *const pDeltaRayTensorTool(dynamic_cast<DeltaRayTensorTool*>(*iter));

        if (!pDeltaRayTensorTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolVector.push_back(pDeltaRayTensorTool);
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, xmlHandle,
        "ClusterRebuilding", m_reclusteringAlgorithmName));    

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterCaloHits", m_minClusterCaloHits));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NMaxTensorToolRepeats", m_nMaxTensorToolRepeats));    
    
    return BaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content

