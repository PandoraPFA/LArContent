/**
 *  @file   LArContent/src/LArThreeDReco/LArShowerFragments/ThreeDRemnantsAlgorithm.cc
 *
 *  @brief  Implementation of the three dimensional remnants algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArThreeDReco/LArShowerFragments/ThreeDRemnantsAlgorithm.h"

#include "LArHelpers/LArGeometryHelper.h"
#include "LArHelpers/LArClusterHelper.h"

using namespace pandora;

namespace lar_content
{

void ThreeDRemnantsAlgorithm::SelectInputClusters(const ClusterList *const pInputClusterList, ClusterList &selectedClusterList) const
{
    for (ClusterList::const_iterator iter = pInputClusterList->begin(), iterEnd = pInputClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (!pCluster->IsAvailable())
            continue;

        if (pCluster->GetNCaloHits() < m_minClusterCaloHits)
            continue;

        selectedClusterList.insert(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDRemnantsAlgorithm::SetPfoParameters(const ProtoParticle &protoParticle, PandoraContentApi::ParticleFlowObject::Parameters &pfoParameters) const
{
    // TODO Correct these placeholder parameters
    pfoParameters.m_particleId = E_MINUS; // Shower
    pfoParameters.m_charge = PdgTable::GetParticleCharge(pfoParameters.m_particleId.Get());
    pfoParameters.m_mass = PdgTable::GetParticleMass(pfoParameters.m_particleId.Get());
    pfoParameters.m_energy = 0.f;
    pfoParameters.m_momentum = CartesianVector(0.f, 0.f, 0.f);
    pfoParameters.m_clusterList.insert(protoParticle.m_clusterListU.begin(), protoParticle.m_clusterListU.end());
    pfoParameters.m_clusterList.insert(protoParticle.m_clusterListV.begin(), protoParticle.m_clusterListV.end());
    pfoParameters.m_clusterList.insert(protoParticle.m_clusterListW.begin(), protoParticle.m_clusterListW.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDRemnantsAlgorithm::CalculateOverlapResult(Cluster *pClusterU, Cluster *pClusterV, Cluster *pClusterW)
{
    try
    {
        // Requirements on X matching
        float xMinU(0.f), xMinV(0.f), xMinW(0.f), xMaxU(0.f), xMaxV(0.f), xMaxW(0.f);
        LArClusterHelper::GetClusterSpanX(pClusterU, xMinU, xMaxU);
        LArClusterHelper::GetClusterSpanX(pClusterV, xMinV, xMaxV);
        LArClusterHelper::GetClusterSpanX(pClusterW, xMinW, xMaxW);

        const float xMin(std::max(xMinU, std::max(xMinV, xMinW)) - m_xOverlapWindow);
        const float xMax(std::min(xMaxU, std::min(xMaxV, xMaxW)) + m_xOverlapWindow);
        const float xOverlap(xMax - xMin);

        if (xOverlap < std::numeric_limits<float>::epsilon())
            return;

        // Requirements on 3D matching
        const float u(LArClusterHelper::GetAverageZ(pClusterU, xMin, xMax));
        const float v(LArClusterHelper::GetAverageZ(pClusterV, xMin, xMax));
        const float w(LArClusterHelper::GetAverageZ(pClusterW, xMin, xMax));

        const float uv2w(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, u, v));
        const float vw2u(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_V, TPC_VIEW_W, v, w));
        const float wu2v(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_W, TPC_VIEW_U, w, u));

        const float pseudoChi2(((u - vw2u) * (u - vw2u) + (v - wu2v) * (v - wu2v) + (w - uv2w) * (w - uv2w)) / 3.f);

        if (pseudoChi2 > m_pseudoChi2Cut)
            return;

        m_overlapTensor.SetOverlapResult(pClusterU, pClusterV, pClusterW, true);
    }
    catch(StatusCodeException &statusCodeException)
    {
        if (STATUS_CODE_NOT_FOUND != statusCodeException.GetStatusCode())
            throw statusCodeException;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDRemnantsAlgorithm::ExamineTensor()
{
    unsigned int repeatCounter(0);

    for (RemnantTensorToolList::const_iterator iter = m_algorithmToolList.begin(), iterEnd = m_algorithmToolList.end(); iter != iterEnd; )
    {
        if ((*iter)->Run(this, m_overlapTensor))
        {
            iter = m_algorithmToolList.begin();

            if (++repeatCounter > m_nMaxTensorToolRepeats)
                break;
        }
        else
        {
            ++iter;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDRemnantsAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    AlgorithmToolList algorithmToolList;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle,
        "TrackTools", algorithmToolList));

    for (AlgorithmToolList::const_iterator iter = algorithmToolList.begin(), iterEnd = algorithmToolList.end(); iter != iterEnd; ++iter)
    {
        RemnantTensorTool *pRemnantTensorTool(dynamic_cast<RemnantTensorTool*>(*iter));

        if (NULL == pRemnantTensorTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolList.push_back(pRemnantTensorTool);
    }

    m_nMaxTensorToolRepeats = 5000;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NMaxTensorToolRepeats", m_nMaxTensorToolRepeats));

    m_minClusterCaloHits = 5;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterCaloHits", m_minClusterCaloHits));

    m_xOverlapWindow = 2.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OverlapWindow", m_xOverlapWindow));

    m_pseudoChi2Cut = 20.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PseudoChi2Cut", m_pseudoChi2Cut));

    return ThreeDBaseAlgorithm<float>::ReadSettings(xmlHandle);
}

} // namespace lar_content
