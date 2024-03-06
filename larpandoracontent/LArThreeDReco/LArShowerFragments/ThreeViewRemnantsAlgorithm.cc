/**
 *  @file   larpandoracontent/LArThreeDReco/LArShowerFragments/ThreeViewRemnantsAlgorithm.cc
 *
 *  @brief  Implementation of the three view remnants algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArThreeDReco/LArShowerFragments/ThreeViewRemnantsAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

using namespace pandora;

namespace lar_content
{

ThreeViewRemnantsAlgorithm::ThreeViewRemnantsAlgorithm() :
    m_nMaxTensorToolRepeats(1000), m_minClusterCaloHits(5), m_xOverlapWindow(2.f), m_pseudoChi2Cut(10.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewRemnantsAlgorithm::SelectInputClusters(const ClusterList *const pInputClusterList, ClusterList &selectedClusterList) const
{
    for (ClusterList::const_iterator iter = pInputClusterList->begin(), iterEnd = pInputClusterList->end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pCluster = *iter;

        if (!pCluster->IsAvailable())
            continue;

        if (pCluster->GetNCaloHits() < m_minClusterCaloHits)
            continue;

        selectedClusterList.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewRemnantsAlgorithm::CalculateOverlapResult(const Cluster *const pClusterU, const Cluster *const pClusterV, const Cluster *const pClusterW)
{
    // Requirements on X matching
    float xMinU(0.f), xMinV(0.f), xMinW(0.f), xMaxU(0.f), xMaxV(0.f), xMaxW(0.f);
    pClusterU->GetClusterSpanX(xMinU, xMaxU);
    pClusterV->GetClusterSpanX(xMinV, xMaxV);
    pClusterW->GetClusterSpanX(xMinW, xMaxW);

    const float xMin(std::max(xMinU, std::max(xMinV, xMinW)) - m_xOverlapWindow);
    const float xMax(std::min(xMaxU, std::min(xMaxV, xMaxW)) + m_xOverlapWindow);
    const float xOverlap(xMax - xMin);

    if (xOverlap < std::numeric_limits<float>::epsilon())
        return;

    // Requirements on 3D matching
    float u(std::numeric_limits<float>::max()), v(std::numeric_limits<float>::max()), w(std::numeric_limits<float>::max());

    if ((STATUS_CODE_SUCCESS != LArClusterHelper::GetAverageZ(pClusterU, xMin, xMax, u)) ||
        (STATUS_CODE_SUCCESS != LArClusterHelper::GetAverageZ(pClusterV, xMin, xMax, v)) ||
        (STATUS_CODE_SUCCESS != LArClusterHelper::GetAverageZ(pClusterW, xMin, xMax, w)))
    {
        return;
    }

    const float uv2w(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, u, v));
    const float vw2u(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_V, TPC_VIEW_W, v, w));
    const float wu2v(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_W, TPC_VIEW_U, w, u));

    const float pseudoChi2(((u - vw2u) * (u - vw2u) + (v - wu2v) * (v - wu2v) + (w - uv2w) * (w - uv2w)) / 3.f);

    if (pseudoChi2 > m_pseudoChi2Cut)
        return;

    // ATTN Essentially a boolean result; actual value matters only so as to ensure that overlap results can be sorted
    const float hackValue(
        pseudoChi2 + pClusterU->GetElectromagneticEnergy() + pClusterV->GetElectromagneticEnergy() + pClusterW->GetElectromagneticEnergy());
    this->GetMatchingControl().GetOverlapTensor().SetOverlapResult(pClusterU, pClusterV, pClusterW, hackValue);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewRemnantsAlgorithm::ExamineOverlapContainer()
{
    unsigned int repeatCounter(0);

    for (RemnantTensorToolVector::const_iterator iter = m_algorithmToolVector.begin(), iterEnd = m_algorithmToolVector.end(); iter != iterEnd;)
    {
        if ((*iter)->Run(this, this->GetMatchingControl().GetOverlapTensor()))
        {
            iter = m_algorithmToolVector.begin();

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

StatusCode ThreeViewRemnantsAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "TrackTools", algorithmToolVector));

    for (AlgorithmToolVector::const_iterator iter = algorithmToolVector.begin(), iterEnd = algorithmToolVector.end(); iter != iterEnd; ++iter)
    {
        RemnantTensorTool *const pRemnantTensorTool(dynamic_cast<RemnantTensorTool *>(*iter));

        if (!pRemnantTensorTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolVector.push_back(pRemnantTensorTool);
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NMaxTensorToolRepeats", m_nMaxTensorToolRepeats));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterCaloHits", m_minClusterCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "OverlapWindow", m_xOverlapWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PseudoChi2Cut", m_pseudoChi2Cut));

    return BaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
