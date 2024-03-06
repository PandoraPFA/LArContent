/**
 *  @file   larpandoracontent/LArThreeDReco/LArPfoRecovery/ParticleRecoveryAlgorithm.cc
 *
 *  @brief  Implementation of the particle recovery algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"

#include "larpandoracontent/LArObjects/LArPointingCluster.h"

#include "larpandoracontent/LArThreeDReco/LArPfoRecovery/ParticleRecoveryAlgorithm.h"

#include <algorithm>

using namespace pandora;

namespace lar_content
{

ParticleRecoveryAlgorithm::ParticleRecoveryAlgorithm() :
    m_checkGaps(true),
    m_minClusterCaloHits(20),
    m_minClusterLengthSquared(5.f * 5.f),
    m_minClusterXSpan(0.25f),
    m_vertexClusterMode(false),
    m_minVertexLongitudinalDistance(-2.5f),
    m_maxVertexTransverseDistance(1.5f),
    m_minXOverlapFraction(0.75f),
    m_minXOverlapFractionGaps(0.75f),
    m_sampleStepSize(0.5f),
    m_slidingFitHalfWindow(10),
    m_pseudoChi2Cut(5.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ParticleRecoveryAlgorithm::Run()
{
    ClusterList inputClusterListU, inputClusterListV, inputClusterListW;
    this->GetInputClusters(inputClusterListU, inputClusterListV, inputClusterListW);

    ClusterList selectedClusterListU, selectedClusterListV, selectedClusterListW;
    this->SelectInputClusters(inputClusterListU, selectedClusterListU);
    this->SelectInputClusters(inputClusterListV, selectedClusterListV);
    this->SelectInputClusters(inputClusterListW, selectedClusterListW);

    SimpleOverlapTensor overlapTensor;
    this->FindOverlaps(selectedClusterListU, selectedClusterListV, overlapTensor);
    this->FindOverlaps(selectedClusterListV, selectedClusterListW, overlapTensor);
    this->FindOverlaps(selectedClusterListW, selectedClusterListU, overlapTensor);
    this->ExamineTensor(overlapTensor);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParticleRecoveryAlgorithm::GetInputClusters(ClusterList &inputClusterListU, ClusterList &inputClusterListV, ClusterList &inputClusterListW) const
{
    for (StringVector::const_iterator iter = m_inputClusterListNames.begin(), iterEnd = m_inputClusterListNames.end(); iter != iterEnd; ++iter)
    {
        const ClusterList *pClusterList(NULL);
        PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, *iter, pClusterList));

        if (!pClusterList || pClusterList->empty())
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "ParticleRecoveryAlgorithm: unable to find cluster list " << *iter << std::endl;

            continue;
        }

        for (ClusterList::const_iterator cIter = pClusterList->begin(), cIterEnd = pClusterList->end(); cIter != cIterEnd; ++cIter)
        {
            const Cluster *const pCluster(*cIter);
            const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));

            if ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType))
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

            ClusterList &clusterList((TPC_VIEW_U == hitType)   ? inputClusterListU
                                     : (TPC_VIEW_V == hitType) ? inputClusterListV
                                                               : inputClusterListW);
            clusterList.push_back(pCluster);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParticleRecoveryAlgorithm::SelectInputClusters(const ClusterList &inputClusterList, ClusterList &selectedClusterList) const
{
    if (m_vertexClusterMode)
    {
        ClusterList vertexClusterList;
        this->VertexClusterSelection(inputClusterList, vertexClusterList);
        this->StandardClusterSelection(vertexClusterList, selectedClusterList);
    }
    else
    {
        this->StandardClusterSelection(inputClusterList, selectedClusterList);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParticleRecoveryAlgorithm::StandardClusterSelection(const ClusterList &inputClusterList, ClusterList &selectedClusterList) const
{
    for (ClusterList::const_iterator iter = inputClusterList.begin(), iterEnd = inputClusterList.end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pCluster = *iter;

        if (!pCluster->IsAvailable())
            continue;

        if (pCluster->GetNCaloHits() < m_minClusterCaloHits)
            continue;

        if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLengthSquared)
            continue;

        float xMin(0.f), xMax(0.f);
        pCluster->GetClusterSpanX(xMin, xMax);

        if ((xMax - xMin) < m_minClusterXSpan)
            continue;

        selectedClusterList.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParticleRecoveryAlgorithm::VertexClusterSelection(const ClusterList &inputClusterList, ClusterList &selectedClusterList) const
{
    CartesianPointVector vertexList;

    for (ClusterList::const_iterator iter = inputClusterList.begin(), iterEnd = inputClusterList.end(); iter != iterEnd; ++iter)
    {
        try
        {
            if (!(*iter)->IsAvailable())
            {
                const LArPointingCluster pointingCluster(*iter);
                vertexList.push_back(pointingCluster.GetInnerVertex().GetPosition());
                vertexList.push_back(pointingCluster.GetOuterVertex().GetPosition());
            }
        }
        catch (StatusCodeException &)
        {
        }
    }

    for (ClusterList::const_iterator iter = inputClusterList.begin(), iterEnd = inputClusterList.end(); iter != iterEnd; ++iter)
    {
        try
        {
            const Cluster *const pCluster = *iter;

            if (!pCluster->IsAvailable())
                continue;

            const LArPointingCluster pointingCluster(pCluster);

            for (CartesianPointVector::const_iterator vIter = vertexList.begin(), vIterEnd = vertexList.end(); vIter != vIterEnd; ++vIter)
            {
                if (LArPointingClusterHelper::IsNode(*vIter, pointingCluster.GetInnerVertex(), m_minVertexLongitudinalDistance, m_maxVertexTransverseDistance) ||
                    LArPointingClusterHelper::IsNode(*vIter, pointingCluster.GetOuterVertex(), m_minVertexLongitudinalDistance, m_maxVertexTransverseDistance))
                {
                    selectedClusterList.push_back(pCluster);
                    break;
                }
            }
        }
        catch (StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParticleRecoveryAlgorithm::FindOverlaps(const ClusterList &clusterList1, const ClusterList &clusterList2, SimpleOverlapTensor &overlapTensor) const
{
    for (ClusterList::const_iterator iter1 = clusterList1.begin(), iter1End = clusterList1.end(); iter1 != iter1End; ++iter1)
    {
        for (ClusterList::const_iterator iter2 = clusterList2.begin(), iter2End = clusterList2.end(); iter2 != iter2End; ++iter2)
        {
            if (this->IsOverlap(*iter1, *iter2))
                overlapTensor.AddAssociation(*iter1, *iter2);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ParticleRecoveryAlgorithm::IsOverlap(const Cluster *const pCluster1, const Cluster *const pCluster2) const
{
    if (LArClusterHelper::GetClusterHitType(pCluster1) == LArClusterHelper::GetClusterHitType(pCluster2))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    if ((0 == pCluster1->GetNCaloHits()) || (0 == pCluster2->GetNCaloHits()))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    float xMin1(0.f), xMax1(0.f), xMin2(0.f), xMax2(0.f);
    pCluster1->GetClusterSpanX(xMin1, xMax1);
    pCluster2->GetClusterSpanX(xMin2, xMax2);

    const float xSpan1(xMax1 - xMin1), xSpan2(xMax2 - xMin2);

    if ((xSpan1 < std::numeric_limits<float>::epsilon()) || (xSpan2 < std::numeric_limits<float>::epsilon()))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const float xOverlap(std::min(xMax1, xMax2) - std::max(xMin1, xMin2));

    float xOverlapFraction1(xOverlap / xSpan1), xOverlapFraction2(xOverlap / xSpan2);

    if (m_checkGaps)
        this->CalculateEffectiveOverlapFractions(pCluster1, xMin1, xMax1, pCluster2, xMin2, xMax2, xOverlapFraction1, xOverlapFraction2);

    if ((xOverlapFraction1 < m_minXOverlapFraction) || (xOverlapFraction2 < m_minXOverlapFraction))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParticleRecoveryAlgorithm::CalculateEffectiveOverlapFractions(const Cluster *const pCluster1, const float xMin1, const float xMax1,
    const Cluster *const pCluster2, const float xMin2, const float xMax2, float &xOverlapFraction1, float &xOverlapFraction2) const
{
    if (PandoraContentApi::GetGeometry(*this)->GetDetectorGapList().empty())
        return;

    const float xMin(std::min(xMin1, xMin2));
    const float xMax(std::max(xMax1, xMax2));
    float xMinEff1(xMin1), xMaxEff1(xMax1), xMinEff2(xMin2), xMaxEff2(xMax2);

    this->CalculateEffectiveSpan(pCluster1, xMin, xMax, xMinEff1, xMaxEff1);
    this->CalculateEffectiveSpan(pCluster2, xMin, xMax, xMinEff2, xMaxEff2);

    const float effectiveXSpan1(xMaxEff1 - xMinEff1), effectiveXSpan2(xMaxEff2 - xMinEff2);
    const float effectiveXOverlapSpan(std::min(xMaxEff1, xMaxEff2) - std::max(xMinEff1, xMinEff2));

    if ((effectiveXSpan1 > std::numeric_limits<float>::epsilon()) && (effectiveXSpan2 > std::numeric_limits<float>::epsilon()))
    {
        xOverlapFraction1 = std::min(1.f, (effectiveXOverlapSpan / effectiveXSpan1));
        xOverlapFraction2 = std::min(1.f, (effectiveXOverlapSpan / effectiveXSpan2));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParticleRecoveryAlgorithm::CalculateEffectiveSpan(
    const pandora::Cluster *const pCluster, const float xMin, const float xMax, float &xMinEff, float &xMaxEff) const
{
    // TODO cache sliding linear fit results and optimise protection against exceptions from TwoDSlidingFitResult and IsXSamplingPointInGap
    try
    {
        const float slidingFitPitch(LArGeometryHelper::GetWirePitch(this->GetPandora(), LArClusterHelper::GetClusterHitType(pCluster)));

        const TwoDSlidingFitResult slidingFitResult(pCluster, m_slidingFitHalfWindow, slidingFitPitch);

        const int nSamplingPointsLeft(1 + static_cast<int>((xMinEff - xMin) / m_sampleStepSize));
        const int nSamplingPointsRight(1 + static_cast<int>((xMax - xMaxEff) / m_sampleStepSize));
        float dxMin(0.f), dxMax(0.f);

        for (int iSample = 1; iSample <= nSamplingPointsLeft; ++iSample)
        {
            const float xSample(std::max(xMin, xMinEff - static_cast<float>(iSample) * m_sampleStepSize));

            if (!LArGeometryHelper::IsXSamplingPointInGap(this->GetPandora(), xSample, slidingFitResult, m_sampleStepSize))
                break;

            dxMin = xMinEff - xSample;
        }

        for (int iSample = 1; iSample <= nSamplingPointsRight; ++iSample)
        {
            const float xSample(std::min(xMax, xMaxEff + static_cast<float>(iSample) * m_sampleStepSize));

            if (!LArGeometryHelper::IsXSamplingPointInGap(this->GetPandora(), xSample, slidingFitResult, m_sampleStepSize))
                break;

            dxMax = xSample - xMaxEff;
        }

        xMinEff = xMinEff - dxMin;
        xMaxEff = xMaxEff + dxMax;
    }
    catch (StatusCodeException &)
    {
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParticleRecoveryAlgorithm::ExamineTensor(const SimpleOverlapTensor &overlapTensor) const
{
    ClusterVector sortedKeyClusters(overlapTensor.GetKeyClusters().begin(), overlapTensor.GetKeyClusters().end());
    std::sort(sortedKeyClusters.begin(), sortedKeyClusters.end(), LArClusterHelper::SortByNHits);

    for (const Cluster *const pKeyCluster : sortedKeyClusters)
    {
        ClusterList clusterListU, clusterListV, clusterListW;

        overlapTensor.GetConnectedElements(pKeyCluster, true, clusterListU, clusterListV, clusterListW);
        const unsigned int nU(clusterListU.size()), nV(clusterListV.size()), nW(clusterListW.size());

        if ((0 == nU * nV) && (0 == nV * nW) && (0 == nW * nU))
            continue;

        ClusterList clusterListAll;
        clusterListAll.insert(clusterListAll.end(), clusterListU.begin(), clusterListU.end());
        clusterListAll.insert(clusterListAll.end(), clusterListV.begin(), clusterListV.end());
        clusterListAll.insert(clusterListAll.end(), clusterListW.begin(), clusterListW.end());

        if ((1 == nU * nV * nW) && this->CheckConsistency(*(clusterListU.begin()), *(clusterListV.begin()), *(clusterListW.begin())))
        {
            this->CreateTrackParticle(clusterListAll);
        }
        else if ((0 == nU * nV * nW) && ((1 == nU && 1 == nV) || (1 == nV && 1 == nW) || (1 == nW && 1 == nU)))
        {
            this->CreateTrackParticle(clusterListAll);
        }
        else
        {
            // TODO - check here whether there is a gap in the 2 in one view when 1:2:0
            // May later choose to resolve simple ambiguities, e.g. of form nU:nV:nW == 1:2:0
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ParticleRecoveryAlgorithm::CheckConsistency(const Cluster *const pClusterU, const Cluster *const pClusterV, const Cluster *const pClusterW) const
{
    // Requirements on X matching
    float xMinU(0.f), xMinV(0.f), xMinW(0.f), xMaxU(0.f), xMaxV(0.f), xMaxW(0.f);
    pClusterU->GetClusterSpanX(xMinU, xMaxU);
    pClusterV->GetClusterSpanX(xMinV, xMaxV);
    pClusterW->GetClusterSpanX(xMinW, xMaxW);

    const float xMin(std::max(xMinU, std::max(xMinV, xMinW)));
    const float xMax(std::min(xMaxU, std::min(xMaxV, xMaxW)));
    const float xOverlap(xMax - xMin);

    if (xOverlap < std::numeric_limits<float>::epsilon())
        return false;

    // Requirements on 3D matching
    float u(std::numeric_limits<float>::max()), v(std::numeric_limits<float>::max()), w(std::numeric_limits<float>::max());

    if ((STATUS_CODE_SUCCESS != LArClusterHelper::GetAverageZ(pClusterU, xMin, xMax, u)) ||
        (STATUS_CODE_SUCCESS != LArClusterHelper::GetAverageZ(pClusterV, xMin, xMax, v)) ||
        (STATUS_CODE_SUCCESS != LArClusterHelper::GetAverageZ(pClusterW, xMin, xMax, w)))
    {
        return false;
    }

    const float uv2w(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, u, v));
    const float vw2u(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_V, TPC_VIEW_W, v, w));
    const float wu2v(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_W, TPC_VIEW_U, w, u));

    const float pseudoChi2(((u - vw2u) * (u - vw2u) + (v - wu2v) * (v - wu2v) + (w - uv2w) * (w - uv2w)) / 3.f);

    if (pseudoChi2 > m_pseudoChi2Cut)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParticleRecoveryAlgorithm::CreateTrackParticle(const ClusterList &clusterList) const
{
    if (clusterList.size() < 2)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const PfoList *pPfoList = NULL;
    std::string pfoListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pPfoList, pfoListName));

    // TODO Correct these placeholder parameters
    PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
    pfoParameters.m_particleId = MU_MINUS;
    pfoParameters.m_charge = PdgTable::GetParticleCharge(pfoParameters.m_particleId.Get());
    pfoParameters.m_mass = PdgTable::GetParticleMass(pfoParameters.m_particleId.Get());
    pfoParameters.m_energy = 0.f;
    pfoParameters.m_momentum = CartesianVector(0.f, 0.f, 0.f);
    pfoParameters.m_clusterList = clusterList;

    const ParticleFlowObject *pPfo(NULL);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pPfo));

    if (!pPfoList->empty())
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, m_outputPfoListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Pfo>(*this, m_outputPfoListName));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

void ParticleRecoveryAlgorithm::SimpleOverlapTensor::AddAssociation(const Cluster *const pCluster1, const Cluster *const pCluster2)
{
    const HitType hitType1(LArClusterHelper::GetClusterHitType(pCluster1));
    const HitType hitType2(LArClusterHelper::GetClusterHitType(pCluster2));

    const Cluster *const pClusterU((TPC_VIEW_U == hitType1) ? pCluster1 : (TPC_VIEW_U == hitType2) ? pCluster2 : NULL);
    const Cluster *const pClusterV((TPC_VIEW_V == hitType1) ? pCluster1 : (TPC_VIEW_V == hitType2) ? pCluster2 : NULL);
    const Cluster *const pClusterW((TPC_VIEW_W == hitType1) ? pCluster1 : (TPC_VIEW_W == hitType2) ? pCluster2 : NULL);

    if (pClusterU && pClusterV && !pClusterW)
    {
        m_clusterNavigationMapUV[pClusterU].push_back(pClusterV);

        if (m_keyClusters.end() == std::find(m_keyClusters.begin(), m_keyClusters.end(), pClusterU))
            m_keyClusters.push_back(pClusterU);
    }
    else if (!pClusterU && pClusterV && pClusterW)
    {
        m_clusterNavigationMapVW[pClusterV].push_back(pClusterW);

        if (m_keyClusters.end() == std::find(m_keyClusters.begin(), m_keyClusters.end(), pClusterV))
            m_keyClusters.push_back(pClusterV);
    }
    else if (pClusterU && !pClusterV && pClusterW)
    {
        m_clusterNavigationMapWU[pClusterW].push_back(pClusterU);

        if (m_keyClusters.end() == std::find(m_keyClusters.begin(), m_keyClusters.end(), pClusterW))
            m_keyClusters.push_back(pClusterW);
    }
    else
    {
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ParticleRecoveryAlgorithm::SimpleOverlapTensor::GetConnectedElements(const Cluster *const pCluster, const bool ignoreUnavailable,
    ClusterList &clusterListU, ClusterList &clusterListV, ClusterList &clusterListW) const
{
    if (ignoreUnavailable && !pCluster->IsAvailable())
        return;

    const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));

    if (!((TPC_VIEW_U == hitType) || (TPC_VIEW_V == hitType) || (TPC_VIEW_W == hitType)))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    ClusterList &clusterList((TPC_VIEW_U == hitType) ? clusterListU : (TPC_VIEW_V == hitType) ? clusterListV : clusterListW);
    const ClusterNavigationMap &navigationMap((TPC_VIEW_U == hitType)   ? m_clusterNavigationMapUV
                                              : (TPC_VIEW_V == hitType) ? m_clusterNavigationMapVW
                                                                        : m_clusterNavigationMapWU);

    if (clusterList.end() != std::find(clusterList.begin(), clusterList.end(), pCluster))
        return;

    clusterList.push_back(pCluster);

    ClusterNavigationMap::const_iterator iter = navigationMap.find(pCluster);

    if (navigationMap.end() == iter)
        return;

    for (ClusterList::const_iterator cIter = iter->second.begin(), cIterEnd = iter->second.end(); cIter != cIterEnd; ++cIter)
        this->GetConnectedElements(*cIter, ignoreUnavailable, clusterListU, clusterListV, clusterListW);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ParticleRecoveryAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputClusterListNames", m_inputClusterListNames));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputPfoListName", m_outputPfoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterCaloHits", m_minClusterCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CheckGaps", m_checkGaps));

    float minClusterLength = std::sqrt(m_minClusterLengthSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterLength", minClusterLength));
    m_minClusterLengthSquared = minClusterLength * minClusterLength;

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterXSpan", m_minClusterXSpan));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexClusterMode", m_vertexClusterMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinVertexLongitudinalDistance", m_minVertexLongitudinalDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxVertexTransverseDistance", m_maxVertexTransverseDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinXOverlapFraction", m_minXOverlapFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinXOverlapFractionGaps", m_minXOverlapFractionGaps));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SampleStepSize", m_sampleStepSize));

    if (m_sampleStepSize < std::numeric_limits<float>::epsilon())
    {
        std::cout << "ParticleRecoveryAlgorithm: Invalid value for SampleStepSize " << m_sampleStepSize << std::endl;
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingFitHalfWindow", m_slidingFitHalfWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PseudoChi2Cut", m_pseudoChi2Cut));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
