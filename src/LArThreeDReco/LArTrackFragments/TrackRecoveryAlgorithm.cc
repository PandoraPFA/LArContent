/**
 *  @file   LArContent/src/LArThreeDReco/LArTrackFragments/TrackRecoveryAlgorithm.cc
 * 
 *  @brief  Implementation of the track recovery algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"

#include "LArThreeDReco/LArTrackFragments/TrackRecoveryAlgorithm.h"

using namespace pandora;

namespace lar_content
{

TrackRecoveryAlgorithm::TrackRecoveryAlgorithm() :
    m_minClusterCaloHits(5),
    m_minClusterLengthSquared(1.f),
    m_minClusterXSpan(0.5f),
    m_minXOverlapFraction(0.5f),
    m_pseudoChi2Cut(20.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackRecoveryAlgorithm::Run()
{
    ClusterList selectedClusterListU, selectedClusterListV, selectedClusterListW;
    this->SelectInputClusters(m_inputClusterListNameU, selectedClusterListU);
    this->SelectInputClusters(m_inputClusterListNameV, selectedClusterListV);
    this->SelectInputClusters(m_inputClusterListNameW, selectedClusterListW);

    SimpleOverlapTensor overlapTensor;
    this->FindOverlaps(selectedClusterListU, selectedClusterListV, overlapTensor);
    this->FindOverlaps(selectedClusterListV, selectedClusterListW, overlapTensor);
    this->FindOverlaps(selectedClusterListW, selectedClusterListU, overlapTensor);

    this->ExamineTensor(overlapTensor);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackRecoveryAlgorithm::SelectInputClusters(const std::string &inputClusterListName, ClusterList &selectedClusterList) const
{
    const ClusterList *pInputClusterList(NULL);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this,
        inputClusterListName, pInputClusterList));

    if (!pInputClusterList)
        return;

    for (ClusterList::const_iterator iter = pInputClusterList->begin(), iterEnd = pInputClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (!pCluster->IsAvailable())
            continue;

        if (pCluster->GetNCaloHits() < m_minClusterCaloHits)
            continue;

        if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLengthSquared)
            continue;

        float xMin(0.f), xMax(0.f);
        LArClusterHelper::GetClusterSpanX(pCluster, xMin, xMax);

        if ((xMax - xMin) < m_minClusterXSpan)
            continue;

        selectedClusterList.insert(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackRecoveryAlgorithm::FindOverlaps(const ClusterList &clusterList1, const ClusterList &clusterList2, SimpleOverlapTensor &overlapTensor) const
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

bool TrackRecoveryAlgorithm::IsOverlap(const Cluster *const pCluster1, const Cluster *const pCluster2) const
{
    if (LArClusterHelper::GetClusterHitType(pCluster1) == LArClusterHelper::GetClusterHitType(pCluster2))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    if ((0 == pCluster1->GetNCaloHits()) || (0 == pCluster2->GetNCaloHits()))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    float xMin1(0.f), xMax1(0.f), xMin2(0.f), xMax2(0.f);
    LArClusterHelper::GetClusterSpanX(pCluster1, xMin1, xMax1);
    LArClusterHelper::GetClusterSpanX(pCluster2, xMin2, xMax2);

    const float xSpan1(xMax1 - xMin1), xSpan2(xMax2 - xMin2);

    if ((xSpan1 < std::numeric_limits<float>::epsilon()) || (xSpan2 < std::numeric_limits<float>::epsilon()))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const float xOverlap(std::min(xMax1, xMax2) - std::max(xMin1, xMin2));

    if (((xOverlap / xSpan1) < m_minXOverlapFraction) || ((xOverlap / xSpan2) < m_minXOverlapFraction))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackRecoveryAlgorithm::ExamineTensor(const SimpleOverlapTensor &overlapTensor) const
{
    for (ClusterList::const_iterator iter = overlapTensor.GetKeyClusters().begin(), iterEnd = overlapTensor.GetKeyClusters().end(); iter != iterEnd; ++iter)
    {
        ClusterList clusterListU, clusterListV, clusterListW;

        overlapTensor.GetConnectedElements(*iter, true, clusterListU, clusterListV, clusterListW);
        const unsigned int nU(clusterListU.size()), nV(clusterListV.size()), nW(clusterListW.size());

        if ((0 == nU * nV) && (0 == nV * nW) && (0 == nW * nU))
            continue;

        ClusterList clusterListAll;
        clusterListAll.insert(clusterListU.begin(), clusterListU.end());
        clusterListAll.insert(clusterListV.begin(), clusterListV.end());
        clusterListAll.insert(clusterListW.begin(), clusterListW.end());

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
            // May later choose to resolve simple ambiguities, e.g. of form nU:nV:nW == 1:2:0
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackRecoveryAlgorithm::CheckConsistency(const Cluster *const pClusterU, const Cluster *const pClusterV, const Cluster *const pClusterW) const
{
    try
    {
        // Requirements on X matching
        float xMinU(0.f), xMinV(0.f), xMinW(0.f), xMaxU(0.f), xMaxV(0.f), xMaxW(0.f);
        LArClusterHelper::GetClusterSpanX(pClusterU, xMinU, xMaxU);
        LArClusterHelper::GetClusterSpanX(pClusterV, xMinV, xMaxV);
        LArClusterHelper::GetClusterSpanX(pClusterW, xMinW, xMaxW);

        const float xMin(std::max(xMinU, std::max(xMinV, xMinW)));
        const float xMax(std::min(xMaxU, std::min(xMaxV, xMaxW)));
        const float xOverlap(xMax - xMin);

        if (xOverlap < std::numeric_limits<float>::epsilon())
            return false;

        // Requirements on 3D matching
        const float u(LArClusterHelper::GetAverageZ(pClusterU, xMin, xMax));
        const float v(LArClusterHelper::GetAverageZ(pClusterV, xMin, xMax));
        const float w(LArClusterHelper::GetAverageZ(pClusterW, xMin, xMax));

        const float uv2w(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, u, v));
        const float vw2u(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_V, TPC_VIEW_W, v, w));
        const float wu2v(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_W, TPC_VIEW_U, w, u));

        const float pseudoChi2(((u - vw2u) * (u - vw2u) + (v - wu2v) * (v - wu2v) + (w - uv2w) * (w - uv2w)) / 3.f);

        if (pseudoChi2 > m_pseudoChi2Cut)
            return false;

        return true;
    }
    catch(StatusCodeException &)
    {
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackRecoveryAlgorithm::CreateTrackParticle(const ClusterList &clusterList) const
{
    if (clusterList.size() < 2)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const PfoList *pPfoList = NULL; std::string pfoListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pPfoList, pfoListName));

    // TODO Correct these placeholder parameters
    PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
    pfoParameters.m_particleId = MU_MINUS; // Track
    pfoParameters.m_charge = PdgTable::GetParticleCharge(pfoParameters.m_particleId.Get());
    pfoParameters.m_mass = PdgTable::GetParticleMass(pfoParameters.m_particleId.Get());
    pfoParameters.m_energy = 0.f;
    pfoParameters.m_momentum = CartesianVector(0.f, 0.f, 0.f);
    pfoParameters.m_clusterList = clusterList;

    ParticleFlowObject *pPfo(NULL);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pPfo));

    if (!pPfoList->empty())
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, m_outputPfoListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Pfo>(*this, m_outputPfoListName));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

void TrackRecoveryAlgorithm::SimpleOverlapTensor::AddAssociation(Cluster *pCluster1, Cluster *pCluster2)
{
    const HitType hitType1(LArClusterHelper::GetClusterHitType(pCluster1));
    const HitType hitType2(LArClusterHelper::GetClusterHitType(pCluster2));

    Cluster *pClusterU((TPC_VIEW_U == hitType1) ? pCluster1 : (TPC_VIEW_U == hitType2) ? pCluster2 : NULL);
    Cluster *pClusterV((TPC_VIEW_V == hitType1) ? pCluster1 : (TPC_VIEW_V == hitType2) ? pCluster2 : NULL);
    Cluster *pClusterW((TPC_VIEW_W == hitType1) ? pCluster1 : (TPC_VIEW_W == hitType2) ? pCluster2 : NULL);

    if (pClusterU && pClusterV && !pClusterW)
    {
        m_clusterNavigationMapUV[pClusterU].insert(pClusterV);
        m_keyClusters.insert(pClusterU);
    }
    else if (!pClusterU && pClusterV && pClusterW)
    {
        m_clusterNavigationMapVW[pClusterV].insert(pClusterW);
        m_keyClusters.insert(pClusterV);
    }
    else if (pClusterU && !pClusterV && pClusterW)
    {
        m_clusterNavigationMapWU[pClusterW].insert(pClusterU);
        m_keyClusters.insert(pClusterW);
    }
    else
    {
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackRecoveryAlgorithm::SimpleOverlapTensor::GetConnectedElements(Cluster *const pCluster, const bool ignoreUnavailable,
    ClusterList &clusterListU, ClusterList &clusterListV, ClusterList &clusterListW) const
{
    if (ignoreUnavailable && !pCluster->IsAvailable())
        return;

    const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));

    if (!((TPC_VIEW_U == hitType) || (TPC_VIEW_V == hitType) || (TPC_VIEW_W == hitType)))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    ClusterList &clusterList((TPC_VIEW_U == hitType) ? clusterListU : (TPC_VIEW_V == hitType) ? clusterListV : clusterListW);
    const ClusterNavigationMap &navigationMap((TPC_VIEW_U == hitType) ? m_clusterNavigationMapUV : (TPC_VIEW_V == hitType) ? m_clusterNavigationMapVW : m_clusterNavigationMapWU);

    if (!clusterList.insert(pCluster).second)
        return;

    ClusterNavigationMap::const_iterator iter = navigationMap.find(pCluster);

    if (navigationMap.end() == iter)
        return;

    for (ClusterList::const_iterator cIter = iter->second.begin(), cIterEnd = iter->second.end(); cIter != cIterEnd; ++cIter)
        this->GetConnectedElements(*cIter, ignoreUnavailable, clusterListU, clusterListV, clusterListW);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackRecoveryAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameU", m_inputClusterListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameV", m_inputClusterListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameW", m_inputClusterListNameW));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputPfoListName", m_outputPfoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterCaloHits", m_minClusterCaloHits));

    float minClusterLength = std::sqrt(m_minClusterLengthSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLength", minClusterLength));
    m_minClusterLengthSquared = minClusterLength * minClusterLength;

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterXSpan", m_minClusterXSpan));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinXOverlapFraction", m_minXOverlapFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PseudoChi2Cut", m_pseudoChi2Cut));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
