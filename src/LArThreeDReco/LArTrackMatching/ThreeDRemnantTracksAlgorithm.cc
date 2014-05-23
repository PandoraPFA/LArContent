/**
 *  @file   LArContent/src/LArThreeDReco/LArTrackMatching/ThreeDRemnantTracksAlgorithm.cc
 *
 *  @brief  Implementation of the three dimensional longitudinal tracks algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArThreeDReco/LArTrackMatching/ThreeDRemnantTracksAlgorithm.h"

#include "LArHelpers/LArGeometryHelper.h"
#include "LArHelpers/LArClusterHelper.h"

using namespace pandora;

namespace lar
{

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDRemnantTracksAlgorithm::CalculateOverlapResult(Cluster *pClusterU, Cluster *pClusterV, Cluster *pClusterW)
{
    // Requirements on overall X overlap
    float xMinU(0.f), xMinV(0.f), xMinW(0.f);
    float xMaxU(0.f), xMaxV(0.f), xMaxW(0.f);
    LArClusterHelper::GetClusterSpanX(pClusterU, xMinU, xMaxU);
    LArClusterHelper::GetClusterSpanX(pClusterV, xMinV, xMaxV);
    LArClusterHelper::GetClusterSpanX(pClusterW, xMinW, xMaxW);

    if (!(xMaxU > xMinV && xMaxV > xMinU) ||
        !(xMaxV > xMinW && xMaxW > xMinV) ||
        !(xMaxW > xMinU && xMaxU > xMinW))
        return;

    const float xSpan(std::max(xMaxU, std::max(xMaxV, xMaxW)) - std::min(xMinU, std::min(xMinV, xMinW)));
    const float xOverlap(std::min(xMaxU, std::min(xMaxV, xMaxW)) - std::max(xMinU, std::max(xMinV, xMinW)));

    if (xOverlap < m_minXOverlap || xOverlap/xSpan < m_minXOverlapFraction)
        return;

    // Match hits with common X coordinates
    CaloHitList clusterHitsU, clusterHitsV, clusterHitsW;
    CaloHitList matchedHitsU, matchedHitsV, matchedHitsW;
    pClusterU->GetOrderedCaloHitList().GetCaloHitList(clusterHitsU);
    pClusterV->GetOrderedCaloHitList().GetCaloHitList(clusterHitsV);
    pClusterW->GetOrderedCaloHitList().GetCaloHitList(clusterHitsW);

    if (clusterHitsU.empty() || clusterHitsV.empty() || clusterHitsW.empty())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    typedef std::map<pandora::CaloHit*, pandora::CaloHitList> HitToHitMap;
    HitToHitMap hitToHitMapUV, hitToHitMapVW, hitToHitMapWU;

    // U-V associations in X coordinate
    for (CaloHitList::const_iterator hIterU = clusterHitsU.begin(), hIterEndU = clusterHitsU.end(); hIterU != hIterEndU; ++hIterU)
    {
        CaloHit *pCaloHitU = *hIterU;

        for (CaloHitList::const_iterator hIterV = clusterHitsV.begin(), hIterEndV = clusterHitsV.end(); hIterV != hIterEndV; ++hIterV)
        {
            CaloHit *pCaloHitV = *hIterV;

            if (std::fabs(pCaloHitU->GetPositionVector().GetX() - pCaloHitV->GetPositionVector().GetX()) >  m_maxXDisplacement)
                continue;

            hitToHitMapUV[pCaloHitU].insert(pCaloHitV);
        }
    }

    // V-W associations in X coordinate
    for (CaloHitList::const_iterator hIterV = clusterHitsV.begin(), hIterEndV = clusterHitsV.end(); hIterV != hIterEndV; ++hIterV)
    {
        CaloHit *pCaloHitV = *hIterV;

        for (CaloHitList::const_iterator hIterW = clusterHitsW.begin(), hIterEndW = clusterHitsW.end(); hIterW != hIterEndW; ++hIterW)
        {
            CaloHit *pCaloHitW = *hIterW;

            if (std::fabs(pCaloHitV->GetPositionVector().GetX() - pCaloHitW->GetPositionVector().GetX()) >  m_maxXDisplacement)
                continue;

            hitToHitMapVW[pCaloHitV].insert(pCaloHitW);
        }
    }

    // W-U associations in X coordinate
    for (CaloHitList::const_iterator hIterW = clusterHitsW.begin(), hIterEndW = clusterHitsW.end(); hIterW != hIterEndW; ++hIterW)
    {
        CaloHit *pCaloHitW = *hIterW;

        for (CaloHitList::const_iterator hIterU = clusterHitsU.begin(), hIterEndU = clusterHitsU.end(); hIterU != hIterEndU; ++hIterU)
        {
            CaloHit *pCaloHitU = *hIterU;

            if (std::fabs(pCaloHitW->GetPositionVector().GetX() - pCaloHitU->GetPositionVector().GetX()) >  m_maxXDisplacement)
                continue;

            hitToHitMapWU[pCaloHitW].insert(pCaloHitU);
        }
    }

    // U-V-W associations in 3D coordinates (TODO: Speed up this calculation...)
    for (HitToHitMap::const_iterator hIterUV = hitToHitMapUV.begin(), hIterEndUV = hitToHitMapUV.end(); hIterUV != hIterEndUV; ++hIterUV)
    {
        CaloHit *pCaloHitU = hIterUV->first;
        const CaloHitList &caloHitListV = hIterUV->second;

        for (CaloHitList::const_iterator hIterV = caloHitListV.begin(), hIterEndV = caloHitListV.end(); hIterV != hIterEndV; ++hIterV)
        {
            CaloHit* pCaloHitV = *hIterV;
            HitToHitMap::const_iterator hIterVW = hitToHitMapVW.find(pCaloHitV);
            if (hitToHitMapVW.end() == hIterVW)
                continue;

            const CaloHitList &caloHitListW = hIterVW->second;

            for (CaloHitList::const_iterator hIterW = caloHitListW.begin(), hIterEndW = caloHitListW.end(); hIterW != hIterEndW; ++hIterW)
            {
                CaloHit* pCaloHitW = *hIterW;
                HitToHitMap::const_iterator hIterWU = hitToHitMapWU.find(pCaloHitW);
                if (hitToHitMapWU.end() == hIterWU)
                    continue;

                const CaloHitList &caloHitListU = hIterWU->second;
                if (!caloHitListU.count(pCaloHitU))
                    continue;

                if (matchedHitsU.count(pCaloHitU) && matchedHitsV.count(pCaloHitV) && matchedHitsW.count(pCaloHitW))
                    continue;

                float chi2(0.f);
                CartesianVector mergedPosition3D(0.f,0.f,0.f);

                LArGeometryHelper::MergeThreePositions3D(TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W,
                    pCaloHitU->GetPositionVector(), pCaloHitV->GetPositionVector(), pCaloHitW->GetPositionVector(),
                    mergedPosition3D, chi2);

                if (chi2 > m_maxMatchedChi2)
                    continue;

                matchedHitsU.insert(pCaloHitU);
                matchedHitsV.insert(pCaloHitV);
                matchedHitsW.insert(pCaloHitW);
            }
        }
    }

    // Calculate matched fractions in each view (use average matched fraction as the overlap result...)
    const float matchedFractionU(static_cast<float>(matchedHitsU.size()) / static_cast<float>(clusterHitsU.size()));
    const float matchedFractionV(static_cast<float>(matchedHitsV.size()) / static_cast<float>(clusterHitsV.size()));
    const float matchedFractionW(static_cast<float>(matchedHitsW.size()) / static_cast<float>(clusterHitsW.size()));

    const float matchedFraction((matchedFractionU + matchedFractionV + matchedFractionW) / 3.f);

    if (matchedFraction < m_minMatchedFraction)
        return;

    m_overlapTensor.SetOverlapResult(pClusterU, pClusterV, pClusterW, matchedFraction);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDRemnantTracksAlgorithm::ExamineTensor()
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

StatusCode ThreeDRemnantTracksAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    AlgorithmToolList algorithmToolList;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle,
        "TrackTools", algorithmToolList));

    for (AlgorithmToolList::const_iterator iter = algorithmToolList.begin(), iterEnd = algorithmToolList.end(); iter != iterEnd; ++iter)
    {
        RemnantTensorTool *pTensorManipulationTool(dynamic_cast<RemnantTensorTool*>(*iter));

        if (NULL == pTensorManipulationTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolList.push_back(pTensorManipulationTool);
    }

    m_nMaxTensorToolRepeats = 5000;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NMaxTensorToolRepeats", m_nMaxTensorToolRepeats));

    m_minXOverlap = 1.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinXOverlap", m_minXOverlap));

    m_minXOverlapFraction = 0.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinXOverlapFraction", m_minXOverlapFraction));

    m_maxXDisplacement = 1.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxXDisplacement", m_maxXDisplacement));

    m_maxMatchedChi2 = 3.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxMatchedChi2", m_maxMatchedChi2));

    m_minMatchedFraction = 0.5f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedFraction", m_minMatchedFraction));

    return ThreeDBaseAlgorithm<float>::ReadSettings(xmlHandle);
}

} // namespace lar
