/**
 *  @file   LArContent/src/LArHelpers/LArParticleIdHelper.cc
 * 
 *  @brief  Implementation of the lar particle id class.
 * 
 *  $Log: $
 */

#include "Helpers/XmlHelper.h"

#include "Objects/Cluster.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArParticleIdHelper.h"

#include <algorithm>
#include <cmath>

namespace lar
{

using namespace pandora;

bool LArParticleIdHelper::LArEmShowerId(const Cluster *const /*pCluster*/)
{
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArParticleIdHelper::LArPhotonId(const Cluster *const /*pCluster*/)
{
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArParticleIdHelper::LArElectronId(const Cluster *const /*pCluster*/)
{
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArParticleIdHelper::LArMuonId(const Cluster *const pCluster)
{
    if (LArClusterHelper::GetLayerOccupancy(pCluster) < m_muonIdMinLayerOccupancy)
        return false;

    const TwoDSlidingFitResult twoDSlidingFitResult(pCluster, m_muonIdLayerFitHalfWindow);

    if (LArParticleIdHelper::GetMuonTrackWidth(twoDSlidingFitResult) > m_muonIdMaxTrackWidth)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArParticleIdHelper::GetMuonTrackWidth(const TwoDSlidingFitResult &twoDSlidingFitResult)
{
    FloatVector residuals;
    const OrderedCaloHitList &orderedCaloHitList(twoDSlidingFitResult.GetCluster()->GetOrderedCaloHitList());
    const LayerFitResultMap &layerFitResultMap(twoDSlidingFitResult.GetLayerFitResultMap());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(); iter != orderedCaloHitList.end(); ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            float rL(0.f), rT(0.f);
            twoDSlidingFitResult.GetLocalPosition((*hitIter)->GetPositionVector(), rL, rT);
            const int layer(twoDSlidingFitResult.GetLayer(rL));

            LayerFitResultMap::const_iterator fitResultIter = layerFitResultMap.find(layer);

            if (layerFitResultMap.end() == fitResultIter)
                continue;

            const double fitT(fitResultIter->second.GetFitT());
            const double gradient(fitResultIter->second.GetGradient());
            const double residualSquared((fitT - rT) * (fitT - rT) / (1. + gradient * gradient)); // angular correction (note: this is cheating!)
            residuals.push_back(residualSquared);
        }
    }

    if (residuals.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    std::sort(residuals.begin(), residuals.end());
    const float theQuantile(residuals[m_muonIdTrackResidualQuantile * residuals.size()]);

    return std::sqrt(theQuantile);
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int LArParticleIdHelper::m_muonIdLayerFitHalfWindow = 20;
float LArParticleIdHelper::m_muonIdMinLayerOccupancy = 0.75f;
float LArParticleIdHelper::m_muonIdMaxTrackWidth = 0.5f;
float LArParticleIdHelper::m_muonIdTrackResidualQuantile = 0.8f;

StatusCode LArParticleIdHelper::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MuonIdLayerFitHalfWindow", m_muonIdLayerFitHalfWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MuonIdMinLayerOccupancy", m_muonIdMinLayerOccupancy));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MuonIdMaxTrackWidth", m_muonIdMaxTrackWidth));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MuonIdTrackResidualQuantile", m_muonIdTrackResidualQuantile));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
