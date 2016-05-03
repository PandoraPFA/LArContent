/**
 *  @file   LArContent/src/LArPlugins/LArParticleIdPlugins.cc
 * 
 *  @brief  Implementation of the lar particle id plugins class.
 * 
 *  $Log: $
 */

#include "Helpers/XmlHelper.h"

#include "Objects/Cluster.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"

#include "LArPlugins/LArParticleIdPlugins.h"

#include "LArObjects/LArTwoDSlidingFitResult.h"

#include <algorithm>
#include <cmath>

namespace lar_content
{

using namespace pandora;

LArParticleIdPlugins::LArMuonId::LArMuonId() :
    m_layerFitHalfWindow(20),
    m_minLayerOccupancy (0.75f),
    m_maxTrackWidth(0.5f),
    m_trackResidualQuantile(0.8f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArParticleIdPlugins::LArMuonId::IsMatch(const Cluster *const pCluster) const
{
    if (LArClusterHelper::GetLayerOccupancy(pCluster) < m_minLayerOccupancy)
        return false;

    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const TwoDSlidingFitResult twoDSlidingFitResult(pCluster, m_layerFitHalfWindow, slidingFitPitch);

    if (this->GetMuonTrackWidth(twoDSlidingFitResult) > m_maxTrackWidth)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArParticleIdPlugins::LArMuonId::GetMuonTrackWidth(const TwoDSlidingFitResult &twoDSlidingFitResult) const
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
    const float theQuantile(residuals[m_trackResidualQuantile * residuals.size()]);

    return std::sqrt(theQuantile);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LArParticleIdPlugins::LArMuonId::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "LayerFitHalfWindow", m_layerFitHalfWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinLayerOccupancy", m_minLayerOccupancy));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTrackWidth", m_maxTrackWidth));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TrackResidualQuantile", m_trackResidualQuantile));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
