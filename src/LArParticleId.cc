/**
 *  @file   LArParticleId.cc
 * 
 *  @brief  Implementation of the lar particle id class.
 * 
 *  $Log: $
 */

#include "Helpers/XmlHelper.h"

#include "Objects/Cluster.h"

#include "LArClusterHelper.h"

#include "LArParticleId.h"
#include "LArPseudoLayerCalculator.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace lar
{

using namespace pandora;

bool LArParticleId::LArEmShowerId(const Cluster *const pCluster)
{
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArParticleId::LArPhotonId(const Cluster *const pCluster)
{
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArParticleId::LArElectronId(const Cluster *const pCluster)
{
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArParticleId::LArMuonId(const Cluster *const pCluster)
{
    if (LArClusterHelper::GetLayerOccupancy(pCluster) < 0.75)
        return false;

    if (LArParticleId::LArTrackWidth(pCluster) > 0.5)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArParticleId::LArTrackWidth(const Cluster *const pCluster)
{
    TwoDSlidingXZFitResult twoDSlidingXZFitResult;
    LArParticleId::LArTwoDSlidingXZFit(pCluster, twoDSlidingXZFitResult);

    return twoDSlidingXZFitResult.GetSlidingFitWidth();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArParticleId::LArTwoDSlidingXZFit(const Cluster *const pCluster, TwoDSlidingXZFitResult &twoDSlidingXZFitResult)
{
    const unsigned int innerLayer(pCluster->GetInnerPseudoLayer());
    const unsigned int outerLayer(pCluster->GetOuterPseudoLayer());
    const unsigned int layerFitHalfWindow(LArParticleId::m_layerFitHalfWindow);

    twoDSlidingXZFitResult.m_pCluster = pCluster;
    twoDSlidingXZFitResult.m_layerFitHalfWindow = layerFitHalfWindow;
    LayerFitResultMap &layerFitResultMap(twoDSlidingXZFitResult.m_layerFitResultMap);
    LayerFitContributionMap &layerFitContributionMap(twoDSlidingXZFitResult.m_layerFitContributionMap);

    // Identify fit contributions
    unsigned int slidingNPoints(0);
    double slidingSumX(0.), slidingSumZ(0.), slidingSumXX(0.), slidingSumZX(0.), slidingSumZZ(0.);
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        const unsigned int iLayer(iter->first);
        TwoDSlidingXZFitResult::LayerFitContribution layerFitContribution(iter->second);
        (void) layerFitContributionMap.insert(LayerFitContributionMap::value_type(iLayer, layerFitContribution));

        if (iLayer < innerLayer + layerFitHalfWindow)
        {
            slidingSumX += layerFitContribution.GetSumX();
            slidingSumZ += layerFitContribution.GetSumZ();
            slidingSumXX += layerFitContribution.GetSumXX();
            slidingSumZX += layerFitContribution.GetSumZX();
            slidingSumZZ += layerFitContribution.GetSumZZ();
            slidingNPoints += layerFitContribution.GetNPoints();
        }
    }

    // Sliding fit
    for (unsigned int iLayer = innerLayer; iLayer <= outerLayer; ++iLayer)
    {
        const unsigned int fwdLayer(iLayer + layerFitHalfWindow);
        LayerFitContributionMap::const_iterator fwdIter = layerFitContributionMap.find(fwdLayer);

        if (layerFitContributionMap.end() != fwdIter)
        {
            slidingSumX += fwdIter->second.GetSumX();
            slidingSumZ += fwdIter->second.GetSumZ();
            slidingSumXX += fwdIter->second.GetSumXX();
            slidingSumZX += fwdIter->second.GetSumZX();
            slidingSumZZ += fwdIter->second.GetSumZZ();
            slidingNPoints += fwdIter->second.GetNPoints();
        }

        const unsigned int bwdLayer(iLayer - layerFitHalfWindow - 1);
        LayerFitContributionMap::const_iterator bwdIter = layerFitContributionMap.find(bwdLayer);

        if (layerFitContributionMap.end() != bwdIter)
        {
            slidingSumX -= bwdIter->second.GetSumX();
            slidingSumZ -= bwdIter->second.GetSumZ();
            slidingSumXX -= bwdIter->second.GetSumXX();
            slidingSumZX -= bwdIter->second.GetSumZX();
            slidingSumZZ -= bwdIter->second.GetSumZZ();
            slidingNPoints -= bwdIter->second.GetNPoints();
        }

        if (slidingNPoints > 0)
        {
            const double denominator(slidingSumZZ - slidingSumZ * slidingSumZ / static_cast<double>(slidingNPoints));

            if (std::fabs(denominator) < std::numeric_limits<double>::epsilon())
                continue;

            const double gradient((slidingSumZX - slidingSumZ * slidingSumX / static_cast<double>(slidingNPoints)) / denominator);
            const double intercept((slidingSumZZ * slidingSumX / static_cast<double>(slidingNPoints) - slidingSumZ * slidingSumZX / static_cast<double>(slidingNPoints)) / denominator);

            const double z(LArPseudoLayerCalculator::GetZCoordinate(iLayer));
            const double fitX(intercept + gradient * z);

            const double variance((slidingSumXX - 2. * intercept * slidingSumX - 2. * gradient * slidingSumZX + intercept * intercept * static_cast<double>(slidingNPoints) + 2. * gradient * intercept * slidingSumZ + gradient * gradient * slidingSumZZ) / (1. + gradient * gradient));
            const double rms(std::sqrt(variance / static_cast<double>(slidingNPoints)));

            (void) layerFitResultMap.insert(LayerFitResultMap::value_type(iLayer, TwoDSlidingXZFitResult::LayerFitResult(z, fitX, gradient, rms)));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int LArParticleId::m_layerFitHalfWindow = 20;
float LArParticleId::m_trackFitMaxRms = 0.25; // cm
float LArParticleId::m_minCosScatteringAngle = std::cos(M_PI * 20.f / 180.f); // radians

StatusCode LArParticleId::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "LayerFitHalfWindow", m_layerFitHalfWindow));    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TrackFitMaxRms", m_trackFitMaxRms));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCosScatteringAngle", m_minCosScatteringAngle));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

  LArParticleId::TwoDSlidingXZFitResult::LayerFitResult::LayerFitResult(const double z, const double fitX, const double gradient, const double rms) :
    m_z(z),
    m_fitX(fitX),
    m_gradient(gradient),
    m_rms(rms)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LArParticleId::TwoDSlidingXZFitResult::LayerFitContribution::LayerFitContribution(const CaloHitList *const pCaloHitList) :
    m_sumX(0.),
    m_sumZ(0.),
    m_sumXX(0.),
    m_sumZX(0.),
    m_sumZZ(0.),
    m_nPoints(0)
{
    for (CaloHitList::const_iterator iter = pCaloHitList->begin(), iterEnd = pCaloHitList->end(); iter != iterEnd; ++iter)
    {
        const CaloHit *pCaloHit = *iter;
        const float x(pCaloHit->GetPositionVector().GetX());
        const float z(pCaloHit->GetPositionVector().GetZ());
        m_sumX += x;
        m_sumZ += z;
        m_sumXX += x * x;
        m_sumZX += z * x;
        m_sumZZ += z * z;
        ++m_nPoints;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

float LArParticleId::TwoDSlidingXZFitResult::GetSlidingFitWidth() const
{
    FloatVector residuals;
    const OrderedCaloHitList &orderedCaloHitList(m_pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(); iter != orderedCaloHitList.end(); ++iter)
    {
        const unsigned int layer(iter->first);
        LayerFitResultMap::const_iterator fitResultIter = m_layerFitResultMap.find(layer);

        if (m_layerFitResultMap.end() == fitResultIter)
            continue;

        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            CaloHit *pCaloHit = *hitIter;
            const double x(pCaloHit->GetPositionVector().GetX());
            const double fitX(fitResultIter->second.GetFitX());
            const double gradient(fitResultIter->second.GetGradient());
            const double residualSquared((fitX - x) * (fitX - x) / (1. + gradient * gradient)); // angular correction (note: this is cheating!)
            residuals.push_back(residualSquared);
        }
    }

    if (residuals.empty())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    std::sort(residuals.begin(), residuals.end());
    static const float m_trackResidualQuantile(0.8f);
    const float theQuantile(residuals[m_trackResidualQuantile * residuals.size()]);
    return std::sqrt(theQuantile);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LArParticleId::TwoDSlidingXZFitResult::FindLargestScatter(unsigned int &largestScatterLayer) const
{
    static const float m_trackFitMaxRms(0.25); // cm
    static const float m_minCosScatteringAngle(std::cos(M_PI * 20.f / 180.f)); // radians

    // Bail out if track is too short
    const unsigned int minLayer(m_pCluster->GetInnerPseudoLayer());
    const unsigned int nClusterLayers(1 + m_pCluster->GetOuterPseudoLayer() - m_pCluster->GetInnerPseudoLayer());

    if (nClusterLayers <= 2 * m_layerFitHalfWindow)
        return STATUS_CODE_NOT_FOUND;

    // Find point of largest scatter
    unsigned int splitLayer(0);
    double splitCosTheta(m_minCosScatteringAngle);

    const OrderedCaloHitList &orderedCaloHitList(m_pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(); iter != orderedCaloHitList.end(); ++iter)
    {
        const unsigned int layer(iter->first);

        if (layer - minLayer >= nClusterLayers - 2 * m_layerFitHalfWindow)
            break;

        LayerFitResultMap::const_iterator fitIter1 = m_layerFitResultMap.find(layer);
        LayerFitResultMap::const_iterator fitIter2 = m_layerFitResultMap.find(layer + 2 * m_layerFitHalfWindow);

        if ((m_layerFitResultMap.end() == fitIter1) || (m_layerFitResultMap.end() == fitIter2))
            continue;

        const double r1(fitIter1->second.GetRms());
        const double r2(fitIter2->second.GetRms());
        const double m1(fitIter1->second.GetGradient());
        const double m2(fitIter2->second.GetGradient());

        const double cosTheta = (1. + m1 * m2) / (std::sqrt(1.+ m1 * m1) * std::sqrt(1. + m2 * m2));

        if ((r1< m_trackFitMaxRms) && (r2 < m_trackFitMaxRms) && (cosTheta < splitCosTheta))
        {
            splitCosTheta = cosTheta;
            splitLayer = layer + m_layerFitHalfWindow;
        }
    }

    if ((splitLayer <= m_pCluster->GetInnerPseudoLayer()) || (splitLayer >= m_pCluster->GetOuterPseudoLayer()))
        return STATUS_CODE_NOT_FOUND;

    largestScatterLayer = splitLayer;
    return STATUS_CODE_SUCCESS;
}

} // namespace lar
