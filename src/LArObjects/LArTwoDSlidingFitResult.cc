/**
 *  @file   LArContent/src/LArObjects/LArTwoDSlidingFitResult.cc
 * 
 *  @brief  Implementation of the lar two dimensional sliding fit result class.
 * 
 *  $Log: $
 */

#include "LArCalculators/LArPseudoLayerCalculator.h"

#include "LArHelpers/LArGeometryHelper.h"

#include "LArObjects/LArTwoDSlidingFitResult.h"

#include "Objects/Cluster.h"

#include <algorithm>
#include <cmath>
#include <limits>

using namespace pandora;

namespace lar
{

TwoDSlidingFitResult::TwoDSlidingFitResult() :
    m_pCluster(NULL),
    m_layerFitHalfWindow(0),
    m_axisIntercept(0.f, 0.f, 0.f),
    m_axisDirection(0.f, 0.f, 0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitResult::GetMinAndMaxX(float &minX, float &maxX) const
{
    CartesianVector startLayerPosition(0.f, 0.f, 0.f);
    LayerFitResultMap::const_iterator startLayerIter = m_layerFitResultMap.begin();
    this->GetGlobalPosition(startLayerIter->second.GetL(), startLayerIter->second.GetFitT(), startLayerPosition);

    CartesianVector endLayerPosition(0.f, 0.f, 0.f);
    LayerFitResultMap::const_reverse_iterator endLayerIter = m_layerFitResultMap.rbegin();
    this->GetGlobalPosition(endLayerIter->second.GetL(), endLayerIter->second.GetFitT(), endLayerPosition);

    minX = std::min(startLayerPosition.GetX(), endLayerPosition.GetX());
    maxX = std::max(startLayerPosition.GetX(), endLayerPosition.GetX());
}

//------------------------------------------------------------------------------------------------------------------------------------------

int TwoDSlidingFitResult::GetLayer(const float rL) const
{
    return std::floor(rL / LArGeometryHelper::GetLArPseudoLayerCalculator()->GetZPitch());
}

//------------------------------------------------------------------------------------------------------------------------------------------

float TwoDSlidingFitResult::GetL(const int layer) const
{
    return (static_cast<float>(layer) + 0.5f) * LArGeometryHelper::GetLArPseudoLayerCalculator()->GetZPitch();
}

//------------------------------------------------------------------------------------------------------------------------------------------

float TwoDSlidingFitResult::GetLayerFitHalfWindowLength() const
{
    return (static_cast<float>(m_layerFitHalfWindow)) * LArGeometryHelper::GetLArPseudoLayerCalculator()->GetZPitch();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitResult::GetLocalPosition(const CartesianVector &position, float &rL, float &rT) const
{
    const CartesianVector displacement(position - m_axisIntercept);
    const CartesianVector crossProduct(displacement.GetCrossProduct(m_axisDirection));

    rL = displacement.GetDotProduct(m_axisDirection);
    rT = (crossProduct.GetY() < 0.f) ? (-1.f * crossProduct.GetMagnitude()) : crossProduct.GetMagnitude();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitResult::GetLocalDirection(const CartesianVector &direction, float&dTdL) const
{
    float pL(0.f), pT(0.f);
    this->GetLocalPosition((direction + m_axisIntercept), pL, pT);

    if (std::fabs(pL) < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    dTdL = pT / pL;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitResult::GetGlobalPosition(const float rL, const float rT, CartesianVector &position) const
{
    const CartesianVector positiveTDirection(m_axisDirection.GetCrossProduct(CartesianVector(0.f, 1.f, 0.f)));
    position = m_axisIntercept + (m_axisDirection * rL) + (positiveTDirection * rT);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitResult::GetGlobalDirection(const float dTdL, CartesianVector &direction) const
{
    const float pL(1.f / std::sqrt(1.f + dTdL * dTdL));
    const float pT(dTdL / std::sqrt(1.f + dTdL * dTdL));

    CartesianVector globalCoordinates(0.f, 0.f, 0.f);
    this->GetGlobalPosition(pL, pT, globalCoordinates);
    direction = (globalCoordinates - m_axisIntercept).GetUnitVector();
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector TwoDSlidingFitResult::GetGlobalMinLayerPosition() const
{
    if (m_layerFitResultMap.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    LayerFitResultMap::const_iterator iter = m_layerFitResultMap.begin();
    CartesianVector position(0.f, 0.f, 0.f);
    this->GetGlobalPosition(iter->second.GetL(), iter->second.GetFitT(), position);
    return position;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector TwoDSlidingFitResult::GetGlobalMaxLayerPosition() const
{
    if (m_layerFitResultMap.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    LayerFitResultMap::const_reverse_iterator iter = m_layerFitResultMap.rbegin();
    CartesianVector position(0.f, 0.f, 0.f);
    this->GetGlobalPosition(iter->second.GetL(), iter->second.GetFitT(), position);
    return position;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector TwoDSlidingFitResult::GetGlobalMinLayerDirection() const
{
    if (m_layerFitResultMap.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    LayerFitResultMap::const_iterator iter = m_layerFitResultMap.begin();
    CartesianVector direction(0.f, 0.f, 0.f);
    this->GetGlobalDirection(iter->second.GetGradient(), direction);
    return direction;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector TwoDSlidingFitResult::GetGlobalMaxLayerDirection() const
{
    if (m_layerFitResultMap.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    LayerFitResultMap::const_reverse_iterator iter = m_layerFitResultMap.rbegin();
    CartesianVector direction(0.f, 0.f, 0.f);
    this->GetGlobalDirection(iter->second.GetGradient(), direction);
    return direction;
}
 
//------------------------------------------------------------------------------------------------------------------------------------------

float TwoDSlidingFitResult::GetMinLayerRms() const
{
    if (m_layerFitResultMap.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    LayerFitResultMap::const_iterator iter = m_layerFitResultMap.begin();
    return iter->second.GetRms();
}

//------------------------------------------------------------------------------------------------------------------------------------------

float TwoDSlidingFitResult::GetMaxLayerRms() const
{
    if (m_layerFitResultMap.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    LayerFitResultMap::const_reverse_iterator iter = m_layerFitResultMap.rbegin();
    return iter->second.GetRms();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitResult::GetGlobalFitPosition(const float rL, CartesianVector &position) const
{
    float firstWeight(0.f), secondWeight(0.f);
    LayerFitResultMap::const_iterator firstLayerIter, secondLayerIter;
    this->GetSurroundingLayerIterators(rL, firstLayerIter, secondLayerIter);
    this->GetLayerInterpolationWeights(rL, firstLayerIter, secondLayerIter, firstWeight, secondWeight);
    this->GetGlobalFitInterpolatedPosition(firstLayerIter, secondLayerIter, firstWeight, secondWeight, position);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitResult::GetGlobalFitDirection(const float rL, CartesianVector &direction) const
{
    float firstWeight(0.f), secondWeight(0.f);
    LayerFitResultMap::const_iterator firstLayerIter, secondLayerIter;
    this->GetSurroundingLayerIterators(rL, firstLayerIter, secondLayerIter);
    this->GetLayerInterpolationWeights(rL, firstLayerIter, secondLayerIter, firstWeight, secondWeight);
    this->GetGlobalFitInterpolatedDirection(firstLayerIter, secondLayerIter, firstWeight, secondWeight, direction);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float TwoDSlidingFitResult::GetFitRms(const float rL) const
{
    float rms(std::numeric_limits<float>::max());
    float firstWeight(0.f), secondWeight(0.f);
    LayerFitResultMap::const_iterator firstLayerIter, secondLayerIter;
    this->GetSurroundingLayerIterators(rL, firstLayerIter, secondLayerIter);
    this->GetLayerInterpolationWeights(rL, firstLayerIter, secondLayerIter, firstWeight, secondWeight);
    this->GetFitInterpolatedRms(firstLayerIter, secondLayerIter, firstWeight, secondWeight, rms);
    return rms;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitResult::GetGlobalFitPosition(const float p, const bool useX, CartesianVector &position) const
{
    float firstWeight(0.f), secondWeight(0.f);
    LayerFitResultMap::const_iterator firstLayerIter, secondLayerIter;
    this->GetSurroundingLayerIterators(p, useX, firstLayerIter, secondLayerIter);
    this->GetLayerInterpolationWeights(p, useX, firstLayerIter, secondLayerIter, firstWeight, secondWeight);
    this->GetGlobalFitInterpolatedPosition(firstLayerIter, secondLayerIter, firstWeight, secondWeight, position);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitResult::GetGlobalFitDirection(const float p, const bool useX, CartesianVector &direction) const
{
    float firstWeight(0.f), secondWeight(0.f);
    LayerFitResultMap::const_iterator firstLayerIter, secondLayerIter;
    this->GetSurroundingLayerIterators(p, useX, firstLayerIter, secondLayerIter);
    this->GetLayerInterpolationWeights(p, useX, firstLayerIter, secondLayerIter, firstWeight, secondWeight);
    this->GetGlobalFitInterpolatedDirection(firstLayerIter, secondLayerIter, firstWeight, secondWeight, direction);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitResult::GetGlobalFitProjection(const CartesianVector &inputPosition, CartesianVector &projectedPosition) const
{
    float rL(0.f), rT(0.f);
    this->GetLocalPosition(inputPosition, rL, rT); 
    this->GetGlobalFitPosition(rL, projectedPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float TwoDSlidingFitResult::GetCosScatteringAngle(const float rL) const
{
    const float halfWindowLength(this->GetLayerFitHalfWindowLength());

    CartesianVector firstDirection(0.f,0.f,0.f);
    CartesianVector secondDirection(0.f,0.f,0.f);

    this->GetGlobalFitDirection(rL - halfWindowLength, firstDirection);
    this->GetGlobalFitDirection(rL + halfWindowLength, secondDirection);

    return firstDirection.GetDotProduct(secondDirection);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitResult::GetGlobalFitInterpolatedPosition(const LayerFitResultMap::const_iterator &firstLayerIter,
    const LayerFitResultMap::const_iterator &secondLayerIter, const float &firstWeight, const float &secondWeight,
    CartesianVector &position) const
{
    if (m_layerFitResultMap.end() == firstLayerIter || m_layerFitResultMap.end() == secondLayerIter)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    CartesianVector firstLayerPosition(0.f, 0.f, 0.f);
    this->GetGlobalPosition(firstLayerIter->second.GetL(), firstLayerIter->second.GetFitT(), firstLayerPosition);

    if (firstLayerIter == secondLayerIter)
    {
        position = firstLayerPosition;
        return;
    }

    CartesianVector secondLayerPosition(0.f, 0.f, 0.f);
    this->GetGlobalPosition(secondLayerIter->second.GetL(), secondLayerIter->second.GetFitT(), secondLayerPosition);

    if (firstWeight + secondWeight < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    position = (firstLayerPosition * firstWeight + secondLayerPosition * secondWeight) * (1.f / (firstWeight + secondWeight));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitResult::GetGlobalFitInterpolatedDirection(const LayerFitResultMap::const_iterator &firstLayerIter,
    const LayerFitResultMap::const_iterator &secondLayerIter, const float &firstWeight, const float &secondWeight,
    CartesianVector &direction) const
{   
    if (m_layerFitResultMap.end() == firstLayerIter || m_layerFitResultMap.end() == secondLayerIter)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    CartesianVector firstLayerDirection(0.f,0.f,0.f);
    this->GetGlobalDirection(firstLayerIter->second.GetGradient(),firstLayerDirection);

    if (firstLayerIter == secondLayerIter)
    {
        direction = firstLayerDirection;
        return;
    }

    CartesianVector secondLayerDirection(0.f,0.f,0.f);
    this->GetGlobalDirection(secondLayerIter->second.GetGradient(),secondLayerDirection);

    if (firstWeight + secondWeight < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    direction = (firstLayerDirection * firstWeight + secondLayerDirection * secondWeight).GetUnitVector();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitResult::GetFitInterpolatedRms(const LayerFitResultMap::const_iterator &firstLayerIter,
    const LayerFitResultMap::const_iterator &secondLayerIter, const float &firstWeight, const float &secondWeight, float &rms) const
{   
    if (m_layerFitResultMap.end() == firstLayerIter || m_layerFitResultMap.end() == secondLayerIter)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const float firstLayerRms(firstLayerIter->second.GetRms());

    if (firstLayerIter == secondLayerIter)
    {
        rms = firstLayerRms;
        return;
    }

    const float secondLayerRms(secondLayerIter->second.GetRms());

    if (firstWeight + secondWeight < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    rms = (firstLayerRms * firstWeight + secondLayerRms * secondWeight) / (firstWeight + secondWeight);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitResult::GetSurroundingLayerIterators(const float rL, LayerFitResultMap::const_iterator &firstLayerIter,
    LayerFitResultMap::const_iterator &secondLayerIter) const
{
    if (m_layerFitResultMap.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED); 

    const int minLayer(m_layerFitResultMap.begin()->first), maxLayer(m_layerFitResultMap.rbegin()->first);
    const int thisLayer(this->GetLayer(rL));
    const int startLayer((thisLayer >= maxLayer) ? thisLayer - 1 : thisLayer);

    // Bail out if the layer is out of range
    if ((startLayer < minLayer) || (startLayer >= maxLayer))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    // First layer iterator
    firstLayerIter = m_layerFitResultMap.end();

    for (int iLayer = startLayer; iLayer >= minLayer; --iLayer)
    {
        firstLayerIter = m_layerFitResultMap.find(iLayer);

        if (m_layerFitResultMap.end() != firstLayerIter)
            break;
    }

    if (m_layerFitResultMap.end() == firstLayerIter)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    // Second layer iterator
    secondLayerIter = m_layerFitResultMap.end();

    for (int iLayer = startLayer + 1; iLayer <= maxLayer; ++iLayer)
    {
        secondLayerIter = m_layerFitResultMap.find(iLayer);

        if (m_layerFitResultMap.end() != secondLayerIter)
            break;
    }

    if (m_layerFitResultMap.end() == secondLayerIter)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND); 
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitResult::GetSurroundingLayerIterators(const float p, const bool useX, LayerFitResultMap::const_iterator &firstLayerIter,
    LayerFitResultMap::const_iterator &secondLayerIter) const
{
    if (m_layerFitResultMap.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    if (useX && (std::fabs(m_axisDirection.GetX()) < std::numeric_limits<float>::epsilon()))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    if (!useX && (std::fabs(m_axisDirection.GetZ()) < std::numeric_limits<float>::epsilon()))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    // Find start layer
    const float startL(useX ? (p - m_axisIntercept.GetX()) / m_axisDirection.GetX() : (p - m_axisIntercept.GetZ()) / m_axisDirection.GetZ());
    const int minLayer(m_layerFitResultMap.begin()->first), maxLayer(m_layerFitResultMap.rbegin()->first);
    const int startLayer(std::max(minLayer, std::min(maxLayer, this->GetLayer(startL))));

    // Find nearest layer iterator to start layer
    LayerFitResultMap::const_iterator startLayerIter = m_layerFitResultMap.end();
    CartesianVector startLayerPosition(0.f, 0.f, 0.f);

    for (int iLayer = startLayer; iLayer <= maxLayer; ++iLayer)
    {
        startLayerIter = m_layerFitResultMap.find(iLayer);

        if (m_layerFitResultMap.end() != startLayerIter)
            break;
    }

    if (m_layerFitResultMap.end() == startLayerIter)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    this->GetGlobalPosition(startLayerIter->second.GetL(), startLayerIter->second.GetFitT(), startLayerPosition);

    const bool startIsAhead(useX ? ((startLayerPosition.GetX() - p) > std::numeric_limits<float>::epsilon()) :
        ((startLayerPosition.GetZ() - p) > std::numeric_limits<float>::epsilon()));
    const bool increasesWithLayers(useX ? (m_axisDirection.GetX() > std::numeric_limits<float>::epsilon()) :
        (m_axisDirection.GetZ() > std::numeric_limits<float>::epsilon()));
    const int increment = ((startIsAhead == increasesWithLayers) ? -1 : +1);

    // Find surrounding layer iterators
    // (Second layer iterator comes immediately after the fit has crossed the target X or Z coordinate
    //  and first layer iterator comes immediately before the second layer iterator).
    firstLayerIter = m_layerFitResultMap.end();
    secondLayerIter = m_layerFitResultMap.end();

    CartesianVector firstLayerPosition(0.f, 0.f, 0.f);
    CartesianVector secondLayerPosition(0.f, 0.f, 0.f);

    for (int iLayer = startLayerIter->first; (iLayer >= minLayer) && (iLayer <= maxLayer); iLayer += increment)
    {
        LayerFitResultMap::const_iterator tempIter = m_layerFitResultMap.find(iLayer);
        if (m_layerFitResultMap.end() == tempIter)
            continue;

        firstLayerIter = secondLayerIter;
        firstLayerPosition = secondLayerPosition;
        secondLayerIter = tempIter;

        this->GetGlobalPosition(secondLayerIter->second.GetL(), secondLayerIter->second.GetFitT(), secondLayerPosition);
        const bool isAhead(useX ? (secondLayerPosition.GetX() > p) : (secondLayerPosition.GetZ() > p));

        if (startIsAhead != isAhead)
            break;

        firstLayerIter = m_layerFitResultMap.end();
    }

    if (m_layerFitResultMap.end() == firstLayerIter || m_layerFitResultMap.end() == secondLayerIter)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitResult::GetLayerInterpolationWeights(const float rL, const LayerFitResultMap::const_iterator &firstLayerIter,
    const LayerFitResultMap::const_iterator &secondLayerIter, float &firstWeight, float &secondWeight) const
{
    if (m_layerFitResultMap.end() == firstLayerIter || m_layerFitResultMap.end() == secondLayerIter)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    
    const float deltaL(rL - firstLayerIter->second.GetL());
    const float deltaLLayers(secondLayerIter->second.GetL() - firstLayerIter->second.GetL());

    if (std::fabs(deltaLLayers) > std::numeric_limits<float>::epsilon())
    {
        firstWeight = 1.f - deltaL / deltaLLayers;
        secondWeight = deltaL / deltaLLayers;
    }
    else
    {
        firstWeight = 0.5f;
        secondWeight = 0.5f;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitResult::GetLayerInterpolationWeights(const float p, const bool useX, const LayerFitResultMap::const_iterator &firstLayerIter,
        const LayerFitResultMap::const_iterator &secondLayerIter, float &firstWeight, float &secondWeight) const
{
    if (m_layerFitResultMap.end() == firstLayerIter || m_layerFitResultMap.end() == secondLayerIter)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    CartesianVector firstLayerPosition(0.f, 0.f, 0.f);
    CartesianVector secondLayerPosition(0.f, 0.f, 0.f);

    this->GetGlobalPosition(firstLayerIter->second.GetL(), firstLayerIter->second.GetFitT(), firstLayerPosition);
    this->GetGlobalPosition(secondLayerIter->second.GetL(), secondLayerIter->second.GetFitT(), secondLayerPosition);

    const float deltaP = useX ? (p - firstLayerPosition.GetX()) : (p - firstLayerPosition.GetZ());
    const float deltaPLayers = useX ? (secondLayerPosition.GetX() - firstLayerPosition.GetX()) : (secondLayerPosition.GetZ() - firstLayerPosition.GetZ());

    if (std::fabs(deltaPLayers) > std::numeric_limits<float>::epsilon())
    {
        firstWeight = 1.f - deltaP / deltaPLayers;
        secondWeight = deltaP / deltaPLayers;
    }
    else
    {
        firstWeight = 0.5f;
        secondWeight = 0.5f;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

TwoDSlidingFitResult::LayerFitResult::LayerFitResult(const double l, const double fitT, const double gradient, const double rms) :
    m_l(l),
    m_fitT(fitT),
    m_gradient(gradient),
    m_rms(rms)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

TwoDSlidingFitResult::LayerFitContribution::LayerFitContribution() :
    m_sumT(0.),
    m_sumL(0.),
    m_sumTT(0.),
    m_sumLT(0.),
    m_sumLL(0.),
    m_nPoints(0)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitResult::LayerFitContribution::AddPoint(const float l, const float t)
{
    const double T = static_cast<double>(t);
    const double L = static_cast<double>(l);

    m_sumT += T;
    m_sumL += L;
    m_sumTT += T * T;
    m_sumLT += L * T;
    m_sumLL += L * L;
    ++m_nPoints;
}

} // namespace lar
