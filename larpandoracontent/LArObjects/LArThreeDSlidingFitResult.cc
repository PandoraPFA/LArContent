/**
 *  @file   larpandoracontent/LArObjects/LArThreeDSlidingFitResult.cc
 *
 *  @brief  Implementation of the lar three dimensional sliding fit result class.
 *
 *  $Log: $
 */

#include "Helpers/ClusterFitHelper.h"

#include "Objects/Cluster.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPCAHelper.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include <algorithm>
#include <cmath>
#include <limits>

using namespace pandora;

namespace lar_content
{

template <typename T>
ThreeDSlidingFitResult::ThreeDSlidingFitResult(const T *const pT, const unsigned int layerWindow, const float layerPitch) :
    m_primaryAxis(ThreeDSlidingFitResult::GetPrimaryAxis(pT)),
    m_axisIntercept(m_primaryAxis.GetPosition()),
    m_axisDirection(m_primaryAxis.GetMomentum()),
    m_firstOrthoDirection(ThreeDSlidingFitResult::GetSeedDirection(m_axisDirection).GetCrossProduct(m_axisDirection).GetUnitVector()),
    m_secondOrthoDirection(m_axisDirection.GetCrossProduct(m_firstOrthoDirection).GetUnitVector()),
    m_firstFitResult(TwoDSlidingFitResult(pT, layerWindow, layerPitch, m_axisIntercept, m_axisDirection, m_firstOrthoDirection)),
    m_secondFitResult(TwoDSlidingFitResult(pT, layerWindow, layerPitch, m_axisIntercept, m_axisDirection, m_secondOrthoDirection)),
    m_minLayer(std::max(m_firstFitResult.GetMinLayer(), m_secondFitResult.GetMinLayer())),
    m_maxLayer(std::min(m_firstFitResult.GetMaxLayer(), m_secondFitResult.GetMaxLayer())),
    m_minLayerPosition(0.f, 0.f, 0.f),
    m_maxLayerPosition(0.f, 0.f, 0.f),
    m_minLayerDirection(0.f, 0.f, 0.f),
    m_maxLayerDirection(0.f, 0.f, 0.f)
{
    if (m_minLayer > m_maxLayer)
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    const float minL(m_firstFitResult.GetL(m_minLayer));
    const float maxL(m_firstFitResult.GetL(m_maxLayer));

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetGlobalFitPosition(minL, m_minLayerPosition));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetGlobalFitPosition(maxL, m_maxLayerPosition));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetGlobalFitDirection(minL, m_minLayerDirection));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetGlobalFitDirection(maxL, m_maxLayerDirection));
}

//------------------------------------------------------------------------------------------------------------------------------------------

const pandora::Cluster *ThreeDSlidingFitResult::GetCluster() const
{
    return m_firstFitResult.GetCluster();
}

//------------------------------------------------------------------------------------------------------------------------------------------

int ThreeDSlidingFitResult::GetMinLayer() const
{
    const int firstLayer(m_firstFitResult.GetMinLayer());
    const int secondLayer(m_secondFitResult.GetMinLayer());

    return std::min(firstLayer, secondLayer);
}

//------------------------------------------------------------------------------------------------------------------------------------------

int ThreeDSlidingFitResult::GetMaxLayer() const
{
    const int firstLayer(m_firstFitResult.GetMaxLayer());
    const int secondLayer(m_secondFitResult.GetMaxLayer());

    return std::max(firstLayer, secondLayer);
}

//------------------------------------------------------------------------------------------------------------------------------------------
 
float ThreeDSlidingFitResult::GetMinLayerRms() const
{
    const float firstRms(m_firstFitResult.GetMinLayerRms());
    const float secondRms(m_secondFitResult.GetMinLayerRms());

    return std::sqrt(firstRms * firstRms + secondRms * secondRms);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ThreeDSlidingFitResult::GetMaxLayerRms() const
{
    const float firstRms(m_firstFitResult.GetMaxLayerRms());
    const float secondRms(m_secondFitResult.GetMaxLayerRms());

    return std::sqrt(firstRms * firstRms + secondRms * secondRms);
}

//------------------------------------------------------------------------------------------------------------------------------------------
   
float ThreeDSlidingFitResult::GetFitRms(const float rL) const
{
    const float firstRms(m_firstFitResult.GetFitRms(rL));
    const float secondRms(m_secondFitResult.GetFitRms(rL));

    return std::sqrt(firstRms * firstRms + secondRms * secondRms);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ThreeDSlidingFitResult::GetLongitudinalDisplacement(const CartesianVector &position) const
{
    return m_axisDirection.GetDotProduct(position - m_axisIntercept);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDSlidingFitResult::GetGlobalFitPosition(const float rL, CartesianVector &position) const
{
    // Check that input coordinates are between first and last layers
    const int layer1(m_firstFitResult.GetLayer(rL));
    const int layer2(m_secondFitResult.GetLayer(rL));

    if (std::min(layer1, layer2) < m_minLayer || std::max(layer1, layer2) > m_maxLayer)
        return STATUS_CODE_INVALID_PARAMETER;

    // Get local positions from each sliding fit (TODO: Make this more efficient)
    CartesianVector firstPosition(0.f, 0.f, 0.f), secondPosition(0.f, 0.f, 0.f);
    const StatusCode statusCode1(m_firstFitResult.GetGlobalFitPosition(rL, firstPosition));

    if (STATUS_CODE_SUCCESS != statusCode1)
        return statusCode1;

    const StatusCode statusCode2(m_secondFitResult.GetGlobalFitPosition(rL, secondPosition));

    if (STATUS_CODE_SUCCESS != statusCode2)
        return statusCode2;

    float rL1(0.f), rT1(0.f), rL2(0.f), rT2(0.f);
    m_firstFitResult.GetLocalPosition(firstPosition, rL1, rT1);
    m_secondFitResult.GetLocalPosition(secondPosition, rL2, rT2);

    // Combine local positions to give an overall global direction
    this->GetGlobalPosition(rL, rT1, rT2, position);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDSlidingFitResult::GetGlobalFitDirection(const float rL, CartesianVector &direction) const
{
    // Check that input coordinates are between first and last layers
    const int layer1(m_firstFitResult.GetLayer(rL));
    const int layer2(m_secondFitResult.GetLayer(rL));

    if (std::min(layer1, layer2) < m_minLayer || std::max(layer1, layer2) > m_maxLayer)
        return STATUS_CODE_INVALID_PARAMETER;

    // Get local directions from each sliding fit (TODO: Make this more efficient)
    CartesianVector firstDirection(0.f, 0.f, 0.f), secondDirection(0.f, 0.f, 0.f);
    const StatusCode statusCode1(m_firstFitResult.GetGlobalFitDirection(rL, firstDirection));

    if (STATUS_CODE_SUCCESS != statusCode1)
        return statusCode1;

    const StatusCode statusCode2(m_secondFitResult.GetGlobalFitDirection(rL, secondDirection));

    if (STATUS_CODE_SUCCESS != statusCode2)
        return statusCode2;

    float dTdL1(0.f), dTdL2(0.f);
    m_firstFitResult.GetLocalDirection(firstDirection, dTdL1);
    m_secondFitResult.GetLocalDirection(secondDirection, dTdL2);

    // Combine local directions to give an overall global direction
    this->GetGlobalDirection(dTdL1, dTdL2, direction);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDSlidingFitResult::GetGlobalPosition(const float rL, const float rT1, const float rT2, CartesianVector &position) const
{
    position = m_axisIntercept + m_axisDirection * rL + m_firstOrthoDirection * rT1 + m_secondOrthoDirection * rT2;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDSlidingFitResult::GetGlobalDirection(const float dTdL1, const float dTdL2, CartesianVector &direction) const
{
    const float pL(1.f / std::sqrt(1.f + dTdL1 * dTdL1 + dTdL2 * dTdL2));
    const float pT1(dTdL1 / std::sqrt(1.f + dTdL1 * dTdL1 + dTdL2 * dTdL2));
    const float pT2(dTdL2 / std::sqrt(1.f + dTdL1 * dTdL1 + dTdL2 * dTdL2));

    CartesianVector globalCoordinates(0.f, 0.f, 0.f);
    this->GetGlobalPosition(pL, pT1, pT2, globalCoordinates);
    direction = (globalCoordinates - m_axisIntercept).GetUnitVector();
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackState ThreeDSlidingFitResult::GetPrimaryAxis(const Cluster *const pCluster)
{
    CartesianPointVector pointVector;
    LArClusterHelper::GetCoordinateVector(pCluster, pointVector);
    return ThreeDSlidingFitResult::GetPrimaryAxis(&pointVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackState ThreeDSlidingFitResult::GetPrimaryAxis(const CartesianPointVector *const pPointVector)
{
    if (pPointVector->size() < 2)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    CartesianVector centroid(0.f, 0.f, 0.f);
    LArPCAHelper::EigenVectors eigenVecs;
    LArPCAHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
    LArPCAHelper::RunPCA(*pPointVector, centroid, eigenValues, eigenVecs);

    float minProjection(std::numeric_limits<float>::max());
    const CartesianVector &fitDirection(eigenVecs.at(0));

    for (const CartesianVector &coordinate : *pPointVector)
        minProjection = std::min(minProjection, fitDirection.GetDotProduct(coordinate - centroid));

    return TrackState(centroid + (fitDirection * minProjection), fitDirection);
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector ThreeDSlidingFitResult::GetSeedDirection(const CartesianVector &axisDirection)
{
    const float px(std::fabs(axisDirection.GetX()));
    const float py(std::fabs(axisDirection.GetY()));
    const float pz(std::fabs(axisDirection.GetZ()));

    if (px < std::min(py,pz) + std::numeric_limits<float>::epsilon())
    {
        return CartesianVector(1.f, 0.f, 0.f);
    }

    if (py < std::min(pz,px) + std::numeric_limits<float>::epsilon())
    {
        return CartesianVector(0.f, 1.f, 0.f);
    }

    if (pz < std::min(px,py) + std::numeric_limits<float>::epsilon())
    {
        return CartesianVector(0.f, 0.f, 1.f);
    }

    throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

template ThreeDSlidingFitResult::ThreeDSlidingFitResult(const pandora::Cluster *const, const unsigned int, const float);
template ThreeDSlidingFitResult::ThreeDSlidingFitResult(const pandora::CartesianPointVector *const, const unsigned int, const float);

} // namespace lar_content
