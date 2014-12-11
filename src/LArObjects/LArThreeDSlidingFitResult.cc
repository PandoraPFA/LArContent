/**
 *  @file   LArContent/src/LArObjects/LArThreeDSlidingFitResult.cc
 *
 *  @brief  Implementation of the lar three dimensional sliding fit result class.
 *
 *  $Log: $
 */

#include "Objects/Cluster.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArObjects/LArThreeDSlidingFitResult.h"

#include <algorithm>
#include <cmath>
#include <limits>

using namespace pandora;

namespace lar_content
{

ThreeDSlidingFitResult::ThreeDSlidingFitResult(const Cluster *const pCluster, const unsigned int layerWindow, const float layerPitch) :
    m_pCluster(pCluster),
    m_primaryAxis(ThreeDSlidingFitResult::GetPrimaryAxis(m_pCluster)),
    m_axisIntercept(m_primaryAxis.GetPosition()),
    m_axisDirection(m_primaryAxis.GetMomentum()),
    m_firstOrthoDirection(ThreeDSlidingFitResult::GetSeedDirection(m_axisDirection).GetCrossProduct(m_axisDirection).GetUnitVector()),
    m_secondOrthoDirection(m_axisDirection.GetCrossProduct(m_firstOrthoDirection).GetUnitVector()),
    m_firstFitResult(TwoDSlidingFitResult(pCluster, layerWindow, layerPitch, m_axisIntercept, m_axisDirection, m_firstOrthoDirection)),
    m_secondFitResult(TwoDSlidingFitResult(pCluster, layerWindow, layerPitch, m_axisIntercept, m_axisDirection, m_secondOrthoDirection))
{

}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackState ThreeDSlidingFitResult::GetPrimaryAxis(const Cluster *pCluster)
{
    CartesianVector innerCoordinate(0.f, 0.f, 0.f), outerCoordinate(0.f, 0.f, 0.f);
    LArClusterHelper::GetExtremalCoordinates(pCluster, innerCoordinate, outerCoordinate);
    return TrackState(innerCoordinate, (outerCoordinate - innerCoordinate).GetUnitVector());
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

} // namespace lar_content
