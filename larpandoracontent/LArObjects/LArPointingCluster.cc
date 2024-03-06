/**
 *  @file   larpandoracontent/LArObjects/LArPointingCluster.cc
 *
 *  @brief  Implementation of the lar pointing cluster class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArObjects/LArPointingCluster.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

using namespace pandora;

namespace lar_content
{

LArPointingCluster::LArPointingCluster(const Cluster *const pCluster, const unsigned int fitHalfLayerWindow, const float fitLayerPitch)
{
    // TODO remove default layer fit window and z pitch values
    if (TPC_3D == LArClusterHelper::GetClusterHitType(pCluster))
    {
        const ThreeDSlidingFitResult slidingFitResult(pCluster, fitHalfLayerWindow, fitLayerPitch);
        this->BuildPointingCluster(slidingFitResult);
    }
    else
    {
        const TwoDSlidingFitResult slidingFitResult(pCluster, fitHalfLayerWindow, fitLayerPitch);
        this->BuildPointingCluster(slidingFitResult);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArPointingCluster::LArPointingCluster(const TwoDSlidingFitResult &slidingFitResult)
{
    this->BuildPointingCluster(slidingFitResult);
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArPointingCluster::LArPointingCluster(const ThreeDSlidingFitResult &slidingFitResult)
{
    this->BuildPointingCluster(slidingFitResult);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPointingCluster::BuildPointingCluster(const TwoDSlidingFitResult &slidingFitResult)
{
    const HitType hitType(LArClusterHelper::GetClusterHitType(slidingFitResult.GetCluster()));

    if (!((TPC_VIEW_U == hitType) || (TPC_VIEW_V == hitType) || (TPC_VIEW_W == hitType)))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    if (slidingFitResult.GetMinLayer() >= slidingFitResult.GetMaxLayer())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const bool isInner((slidingFitResult.GetGlobalMinLayerPosition().GetZ() < slidingFitResult.GetGlobalMaxLayerPosition().GetZ()));

    m_pCluster = slidingFitResult.GetCluster();

    const Vertex minVertex(m_pCluster, slidingFitResult.GetGlobalMinLayerPosition(), slidingFitResult.GetGlobalMinLayerDirection(),
        slidingFitResult.GetMinLayerRms(), isInner);
    const Vertex maxVertex(m_pCluster, slidingFitResult.GetGlobalMaxLayerPosition(), slidingFitResult.GetGlobalMaxLayerDirection() * -1.f,
        slidingFitResult.GetMaxLayerRms(), !isInner);

    m_innerVertex = (isInner ? minVertex : maxVertex);
    m_outerVertex = (isInner ? maxVertex : minVertex);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPointingCluster::BuildPointingCluster(const ThreeDSlidingFitResult &slidingFitResult)
{
    if (TPC_3D != LArClusterHelper::GetClusterHitType(slidingFitResult.GetCluster()))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    if (slidingFitResult.GetMinLayer() >= slidingFitResult.GetMaxLayer())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const bool isInner((slidingFitResult.GetGlobalMinLayerPosition().GetZ() < slidingFitResult.GetGlobalMaxLayerPosition().GetZ()) &&
                       (slidingFitResult.GetMinLayer() < slidingFitResult.GetMaxLayer()));

    m_pCluster = slidingFitResult.GetCluster();

    const Vertex minVertex(m_pCluster, slidingFitResult.GetGlobalMinLayerPosition(), slidingFitResult.GetGlobalMinLayerDirection(),
        slidingFitResult.GetMinLayerRms(), isInner);
    const Vertex maxVertex(m_pCluster, slidingFitResult.GetGlobalMaxLayerPosition(), slidingFitResult.GetGlobalMaxLayerDirection() * -1.f,
        slidingFitResult.GetMaxLayerRms(), !isInner);

    m_innerVertex = (isInner ? minVertex : maxVertex);
    m_outerVertex = (isInner ? maxVertex : minVertex);
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArPointingCluster::Vertex::Vertex() :
    m_pCluster(NULL), m_position(0.f, 0.f, 0.f), m_direction(0.f, 0.f, 0.f), m_rms(std::numeric_limits<float>::max()), m_isInner(false), m_isInitialized(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArPointingCluster::Vertex::Vertex(
    const Cluster *const pCluster, const CartesianVector &position, const CartesianVector &direction, const float rms, const bool isInner) :
    m_pCluster(pCluster), m_position(position), m_direction(direction), m_rms(rms), m_isInner(isInner), m_isInitialized(true)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArPointingCluster::Vertex::Vertex(const Vertex &rhs) :
    m_pCluster(rhs.m_pCluster), m_position(rhs.m_position), m_direction(rhs.m_direction), m_rms(rhs.m_rms), m_isInner(rhs.m_isInner), m_isInitialized(rhs.m_isInitialized)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArPointingCluster::Vertex::~Vertex()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArPointingCluster::Vertex &LArPointingCluster::Vertex::operator=(const Vertex &rhs)
{
    m_pCluster = rhs.m_pCluster;
    m_position = rhs.m_position;
    m_direction = rhs.m_direction;
    m_rms = rhs.m_rms;
    m_isInner = rhs.m_isInner;
    m_isInitialized = rhs.m_isInitialized;

    return *this;
}

} // namespace lar_content
