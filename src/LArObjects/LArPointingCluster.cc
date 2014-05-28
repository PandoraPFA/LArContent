/**
 *  @file   LArContent/src/LArObjects/LArPointingCluster.cc
 *
 *  @brief  Implementation of the lar pointing cluster class.
 *
 *  $Log: $
 */

#include "Helpers/ClusterHelper.h"

#include "Objects/CartesianVector.h"
#include "Objects/Cluster.h"

#include "LArObjects/LArPointingCluster.h"
#include "LArHelpers/LArClusterHelper.h"

using namespace pandora;

namespace lar
{

LArPointingCluster::LArPointingCluster(pandora::Cluster *const pCluster) :
    m_pCluster(pCluster)
{
    TwoDSlidingFitResult slidingFitResult;
    LArClusterHelper::LArTwoDSlidingFit(pCluster, 10, slidingFitResult);

    const float minLayerZ(slidingFitResult.GetGlobalMinLayerPosition().GetZ());
    const float maxLayerZ(slidingFitResult.GetGlobalMaxLayerPosition().GetZ());

    const int minLayer(slidingFitResult.GetMinLayer());
    const int maxLayer(slidingFitResult.GetMaxLayer());

    if (minLayer >= maxLayer)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const int innerLayer((minLayerZ < maxLayerZ) ? minLayer : maxLayer);
    const int outerLayer((minLayerZ < maxLayerZ) ? maxLayer : minLayer);

    if (innerLayer == minLayer)
    {
        m_innerVertex = Vertex(pCluster, slidingFitResult.GetGlobalMinLayerPosition(), slidingFitResult.GetGlobalMinLayerDirection(),
            slidingFitResult.GetMinLayerRms(), true);
    }
    else
    {
        m_innerVertex = Vertex(pCluster, slidingFitResult.GetGlobalMaxLayerPosition(), slidingFitResult.GetGlobalMaxLayerDirection() * -1.f,
            slidingFitResult.GetMaxLayerRms(), true);
    }

    if (outerLayer == minLayer)
    {
        m_outerVertex = Vertex(pCluster, slidingFitResult.GetGlobalMinLayerPosition(), slidingFitResult.GetGlobalMinLayerDirection(),
            slidingFitResult.GetMinLayerRms(), false);
    }
    else
    {
        m_outerVertex = Vertex(pCluster, slidingFitResult.GetGlobalMaxLayerPosition(), slidingFitResult.GetGlobalMaxLayerDirection() * -1.f,
            slidingFitResult.GetMaxLayerRms(), false);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LArPointingCluster::Vertex::Vertex() :
    m_pCluster(NULL),
    m_position(0.f, 0.f, 0.f),
    m_direction(0.f, 0.f, 0.f),
    m_rms(std::numeric_limits<float>::max()),
    m_isInner(false),
    m_isInitialized(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArPointingCluster::Vertex::Vertex(Cluster *const pCluster, const CartesianVector &position, const CartesianVector &direction,
        const float rms, const bool isInner) :
    m_pCluster(pCluster),
    m_position(position),
    m_direction(direction),
    m_rms(rms),
    m_isInner(isInner),
    m_isInitialized(true)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArPointingCluster::Vertex::Vertex(const Vertex &rhs) :
    m_pCluster(rhs.m_pCluster),
    m_position(rhs.m_position),
    m_direction(rhs.m_direction),
    m_rms(rhs.m_rms),
    m_isInner(rhs.m_isInner),
    m_isInitialized(rhs.m_isInitialized)
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

} // namespace lar
