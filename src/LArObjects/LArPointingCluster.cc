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

LArPointingCluster::LArPointingCluster::Vertex::Vertex(Cluster *const pCluster, const bool useInnerVertex) :
    m_pCluster(pCluster),
    m_position(0.f, 0.f, 0.f),
    m_direction(0.f, 0.f, 0.f),
    m_rms(std::numeric_limits<float>::max()),
    m_isInner(useInnerVertex)
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
    const int fitLayer(useInnerVertex ? innerLayer : outerLayer);

    if (fitLayer == minLayer)
    {
        m_position  = slidingFitResult.GetGlobalMinLayerPosition();
        m_direction = slidingFitResult.GetGlobalMinLayerDirection();
        m_rms       = slidingFitResult.GetGlobalMinLayerRms();
    }
    else
    {
        m_position  = slidingFitResult.GetGlobalMaxLayerPosition();
        m_direction = slidingFitResult.GetGlobalMaxLayerDirection() * -1.f;
        m_rms       = slidingFitResult.GetGlobalMaxLayerRms();
    }
}

} // namespace lar
