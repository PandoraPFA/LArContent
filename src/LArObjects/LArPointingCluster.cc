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

using namespace pandora;

namespace lar
{

LArPointingCluster::LArPointingCluster::Vertex::Vertex(Cluster *const pCluster, const bool useInnerVertex, const unsigned int nLayersToSkip) :
    m_pCluster(pCluster),
    m_position(0.f, 0.f, 0.f),
    m_direction(0.f, 0.f, 0.f),
    m_rms(std::numeric_limits<float>::max()),
    m_isInner(useInnerVertex),
    m_nSkippedLayers(nLayersToSkip)
{
    static const unsigned int m_numberOfLayersToFit = 30; // TODO, read from settings, probably via LArClusterHelper
    const unsigned int innerLayer(pCluster->GetInnerPseudoLayer());
    const unsigned int outerLayer(pCluster->GetOuterPseudoLayer());

    if (useInnerVertex)
    {
        const unsigned int innerFitLayer(std::min(innerLayer + nLayersToSkip, outerLayer));
        const unsigned int outerFitLayer(std::min(innerLayer + nLayersToSkip + m_numberOfLayersToFit, outerLayer));

        if (outerFitLayer <= innerFitLayer)
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        ClusterHelper::ClusterFitResult clusterFitResult;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, ClusterHelper::FitLayers(pCluster, innerFitLayer, outerFitLayer, clusterFitResult)); 

        const CartesianVector &fitIntercept(clusterFitResult.GetIntercept());
        const CartesianVector fitDirection(clusterFitResult.GetDirection());
        const CartesianVector innerCentroid(pCluster->GetCentroid(pCluster->GetInnerPseudoLayer()));

        m_rms = clusterFitResult.GetRms();
        m_position = fitIntercept + fitDirection * (fitDirection.GetDotProduct(innerCentroid - fitIntercept));
        m_direction = fitDirection;
    }
    else
    {
        const int innerFitLayer(std::max(static_cast<int>(outerLayer) - static_cast<int>(nLayersToSkip) - static_cast<int>(m_numberOfLayersToFit), static_cast<int>(innerLayer)));
        const int outerFitLayer(std::max(static_cast<int>(outerLayer) - static_cast<int>(nLayersToSkip), static_cast<int>(innerLayer)));

        if (outerFitLayer <= innerFitLayer)
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        ClusterHelper::ClusterFitResult clusterFitResult;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, ClusterHelper::FitLayers(pCluster, innerFitLayer, outerFitLayer, clusterFitResult)); 

        const CartesianVector &fitIntercept(clusterFitResult.GetIntercept());
        const CartesianVector fitDirection(clusterFitResult.GetDirection() * -1.f);
        const CartesianVector outerCentroid(pCluster->GetCentroid(pCluster->GetOuterPseudoLayer()));

        m_rms = clusterFitResult.GetRms();
        m_position = fitIntercept + fitDirection * (fitDirection.GetDotProduct(outerCentroid - fitIntercept));
        m_direction = fitDirection;
    }
}

} // namespace lar
