/**
 *  @file   LArContent/src/Objects/LArPointingCluster.cc
 * 
 *  @brief  Implementation of the lar pointing cluster class.
 * 
 *  $Log: $
 */

#include "Helpers/ClusterHelper.h"

#include "Objects/CartesianVector.h"
#include "Objects/Cluster.h"

#include "Objects/LArPointingCluster.h"

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
    this->GetVertex( 0, m_position, m_direction, m_rms );

  /*
    // KEEP OLD IMPLEMENTATION FOR NOW
    const unsigned int m_numberOfLayersToFit = 30; // TODO, read from settings, probably via LArClusterHelper

    if (m_isInner)
    {
        ClusterHelper::ClusterFitResult clusterFitResult;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, ClusterHelper::FitStart(pCluster, m_numberOfLayersToFit, clusterFitResult)); 

        const CartesianVector &fitIntercept(clusterFitResult.GetIntercept());
        const CartesianVector fitDirection(clusterFitResult.GetDirection());
        const CartesianVector innerCentroid(pCluster->GetCentroid(pCluster->GetInnerPseudoLayer()));

        m_rms = clusterFitResult.GetRms();
        m_position = fitIntercept + fitDirection * (fitDirection.GetDotProduct(innerCentroid - fitIntercept));
        m_direction = fitDirection;
    }
    else
    {
        ClusterHelper::ClusterFitResult clusterFitResult;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, ClusterHelper::FitEnd(pCluster, m_numberOfLayersToFit, clusterFitResult)); 

        const CartesianVector &fitIntercept(clusterFitResult.GetIntercept());
        const CartesianVector fitDirection(clusterFitResult.GetDirection() * -1.f);
        const CartesianVector outerCentroid(pCluster->GetCentroid(pCluster->GetOuterPseudoLayer()));

        m_rms = clusterFitResult.GetRms();
        m_position = fitIntercept + fitDirection * (fitDirection.GetDotProduct(outerCentroid - fitIntercept));
        m_direction = fitDirection;
    }
  */
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPointingCluster::Vertex::GetVertex( const unsigned int nLayersToSkip, pandora::CartesianVector& position, pandora::CartesianVector& direction, float& rms ) const
{
    const unsigned int m_numberOfLayersToFit = 30;    // TODO, read from settings, probably via LArClusterHelper
    const unsigned int m_minNumberOfLayersToFit = 20; // TODO, read from settings, probably via LArClusterHelper

    const unsigned int innerLayer = m_pCluster->GetInnerPseudoLayer();
    const unsigned int outerLayer = m_pCluster->GetOuterPseudoLayer();

    if ( m_isInner )
    {
        // Calculate inner and outer layers 
        unsigned int innerFitLayer = innerLayer;
        unsigned int outerFitLayer = outerLayer;

        if ( innerLayer + nLayersToSkip < outerLayer )
	    innerFitLayer = innerLayer + nLayersToSkip;
        else innerFitLayer = outerLayer;

        if ( innerFitLayer + (m_numberOfLayersToFit-1) < outerLayer )
	    outerFitLayer = innerFitLayer + (m_numberOfLayersToFit-1);
        else outerFitLayer = outerLayer;

        if ( innerFitLayer + (m_minNumberOfLayersToFit-1) > outerFitLayer )
	{
	    if ( innerLayer + (m_minNumberOfLayersToFit-1) < outerFitLayer )
	        innerFitLayer = outerFitLayer - (m_minNumberOfLayersToFit-1);
            else innerFitLayer = innerLayer;
	}

        // Apply a linear fit to these layers
        ClusterHelper::ClusterFitResult clusterFitResult;

        if ( STATUS_CODE_SUCCESS == ClusterHelper::FitLayers(m_pCluster, innerFitLayer, outerFitLayer, clusterFitResult) )
	{
            const CartesianVector &fitIntercept(clusterFitResult.GetIntercept());
            const CartesianVector fitDirection(clusterFitResult.GetDirection());
            const CartesianVector innerCentroid(m_pCluster->GetCentroid(innerLayer));

            rms = clusterFitResult.GetRms();
            position = fitIntercept + fitDirection * (fitDirection.GetDotProduct(innerCentroid - fitIntercept));
            direction = fitDirection;
	}
        else
	{
	    // Try to use the inner and outer centroids
	    rms = std::numeric_limits<float>::max();

            const CartesianVector innerCentroid(m_pCluster->GetCentroid(innerLayer));
            const CartesianVector outerCentroid(m_pCluster->GetCentroid(outerLayer));

            if ( (outerCentroid - innerCentroid).GetMagnitudeSquared() > std::numeric_limits<float>::epsilon() )
	    {
	        position = innerCentroid;
                direction = (outerCentroid - innerCentroid).GetUnitVector();
	    }
            else 
            {
                // Single-layer clusters aren't pointing clusters!
                throw pandora::StatusCodeException(STATUS_CODE_NOT_ALLOWED);
	    }
	}
    }

    else
    {   
        // Calculate inner and outer layers
        unsigned int innerFitLayer = innerLayer;
        unsigned int outerFitLayer = outerLayer;

        if ( outerLayer > innerLayer + nLayersToSkip ) 
            outerFitLayer = outerLayer - nLayersToSkip; 
        else outerFitLayer = innerLayer;

        if ( innerLayer + (m_numberOfLayersToFit-1) < outerFitLayer )
	    innerFitLayer = outerFitLayer - (m_numberOfLayersToFit-1);
        else innerFitLayer = innerLayer;

        if ( innerFitLayer + (m_minNumberOfLayersToFit-1) > outerFitLayer )
	{
	    if ( innerFitLayer + (m_minNumberOfLayersToFit-1) < outerLayer )
	        outerFitLayer = innerFitLayer + (m_minNumberOfLayersToFit-1);
	    else outerFitLayer = outerLayer;
	}

        // Apply a linear fit to these layers
        ClusterHelper::ClusterFitResult clusterFitResult;

        if ( STATUS_CODE_SUCCESS == ClusterHelper::FitLayers(m_pCluster, innerFitLayer, outerFitLayer, clusterFitResult) )
	{
            const CartesianVector &fitIntercept(clusterFitResult.GetIntercept());
            const CartesianVector fitDirection(clusterFitResult.GetDirection() * -1.f);
            const CartesianVector outerCentroid(m_pCluster->GetCentroid(outerLayer));

            rms = clusterFitResult.GetRms();
            position = fitIntercept + fitDirection * (fitDirection.GetDotProduct(outerCentroid - fitIntercept));
            direction = fitDirection;
	}
        else
	{
            // Try to use the inner and outer centroids
            rms = std::numeric_limits<float>::max();

            const CartesianVector innerCentroid(m_pCluster->GetCentroid(innerLayer));
            const CartesianVector outerCentroid(m_pCluster->GetCentroid(outerLayer));

            if ( (outerCentroid - innerCentroid).GetMagnitudeSquared() > std::numeric_limits<float>::epsilon() )
	    {
                position = outerCentroid;
                direction = (innerCentroid - outerCentroid).GetUnitVector();
	    }
            else
	    { 
                // Single-layer clusters aren't pointing clusters!
                throw pandora::StatusCodeException(STATUS_CODE_NOT_ALLOWED);
	    }
	}
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArPointingCluster::GetLengthSquared() const
{
    return (m_outerVertex.GetPosition() - m_innerVertex.GetPosition()).GetMagnitudeSquared();
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArPointingCluster::GetLength() const
{
    return std::sqrt( GetLengthSquared() );
}

} // namespace lar
