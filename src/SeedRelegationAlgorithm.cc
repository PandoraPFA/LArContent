/**
 *  @file   SeedRelegationAlgorithm.cc
 * 
 *  @brief  Implementation of the seed relegation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "SeedRelegationAlgorithm.h"

#include "LArVertexHelper.h"
#include "LArParticleId.h"

using namespace pandora;

namespace lar
{

StatusCode SeedRelegationAlgorithm::Run()
{
    // Read in Cluster Lists
    const ClusterList *pSeedClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetClusterList(*this, m_seedClusterListName, pSeedClusterList));

    const ClusterList *pNonSeedClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetClusterList(*this, m_nonSeedClusterListName, pNonSeedClusterList));

    ClusterList clustersToDemote;

    for( ClusterList::const_iterator iterI = pSeedClusterList->begin(), iterEndI = pSeedClusterList->end(); iterI != iterEndI; ++iterI )
      {
        Cluster* seedClusterI = *iterI; 

        // Simple checks on length and proximity to vertex
        if( IsSeedClusterLongEnough(seedClusterI)==true ) continue;

        if( IsSeedClusterConnectedToVertex(seedClusterI)==true ) continue;
       
        // Check against the other Seed Clusters
        bool relegateThisCluster(false);

        for ( ClusterList::const_iterator iterJ = pSeedClusterList->begin(), iterEndJ = pSeedClusterList->end(); !relegateThisCluster && iterJ != iterEndJ; ++iterJ )
        {
            Cluster* seedClusterJ = *iterJ; 

            if( seedClusterI==seedClusterJ ) continue;

            if( IsSeedClusterConnectedToLargerCluster(seedClusterI,seedClusterJ)==true )
	        relegateThisCluster = true;
	}
     
        if ( relegateThisCluster )  clustersToDemote.insert(seedClusterI);  
      }

    // Move demoted clusters from Seed List to Non-Seed List
    if( clustersToDemote.empty()==false )
      {

// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&clustersToDemote, "Demoted", RED); // Relegation colours are red and white...
// PandoraMonitoringApi::ViewEvent();

         PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveClusterList(*this, 
             m_seedClusterListName, m_nonSeedClusterListName, clustersToDemote));
      }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SeedRelegationAlgorithm::IsSeedClusterLongEnough( const Cluster* const pCluster ) const
{
    const CartesianVector innerCentroid  = pCluster->GetCentroid( pCluster->GetInnerPseudoLayer() );
    const CartesianVector outerCentroid  = pCluster->GetCentroid( pCluster->GetOuterPseudoLayer() );
    const float lengthSquared  = (outerCentroid-innerCentroid).GetMagnitudeSquared();

    if (lengthSquared > m_minClusterLengthSquared)  return true;  else return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SeedRelegationAlgorithm::IsSeedClusterConnectedToVertex( const Cluster* const pCluster ) const
{
    if( LArVertexHelper::DoesCurrentVertexExist()==false ) return false;

    const float clusterImpactParameter( LArVertexHelper::GetImpactParameterToCurrentVertex(pCluster) );
    const float clusterRadialDistanceSquared( LArVertexHelper::GetDistanceSquaredToCurrentVertex(pCluster) );

    if (clusterRadialDistanceSquared < m_vertexInnerRadialDistanceSquared) 
          return true;

    if (clusterRadialDistanceSquared < m_vertexOuterRadialDistanceSquared && clusterImpactParameter<m_vertexImpactParameter) 
          return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SeedRelegationAlgorithm::IsSeedClusterConnectedToLargerCluster( const Cluster* const pClusterI, const Cluster* const pClusterJ )
{
    // Bail out if Cluster I is sizeable compared with Cluster J
    if ( static_cast<float>(pClusterI->GetNCaloHits()) > 0.15*static_cast<float>(pClusterJ->GetNCaloHits()) )
        return false;

    // Calculate inner and outer positions for each cluster 
    const unsigned int innerLayerI = pClusterI->GetInnerPseudoLayer();
    const unsigned int outerLayerI = pClusterI->GetOuterPseudoLayer();
 
    const unsigned int innerLayerJ = pClusterJ->GetInnerPseudoLayer();
    const unsigned int outerLayerJ = pClusterJ->GetOuterPseudoLayer();

    const CartesianVector innerCentroidI = pClusterI->GetCentroid( innerLayerI );
    const CartesianVector outerCentroidI = pClusterI->GetCentroid( outerLayerI );

    const CartesianVector innerCentroidJ = pClusterJ->GetCentroid( innerLayerJ );
    const CartesianVector outerCentroidJ = pClusterJ->GetCentroid( outerLayerJ );

    const float clusterLengthSquaredI = (outerCentroidI-innerCentroidI).GetMagnitudeSquared();
    const float clusterLengthSquaredJ = (outerCentroidJ-innerCentroidJ).GetMagnitudeSquared();
    const float clusterHalfLengthSquaredJ = 0.25*clusterLengthSquaredJ;

    // Apply linear fit to each cluster
    ClusterHelper::ClusterFitResult clusterFitI( pClusterI->GetFitToAllHitsResult() );
    ClusterHelper::ClusterFitResult clusterFitJ( pClusterJ->GetFitToAllHitsResult() );

    const CartesianVector clusterDirectionI = clusterFitI.GetDirection();
    const CartesianVector clusterDirectionJ = clusterFitJ.GetDirection();
   
    // Quick check: transverse separation of Cluster I and Cluster J 
    if ( std::min(clusterDirectionJ.GetCrossProduct(innerCentroidI-innerCentroidJ).GetMagnitudeSquared(),
                  clusterDirectionJ.GetCrossProduct(outerCentroidI-innerCentroidJ).GetMagnitudeSquared()) 
                > m_minClusterSeparationSquared )
        return false;

    // Quick check: longitudinal separation of Cluster I and Cluster J
    if ( std::max((outerCentroidI-innerCentroidJ).GetMagnitudeSquared(),(outerCentroidJ-innerCentroidI).GetMagnitudeSquared())
	  > 3.0*clusterLengthSquaredJ )
        return false;

    // Quick check: relative angle of Cluster I and Cluster J
    if( clusterDirectionI.GetDotProduct(clusterDirectionJ) < m_minClusterDotProduct )
        return false;

    // Look up the direction of each cluster
    const bool isForwardI  = LArVertexHelper::DoesCurrentVertexExist() ? LArVertexHelper::IsForwardInZ(pClusterI) : false;
    const bool isBackwardI = LArVertexHelper::DoesCurrentVertexExist() ? LArVertexHelper::IsBackwardInZ(pClusterI) : false;
 
    const bool checkForwardI  = ( true==isForwardI  || false==isBackwardI );
    const bool checkBackwardI = ( true==isBackwardI || false==isForwardI  );
    
    const bool isForwardJ  = LArVertexHelper::DoesCurrentVertexExist() ? LArVertexHelper::IsForwardInZ(pClusterJ) : false;
    const bool isBackwardJ = LArVertexHelper::DoesCurrentVertexExist() ? LArVertexHelper::IsBackwardInZ(pClusterJ) : false;
  
    const bool checkForwardJ  = ( true==isForwardJ  || false==isBackwardJ );
    const bool checkBackwardJ = ( true==isBackwardJ || false==isForwardJ  );

    // Check separation of inner and outer centroids of Cluster I to Cluster J 
    if ( ( true==checkForwardI  && true==checkForwardJ  && (innerCentroidI-innerCentroidJ).GetMagnitudeSquared() < m_minClusterVertexSeparationSquared )
      || ( true==checkBackwardI && true==checkForwardJ  && (outerCentroidI-innerCentroidJ).GetMagnitudeSquared() < m_minClusterVertexSeparationSquared )
      || ( true==checkForwardI  && true==checkBackwardJ && (innerCentroidI-outerCentroidJ).GetMagnitudeSquared() < m_minClusterVertexSeparationSquared )
      || ( true==checkBackwardI && true==checkBackwardJ && (outerCentroidI-outerCentroidJ).GetMagnitudeSquared() < m_minClusterVertexSeparationSquared ) )
        return false;

    // Now relegate small seed clusters that subtend small angles to large clusters
    if ( true==checkForwardI && true==checkForwardJ
      && (innerCentroidI-innerCentroidJ).GetMagnitudeSquared() > std::max(m_minClusterRadialSeparationSquared,clusterHalfLengthSquaredJ)
      && this->IsSeedClusterAtSmallAngleToLargerCluster(innerCentroidI, clusterDirectionI, innerCentroidJ, clusterDirectionJ) )
        return true;	

    if ( true==checkBackwardI && true==checkForwardJ
       && (outerCentroidI-innerCentroidJ).GetMagnitudeSquared() > std::max(m_minClusterRadialSeparationSquared,clusterHalfLengthSquaredJ)
       && this->IsSeedClusterAtSmallAngleToLargerCluster(outerCentroidI, clusterDirectionI*-1.0, innerCentroidJ, clusterDirectionJ) )
        return true;    
 
    if ( true==checkForwardI && true==checkBackwardJ
      && (innerCentroidI-outerCentroidJ).GetMagnitudeSquared() > std::max(m_minClusterRadialSeparationSquared,clusterHalfLengthSquaredJ) 
      && this->IsSeedClusterAtSmallAngleToLargerCluster(innerCentroidI, clusterDirectionI, outerCentroidJ, clusterDirectionJ*-1.0) )
        return true;

    if ( true==checkBackwardI && true==checkBackwardJ
     && (outerCentroidI-outerCentroidJ).GetMagnitudeSquared() > std::max(m_minClusterRadialSeparationSquared,clusterHalfLengthSquaredJ)
     && this->IsSeedClusterAtSmallAngleToLargerCluster(outerCentroidI, clusterDirectionI*-1.0, outerCentroidJ, clusterDirectionJ*-1.0) )
        return true;
	    
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SeedRelegationAlgorithm::IsSeedClusterAtSmallAngleToLargerCluster( const pandora::CartesianVector positionI, const pandora::CartesianVector directionI, const pandora::CartesianVector positionJ, const pandora::CartesianVector directionJ )
{
    const CartesianVector p0(directionJ);
    const CartesianVector p1((positionI-positionJ).GetUnitVector()); 
    const CartesianVector p2(directionI);

    if ( (180.0/M_PI) * p0.GetOpeningAngle(p1) < 7.5 && (180.0/M_PI) * p1.GetOpeningAngle(p2) < 15.0 )   // TODO: move to configuration
        return true;

    else return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SeedRelegationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{ 
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "SeedClusterListName", m_seedClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NonSeedClusterListName", m_nonSeedClusterListName));

    m_vertexInnerRadialDistance = 25.0; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "VertexInnerRadialDistance", m_vertexInnerRadialDistance));
    m_vertexInnerRadialDistanceSquared = m_vertexInnerRadialDistance*m_vertexInnerRadialDistance;

    m_vertexOuterRadialDistance = 75.0; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "VertexOuterRadialDistance", m_vertexOuterRadialDistance));
    m_vertexOuterRadialDistanceSquared = m_vertexOuterRadialDistance*m_vertexOuterRadialDistance;

    m_vertexImpactParameter = 5.0; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "VertexImpactParameter", m_vertexImpactParameter));

    m_minClusterLength = 50.0; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "GoodClusterLength", m_minClusterLength));
    m_minClusterLengthSquared = m_minClusterLength*m_minClusterLength;

    m_minClusterSeparation = 50.0; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "GoodClusterSeparation", m_minClusterSeparation));
    m_minClusterSeparationSquared = m_minClusterSeparation*m_minClusterSeparation;

    m_minClusterRadialSeparation = 50.0; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "GoodClusterRadialSeparation", m_minClusterRadialSeparation));
    m_minClusterRadialSeparationSquared = m_minClusterRadialSeparation*m_minClusterRadialSeparation;

    m_minClusterVertexSeparation = 5.0; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "GoodClusterVertexSeparation", m_minClusterVertexSeparation));
    m_minClusterVertexSeparationSquared = m_minClusterVertexSeparation*m_minClusterVertexSeparation;

    m_minClusterDotProduct = 0.866; //  theta = 30 degrees
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "GoodClusterDotProduct", m_minClusterDotProduct));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
