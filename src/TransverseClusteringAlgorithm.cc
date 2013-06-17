/**
 *  @file   TranverseClusteringAlgorithm.cc
 * 
 *  @brief  Implementation of the transverse clustering algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArClusterHelper.h"

#include "TransverseClusteringAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode TransverseClusteringAlgorithm::Run()
{   

    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentClusterList(*this, pClusterList));

    // Cluster Vectors
    ClusterVector inputClusters, transverseClusters, longitudinalClusters, seedClusters, nonSeedClusters;


    // Sort clusters into a well-defined order
    // (using ClusterHelper tools)
    this->GetSortedClusters( pClusterList, inputClusters );
    

    // Separate transverse and longitudinal clusters
    // (transverse clusters care either short or non-longitudinal)
    this->GetTransverseClusters( inputClusters, transverseClusters, longitudinalClusters );


    // Separate transverse clusters into seed and non-seed types
    // (seed clusters are not associated with a longitudinal cluster)
    this->GetSeedClusters( transverseClusters, longitudinalClusters, seedClusters, nonSeedClusters );


    // Use seed clusters to create transverse cluster objects
    // (these are protoclusters with a direction and inner/outer vertices).
    LArTransverseClusterMap transverseClusterMap;
    this->FillTransverseClusterMap( seedClusters, transverseClusters, transverseClusterMap );


    // Form associations between transverse clusters and prepare merges
    LArClusterMergeMap clusterMergeMap;
    this->FillClusterMergeMap( transverseClusterMap, clusterMergeMap );
  



    // ---> begin display


// ClusterList tempList1, tempList2, tempList3, tempList4;
// tempList1.insert(transverseClusters.begin(), transverseClusters.end());
// tempList2.insert(longitudinalClusters.begin(), longitudinalClusters.end());
// tempList3.insert(seedClusters.begin(), seedClusters.end());
// tempList4.insert(nonSeedClusters.begin(), nonSeedClusters.end());

// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f); 
// PandoraMonitoringApi::VisualizeClusters(&tempList1, "TransverseClusters", GREEN);
// PandoraMonitoringApi::VisualizeClusters(&tempList2, "LongitudinalClusters", BLUE);
// PandoraMonitoringApi::VisualizeClusters(&tempList3, "SeedClusters", RED);
// PandoraMonitoringApi::VisualizeClusters(&tempList4, "NonSeedClusters", YELLOW);

// for ( LArTransverseClusterMap::const_iterator iterMarker = transverseClusterMap.begin(), iterEndMarker = transverseClusterMap.end(); iterMarker != iterEndMarker; ++iterMarker )
// {
//   LArTransverseCluster transverseCluster = iterMarker->second;

//   float minX = transverseCluster.GetInnerVertex().GetX();
//   float maxX = transverseCluster.GetOuterVertex().GetX();
//   float minZ = transverseCluster.GetInnerVertex().GetZ();
//   float maxZ = transverseCluster.GetOuterVertex().GetZ();

//   unsigned numMarkers = 20;

//   for( unsigned int n=0; n<numMarkers; ++n ){
//     float posX = minX + (((float)n+0.5)/((float)numMarkers))*(maxX-minX);
//     float posZ = minZ + (((float)n+0.5)/((float)numMarkers))*(maxZ-minZ);
//     CartesianVector marker(posX,0.f,posZ);
//     PANDORA_MONITORING_API(AddMarkerToVisualization(&marker, "hit", RED, 1.));
//   }    
// }

// PandoraMonitoringApi::ViewEvent();


// for (LArClusterMergeMap::const_iterator iter1 = clusterMergeMap.begin(), iterEnd1 = clusterMergeMap.end(); iter1 != iterEnd1; ++iter1)
// {
//   Cluster *pSeedCluster = iter1->first;

//   ClusterList tempList;
//   tempList.insert(pSeedCluster);

//   for( ClusterList::iterator iter2 = iter1->second.begin(), iterEnd2 = iter1->second.end(); iter2 != iterEnd2; ++iter2 )
//   {
//     Cluster *pAssociatedCluster = *iter2;
//     tempList.insert(pAssociatedCluster);
//   }
          
//   PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f); 
//   PandoraMonitoringApi::VisualizeClusters(&tempList, "NewCluster",  BLACK);
//   PandoraMonitoringApi::ViewEvent();	
// }

    // ---> end display



    // Make cluster merges
    for (LArClusterMergeMap::const_iterator iter1 = clusterMergeMap.begin(), iterEnd1 = clusterMergeMap.end(); iter1 != iterEnd1; ++iter1)
    {
        Cluster *pSeedCluster = iter1->first;
        for( ClusterList::iterator iter2 = iter1->second.begin(), iterEnd2 = iter1->second.end(); iter2 != iterEnd2; ++iter2 )
        {
            Cluster *pAssociatedCluster = *iter2;

            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pSeedCluster, pAssociatedCluster));            
        }
    }

   

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseClusteringAlgorithm::GetSortedClusters( const ClusterList *const pClusterList, ClusterVector &clusterVector )
{
    // clear vector
    clusterVector.clear();

    // loop over input cluster list
    for ( ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter )
        clusterVector.push_back(*iter);

    // sort clusters
    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNOccupiedLayers);
}
   
//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseClusteringAlgorithm::GetTransverseClusters( const ClusterVector &inputVector, ClusterVector &transverseVector, ClusterVector &longitudinalVector )
{
    // clear vectors  
    transverseVector.clear();
    longitudinalVector.clear();

    // loop over input cluster list
    for ( ClusterVector::const_iterator iter = inputVector.begin(), iterEnd = inputVector.end(); iter != iterEnd; ++iter )
    {
        Cluster* pCluster = *iter;  

        // separate short and long clusters
        if ( 1 + pCluster->GetOuterPseudoLayer() - pCluster->GetInnerPseudoLayer() <= m_clusterLayers ) 
	{
	    transverseVector.push_back(pCluster);
	}

        // separate transverse and longitudinal clusters
        else
	{
            const CartesianVector longitudinalDirection( 0.f, 0.f, 1.f );
 
            const ClusterHelper::ClusterFitResult &clusterFitResult( pCluster->GetFitToAllHitsResult() );

	    if ( clusterFitResult.IsFitSuccessful() )
	    {            
                const CartesianVector fitDirection( clusterFitResult.GetDirection() );

                // calculate dot product with respect to forward direction
                if ( std::fabs( fitDirection.GetDotProduct(longitudinalDirection) ) < m_clusterCosAngle )
                    transverseVector.push_back(pCluster);
                else
                    longitudinalVector.push_back(pCluster);
	    }
	}
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseClusteringAlgorithm::GetSeedClusters( const pandora::ClusterVector &shortVector, const pandora::ClusterVector &longVector, pandora::ClusterVector &seedVector,  pandora::ClusterVector &nonSeedVector )
{
    // clear vectors
    seedVector.clear();
    nonSeedVector.clear();

    // loop over short cluster list
    for ( ClusterVector::const_iterator iterShort = shortVector.begin(), iterEndShort = shortVector.end(); iterShort != iterEndShort; ++iterShort )
    {
        Cluster* pClusterShort = *iterShort;  

        bool isSeedCluster(true);

        if ( 1 + pClusterShort->GetOuterPseudoLayer() - pClusterShort->GetInnerPseudoLayer() <= m_clusterLayers ) 
	{
            CartesianVector position( ( pClusterShort->GetCentroid(pClusterShort->GetInnerPseudoLayer())
				      + pClusterShort->GetCentroid(pClusterShort->GetOuterPseudoLayer()) ) * 0.5 );
                
            // loop over long cluster list
            for ( ClusterVector::const_iterator iterLong = longVector.begin(), iterEndLong = longVector.end(); iterLong != iterEndLong; ++iterLong )
            {
                Cluster* pClusterLong = *iterLong;  

	        if( LArClusterHelper::GetClosestDistance( position, pClusterLong ) < m_clusterWindow )
	        {
		    isSeedCluster = false;
                    break;
		}
	    }
	}

        if ( isSeedCluster ) seedVector.push_back( pClusterShort );
        else                 nonSeedVector.push_back( pClusterShort );
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseClusteringAlgorithm::GetTransverseAssociatedClusters( const Cluster* inputCluster, const ClusterVector &inputClusterVector, ClusterVector &outputClusterVector )
{  
    // clear vector 
    outputClusterVector.clear();

    // input cluster
    Cluster* pCluster = (Cluster*)inputCluster;  // John won't approve of this...

    // long clusters 
    if ( 1 + pCluster->GetOuterPseudoLayer() - pCluster->GetInnerPseudoLayer() > m_clusterLayers ) 
    {
        outputClusterVector.push_back(pCluster);
        return;
    }
    
    // short clusters
    for ( ClusterVector::const_iterator iter = inputClusterVector.begin(), iterEnd = inputClusterVector.end(); iter != iterEnd; ++iter )
    {
        Cluster* pNearbyCluster = *iter;  

        if ( 1 + pCluster->GetOuterPseudoLayer() - pCluster->GetInnerPseudoLayer() > m_clusterLayers
          || 1 + pNearbyCluster->GetOuterPseudoLayer() - pNearbyCluster->GetInnerPseudoLayer() > m_clusterLayers ) 
	    continue;

        if ( pCluster == pNearbyCluster )
	    continue;

        if ( this->IsTransverseAssociated( pCluster, pNearbyCluster ) ) 
            outputClusterVector.push_back( pNearbyCluster );
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseClusteringAlgorithm::CollectAssociatedClusters(Cluster *pSeedCluster, Cluster *pCurrentCluster, LArClusterMergeMap &clusterMergeMap, LArClusterVetoMap &clusterVetoMap, ClusterList &associatedClusterList)
{
    LArClusterVetoMap::const_iterator iter0 = clusterVetoMap.find(pCurrentCluster);

    if (iter0 != clusterVetoMap.end())
        return;

    LArClusterMergeMap::const_iterator iter1 = clusterMergeMap.find(pCurrentCluster);

    if (iter1 == clusterMergeMap.end())
        return;

    for (ClusterList::const_iterator iter2 = iter1->second.begin(), iterEnd2 = iter1->second.end(); iter2 != iterEnd2; ++iter2)
    {
        Cluster* pAssociatedCluster = *iter2; 

        if (pAssociatedCluster == pSeedCluster)
            continue;

        if (!associatedClusterList.insert(pAssociatedCluster).second)
            continue;

        this->CollectAssociatedClusters(pSeedCluster, pAssociatedCluster, clusterMergeMap, clusterVetoMap, associatedClusterList);
    }

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseClusteringAlgorithm::FillTransverseClusterMap( const ClusterVector &seedClusters, const ClusterVector &transverseClusters, LArTransverseClusterMap &transverseClusterMap )
{ 
    ClusterVector associatedClusters;

    for ( ClusterVector::const_iterator iter = seedClusters.begin(), iterEnd = seedClusters.end(); iter != iterEnd; ++iter )
    {  
        const Cluster* pCluster = *iter;

        associatedClusters.clear();

        this->GetTransverseAssociatedClusters( pCluster, transverseClusters, associatedClusters );

        if ( associatedClusters.size()==0 ) continue;

        transverseClusterMap.insert( std::pair<const Cluster*,LArTransverseCluster>(pCluster, LArTransverseCluster(pCluster,associatedClusters)) );
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseClusteringAlgorithm::FillClusterMergeMap( const LArTransverseClusterMap &transverseClusterMap, LArClusterMergeMap &outputMergeMap )
{

    // Find associations between clusters
    LArClusterMergeMap firstMergeMap;

    for ( LArTransverseClusterMap::const_iterator iter1 = transverseClusterMap.begin(), iterEnd1 = transverseClusterMap.end(); iter1 != iterEnd1; ++iter1 )
    {
        LArTransverseCluster transCluster1 = iter1->second;

        for ( LArTransverseClusterMap::const_iterator iter2 = iter1, iterEnd2 = iterEnd1; iter2 != iterEnd2; ++iter2 )
        {
            LArTransverseCluster transCluster2 = iter2->second;

            if ( transCluster1.GetCluster() == transCluster2.GetCluster() ) 
                continue;

            if ( this->IsTransverseAssociated( transCluster1, transCluster2 ) == false ) 
                continue;
	    
            for( unsigned int n1 = 0; n1 < transCluster1.GetNumClusters(); ++n1 )
	    {
                Cluster* pCluster1 = (Cluster*)(transCluster1.GetCluster(n1));

                for( unsigned int n2 = 0; n2 < transCluster2.GetNumClusters(); ++n2 )
		{		        
		    Cluster* pCluster2 = (Cluster*)(transCluster2.GetCluster(n2));

                    if( pCluster1 == pCluster2 ) continue;
		    
                    if( firstMergeMap[pCluster1].find(pCluster2) == firstMergeMap[pCluster1].end() ) 
                        firstMergeMap[pCluster1].insert(pCluster2);

                    if( firstMergeMap[pCluster2].find(pCluster1) == firstMergeMap[pCluster2].end() ) 
                        firstMergeMap[pCluster2].insert(pCluster1);
		}
	    }
	}
    }   
 
    // Collect sets of associated clusters
    LArClusterMergeMap secondMergeMap;
    LArClusterVetoMap clusterVetoMap;
    
    for (LArClusterMergeMap::const_iterator iter1 = firstMergeMap.begin(), iterEnd1 = firstMergeMap.end(); iter1 != iterEnd1; ++iter1)
    {
        Cluster *pSeedCluster = iter1->first;

        ClusterList associationList;
        this->CollectAssociatedClusters(pSeedCluster, pSeedCluster, firstMergeMap, clusterVetoMap, associationList);

        for( ClusterList::iterator iter2 = associationList.begin(), iterEnd2 = associationList.end(); iter2 != iterEnd2; ++iter2 )
        {
            Cluster *pAssociatedCluster = *iter2;

            secondMergeMap[pSeedCluster].insert(pAssociatedCluster);
	}

        for( ClusterList::iterator iter2 = associationList.begin(), iterEnd2 = associationList.end(); iter2 != iterEnd2; ++iter2 )
        {
            Cluster *pAssociatedCluster = *iter2;

            clusterVetoMap[pAssociatedCluster] = true;
        }
    }


    // Select good transverse clusters
    for (LArClusterMergeMap::const_iterator iter1 = secondMergeMap.begin(), iterEnd1 = secondMergeMap.end(); iter1 != iterEnd1; ++iter1)
    {
        Cluster *pSeedCluster = iter1->first;

        if( iter1->second.size()<m_minTransverseLayers ) 
            continue;

        ClusterVector associatedClusters;
        for( ClusterList::iterator iter2 = iter1->second.begin(), iterEnd2 = iter1->second.end(); iter2 != iterEnd2; ++iter2 )
	    associatedClusters.push_back(*iter2);

        LArTransverseCluster possibleNewCluster(pSeedCluster,associatedClusters);

        if( (possibleNewCluster.GetOuterVertex()-possibleNewCluster.GetInnerVertex()).GetMagnitudeSquared() < m_minTransverseLength*m_minTransverseLength )
	    continue;


	for( ClusterList::iterator iter2 = iter1->second.begin(), iterEnd2 = iter1->second.end(); iter2 != iterEnd2; ++iter2 )
	{
            Cluster *pAssociatedCluster = *iter2;

            outputMergeMap[pSeedCluster].insert(pAssociatedCluster);
	}
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TransverseClusteringAlgorithm::IsTransverseAssociated( const pandora::Cluster* pCluster1, const pandora::Cluster* pCluster2 )
{

    //
    // compare centroids
    //

    float aveX1( 0.5 * ( pCluster1->GetCentroid( pCluster1->GetInnerPseudoLayer() ).GetX()
		       + pCluster1->GetCentroid( pCluster1->GetOuterPseudoLayer() ).GetX() ) );

    float aveZ1( 0.5 * ( pCluster1->GetCentroid( pCluster1->GetInnerPseudoLayer() ).GetZ()
		       + pCluster1->GetCentroid( pCluster1->GetOuterPseudoLayer() ).GetZ() ) );

    float aveX2( 0.5 * ( pCluster2->GetCentroid( pCluster2->GetInnerPseudoLayer() ).GetX()
		       + pCluster2->GetCentroid( pCluster2->GetOuterPseudoLayer() ).GetX() ) );

    float aveZ2( 0.5 * ( pCluster2->GetCentroid( pCluster2->GetInnerPseudoLayer() ).GetZ()
		       + pCluster2->GetCentroid( pCluster2->GetOuterPseudoLayer() ).GetZ() ) );


    if ( std::fabs(aveX2 - aveX1) < m_clusterWindow && std::fabs(aveZ2 - aveZ1) < m_clusterWindow
      && std::fabs(aveZ2 - aveZ1) < std::fabs(aveX2 - aveX1) * std::fabs(m_clusterTanAngle) )
        return true;
        

    return false;
}



//------------------------------------------------------------------------------------------------------------------------------------------

bool TransverseClusteringAlgorithm::IsTransverseAssociated( const LArTransverseCluster &transCluster1, const LArTransverseCluster &transCluster2 )
{
    // relative angle  
    const CartesianVector& direction1(transCluster1.GetDirection());
    const CartesianVector& direction2(transCluster2.GetDirection());

    if( direction1.GetDotProduct(direction2) < m_minCosRelativeAngle )
        return false;

    // relative distance
    const CartesianVector& innerVertex1(transCluster1.GetInnerVertex());
    const CartesianVector& innerVertex2(transCluster2.GetInnerVertex());

    if( this->IsTransverseAssociated( transCluster1, innerVertex2 ) == false
     && this->IsTransverseAssociated( transCluster2, innerVertex1 ) == false )
        return false;

    const CartesianVector& outerVertex1(transCluster1.GetOuterVertex());
    const CartesianVector& outerVertex2(transCluster2.GetOuterVertex());

    if( this->IsTransverseAssociated( transCluster1, outerVertex2 ) == false
     && this->IsTransverseAssociated( transCluster2, outerVertex1 ) == false )
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TransverseClusteringAlgorithm::IsTransverseAssociated( const LArTransverseCluster &transCluster, const pandora::CartesianVector& testVertex )
{
  

    const CartesianVector& innerVertex = transCluster.GetInnerVertex();
    const CartesianVector& outerVertex = transCluster.GetOuterVertex();
    const CartesianVector& direction   = transCluster.GetDirection();
    
    if( std::fabs( direction.GetCrossProduct(testVertex - innerVertex).GetMagnitudeSquared() ) > m_maxTransverseSeparation*m_maxTransverseSeparation )
        return false;
    
    if( direction.GetDotProduct(testVertex - innerVertex) < -m_maxLongitudinalSeparation 
     || direction.GetDotProduct(testVertex - outerVertex) > +m_maxLongitudinalSeparation )
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

TransverseClusteringAlgorithm::LArTransverseCluster::LArTransverseCluster(const Cluster *pSeedCluster, const ClusterVector &associatedClusters) :
    m_innerVertex(0.f,0.f,0.f),
    m_outerVertex(0.f,0.f,0.f),
    m_direction(0.f,0.f,1.f),
    m_rms(std::numeric_limits<float>::max())
{
    float Swzz(0.f);  
    float Swxx(0.f);  
    float Swzx(0.f);
    float Swz(0.f);
    float Swx(0.f);
    float Sw(0.f);

    float minX(-999.f);
    float maxX(-999.f);

    float aveX(0.f); 
    float aveZ(0.f);

    Cluster* pCluster = (Cluster*)pSeedCluster;  // John won't approve of this...
    m_clusterVector.push_back(pCluster);

    for ( ClusterVector::const_iterator iter = associatedClusters.begin(), iterEnd = associatedClusters.end(); iter != iterEnd; ++iter )
    {
        Cluster* pCluster = *iter;

        if( pCluster == pSeedCluster ) continue;

        m_clusterVector.push_back(*iter);
    }

    for ( ClusterVector::const_iterator iterI = m_clusterVector.begin(), iterEndI = m_clusterVector.end(); iterI != iterEndI; ++iterI )
    {
        const Cluster* pInputCluster = *iterI;  

        const OrderedCaloHitList &orderedCaloHitList(pInputCluster->GetOrderedCaloHitList());

        for (OrderedCaloHitList::const_iterator iterJ = orderedCaloHitList.begin(), iterEndJ = orderedCaloHitList.end(); iterJ != iterEndJ; ++iterJ )
        {
            const CaloHitList *pCaloHitList = iterJ->second;

            for (CaloHitList::const_iterator iterK = pCaloHitList->begin(), iterEndK = pCaloHitList->end(); iterK != iterEndK; ++iterK )
            {
	        const CaloHit* pCaloHit = *iterK;

                if( minX<0.f || pCaloHit->GetPositionVector().GetX()<minX ) 
                    minX = pCaloHit->GetPositionVector().GetX();

                if( maxX<0.f || pCaloHit->GetPositionVector().GetX()>maxX ) 
                    maxX = pCaloHit->GetPositionVector().GetX();

                Swzz += pCaloHit->GetPositionVector().GetZ() * pCaloHit->GetPositionVector().GetZ();
                Swxx += pCaloHit->GetPositionVector().GetX() * pCaloHit->GetPositionVector().GetX();
                Swzx += pCaloHit->GetPositionVector().GetZ() * pCaloHit->GetPositionVector().GetX();
                Swz  += pCaloHit->GetPositionVector().GetZ();
                Swx  += pCaloHit->GetPositionVector().GetX();
                Sw   += 1.f;
            }
        }
    }

    if ( Sw > 0.f )
    { 
        aveX = Swx/Sw; 
        aveZ = Swz/Sw;
     
        if ( Sw * Swxx - Swx * Swx > 0.f )
        {
            float m( (Sw * Swzx - Swx * Swz) / (Sw * Swxx - Swx * Swx) );
            float c( (Swz - m * Swx) / Sw );

            float px( 1.f / std::sqrt(1.f + m * m) );
            float pz( m / std::sqrt(1.f + m * m) );

            m_innerVertex.SetValues( minX, 0.f, aveZ + m*(minX-aveX) );
            m_outerVertex.SetValues( maxX, 0.f, aveZ + m*(maxX-aveX) );
            m_direction.SetValues( px, 0.f, pz );

            m_rms = std::sqrt( ((Swzz + m * m * Swxx + c * c * Sw) - 2.f * (m * Swzx + c * Swz - m * c * Swx)) / Sw );
	}
        else
	{
            m_innerVertex.SetValues( aveX, 0.f, aveZ );
	    m_outerVertex.SetValues( aveX, 0.f, aveZ );
            m_direction.SetValues( 1.f, 0.f, 0.f);
            //m_rms = std::numeric_limits<float>::max();
	}
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TransverseClusteringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_clusterWindow = 3.0;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterWindow", m_clusterWindow));

    m_clusterLayers = 5;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterLayers", m_clusterLayers));

    m_clusterAngle = 45.0;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterAngle", m_clusterAngle));

    m_clusterCosAngle = std::cos( m_clusterAngle * M_PI / 180.f );
    m_clusterTanAngle = std::tan( m_clusterAngle * M_PI / 180.f );

    m_minCosRelativeAngle = 0.866;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCosRelativeAngle", m_minCosRelativeAngle));

    m_maxTransverseSeparation = 1.5;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTransverseSeparation", m_maxTransverseSeparation));

    m_maxLongitudinalSeparation = 7.5;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxLongitudinalSeparation", m_maxLongitudinalSeparation));

    m_minTransverseLength = 5.0;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinTransverseLength", m_minTransverseLength));

    m_minTransverseLayers = 5;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinTransverseLayers", m_minTransverseLayers));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
