/**
 *  @file   PandoraPFANew/Framework/src/Pandora/VertexFindingAlgorithm.cc
 * 
 *  @brief  Implementation of the cluster creation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "VertexFindingAlgorithm.h"

#include "LArVertexHelper.h"
#include "LArGeometryHelper.h"
#include "LArClusterHelper.h"
#include "LArPointingClusterHelper.h"

#include <fstream>
#include <cmath>

using namespace pandora;

namespace lar
{

StatusCode VertexFindingAlgorithm::Run()
{
    // Cheat the vertex
    if ( m_useTrueVertex ) return SetTrueVertex();


    // Get the cluster lists for each view
    const ClusterList *pClusterListU = NULL;    
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetClusterList(*this, m_clusterListNameU, pClusterListU));

    const ClusterList *pClusterListV = NULL;    
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetClusterList(*this, m_clusterListNameV, pClusterListV));

    const ClusterList *pClusterListW = NULL;    
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetClusterList(*this, m_clusterListNameW, pClusterListW));


    // Generate list of clean pointing clusters
    LArPointingClusterMap pointingClusterMapU, pointingClusterMapV, pointingClusterMapW;

    this->GetListOfCleanPointingClusters( pClusterListU, pointingClusterMapU );
    this->GetListOfCleanPointingClusters( pClusterListV, pointingClusterMapV );
    this->GetListOfCleanPointingClusters( pClusterListW, pointingClusterMapW );

    
    // Select a list of vertex clusters
    LArPointingClusterVertexList pointingVertexListU, pointingVertexListV, pointingVertexListW;
    
    this->GetListOfCleanVertexClusters( pointingClusterMapU, pointingVertexListU );
    this->GetListOfCleanVertexClusters( pointingClusterMapV, pointingVertexListV );
    this->GetListOfCleanVertexClusters( pointingClusterMapW, pointingVertexListW );


    // Use selected vertex clusters to generate candidate vertex positions
    CartesianPointList candidateVertexListU, candidateVertexListV,  candidateVertexListW; 

    this->GetListOfCandidateVertexPositions( pointingClusterMapU, pointingVertexListU, candidateVertexListU );
    this->GetListOfCandidateVertexPositions( pointingClusterMapV, pointingVertexListV, candidateVertexListV );
    this->GetListOfCandidateVertexPositions( pointingClusterMapW, pointingVertexListW, candidateVertexListW );



   
   





    
    CartesianPointList matchedVertexListU, matchedVertexListV,  matchedVertexListW; 

    this->GetListOfMatchedVertexPositions( pointingClusterMapU,  pointingClusterMapV,  pointingClusterMapW, 
                                           candidateVertexListU, candidateVertexListV, candidateVertexListW,
                                           matchedVertexListU,   matchedVertexListV,   matchedVertexListW );





 
   



// 
// 
// =====================================================


    // Process individual views
    VertexFigureOfMeritMap theFigureOfMeritMapU;
    VertexFigureOfMeritMap theFigureOfMeritMapV;
    VertexFigureOfMeritMap theFigureOfMeritMapW;

    ProcessSingleView( pointingClusterMapU, pointingVertexListU, theFigureOfMeritMapU ); 

    ProcessSingleView( pointingClusterMapV, pointingVertexListV, theFigureOfMeritMapV ); 

    ProcessSingleView( pointingClusterMapW, pointingVertexListW, theFigureOfMeritMapW ); 




    // Process individual vertices
    CartesianVector recoVertexU(0.f,0.f,0.f);
    CartesianVector recoVertexV(0.f,0.f,0.f);
    CartesianVector recoVertexW(0.f,0.f,0.f);

    ProcessSingleVertex( pClusterListU, theFigureOfMeritMapU, recoVertexU ); 

    ProcessSingleVertex( pClusterListV, theFigureOfMeritMapV, recoVertexV ); 

    ProcessSingleVertex( pClusterListW, theFigureOfMeritMapW, recoVertexW ); 




    // Big Nested Loop 
    CartesianVector mergedVertexU(0.f,0.f,0.f);
    CartesianVector mergedVertexV(0.f,0.f,0.f);
    CartesianVector mergedVertexW(0.f,0.f,0.f);

    CartesianVector bestMergedVertexU(0.f,0.f,0.f);
    CartesianVector bestMergedVertexV(0.f,0.f,0.f);
    CartesianVector bestMergedVertexW(0.f,0.f,0.f);

    float chiSquared(0.f);
    float mergedFigureOfMerit(0.f);
    float bestMergedFigureOfMerit(-99999.f);

    float theMagicNumber(0.1);

    for( VertexFigureOfMeritMap::const_iterator iterU = theFigureOfMeritMapU.begin(), iterEndU = theFigureOfMeritMapU.end(); iterU != iterEndU; ++iterU )
    {
        const CartesianVector vertexU = (iterU->first)->GetPosition();
        float          figureOfMeritU = iterU->second;

        for( VertexFigureOfMeritMap::const_iterator iterV = theFigureOfMeritMapV.begin(), iterEndV = theFigureOfMeritMapV.end(); iterV != iterEndV; ++iterV )
        {
            const CartesianVector vertexV = (iterV->first)->GetPosition();
            float          figureOfMeritV = iterV->second;

            for( VertexFigureOfMeritMap::const_iterator iterW = theFigureOfMeritMapW.begin(), iterEndW = theFigureOfMeritMapW.end(); iterW != iterEndW; ++iterW )
            {
                const CartesianVector vertexW = (iterW->first)->GetPosition();
                float          figureOfMeritW = iterW->second;

                LArGeometryHelper::MergeThreeViews( vertexU, vertexV, vertexW,
                                                    mergedVertexU, mergedVertexV, mergedVertexW,
                                                    chiSquared );

                mergedFigureOfMerit = figureOfMeritU + figureOfMeritV + figureOfMeritW - theMagicNumber * chiSquared;

                if( mergedFigureOfMerit>bestMergedFigureOfMerit )
                {
                    bestMergedVertexU = mergedVertexU;
                    bestMergedVertexV = mergedVertexV;
                    bestMergedVertexW = mergedVertexW;
                    bestMergedFigureOfMerit = mergedFigureOfMerit;
		}
	    }
	}
    }
        

    // Clean up
    this->CleanUp( theFigureOfMeritMapU ); 
    
    this->CleanUp( theFigureOfMeritMapV ); 

    this->CleanUp( theFigureOfMeritMapW ); 

   
    // Set vertices
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, SetVertex( bestMergedVertexU, m_vertexNameU ));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, SetVertex( bestMergedVertexV, m_vertexNameV ));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, SetVertex( bestMergedVertexW, m_vertexNameW ));



    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexFindingAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector)
{
    for ( ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter ) 
    {
        Cluster *pCluster = *iter;

        if( LArClusterHelper::GetLayerSpan(pCluster) < 10 )
            continue;

        if( LArClusterHelper::GetLengthSquared(pCluster) < 15.f )
            continue;

        if ( LArClusterHelper::GetLayerOccupancy(pCluster) < 0.75 ) 
            continue;

        clusterVector.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexFindingAlgorithm::GetListOfCleanPointingClusters(const ClusterList *const pClusterList, LArPointingClusterMap& pointingClusterMap )
{
    // Select a set of clean clusters
    ClusterVector clusterVector;
    this->GetListOfCleanClusters(pClusterList, clusterVector);
    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNOccupiedLayers);

    // Generate a list of clean pointing clusters
    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
        pointingClusterMap.insert( std::pair<Cluster*,LArPointingCluster>(*iter,LArPointingCluster(*iter)) );
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexFindingAlgorithm::GetListOfCleanVertexClusters( const LArPointingClusterMap& pointingClusterMap, LArPointingClusterVertexList& cleanVertexList )
{ 
    // Build lists of pointing clusters
    LArPointingClusterList cleanClusterList;
    LArPointingClusterList fullClusterList;
    LArPointingClusterList possibleClusterList;
    
    float totalLength(0.f);

    for (LArPointingClusterMap::const_iterator iter = pointingClusterMap.begin(), iterEnd = pointingClusterMap.end(); iter != iterEnd; ++iter)
    {
        const LArPointingCluster& cluster = iter->second;

        totalLength += LArClusterHelper::GetLength( cluster.GetCluster() );
    }

    for (LArPointingClusterMap::const_iterator iterI = pointingClusterMap.begin(), iterEndI = pointingClusterMap.end(); iterI != iterEndI; ++iterI)
    {
        const LArPointingCluster& clusterI = iterI->second;

        // select clean clusters
        if( clusterI.GetLength() > 10.f 
         || LArClusterHelper::GetLength( clusterI.GetCluster() ) > 10.f
	 || LArClusterHelper::GetLength( clusterI.GetCluster() ) > 0.1 * totalLength )
        {
            fullClusterList.push_back(clusterI);
            cleanClusterList.push_back(clusterI);

            continue;
	}

        // separate vertex-like and shower-like clusters
        bool isVertexLike(false), isShowerLike(false);

        for (LArPointingClusterMap::const_iterator iterJ = pointingClusterMap.begin(), iterEndJ = pointingClusterMap.end(); iterJ != iterEndJ; ++iterJ)
        {
            const LArPointingCluster& clusterJ = iterJ->second;

            if ( clusterI.GetCluster() == clusterJ.GetCluster() ) continue; 
            
            if ( this->IsAdjacent(clusterJ.GetInnerVertex(), clusterI.GetInnerVertex().GetPosition()) == true
	      || this->IsAdjacent(clusterJ.GetOuterVertex(), clusterI.GetOuterVertex().GetPosition()) == true )
	    {
	        isVertexLike = true;
            }

	    if ( this->IsAdjacent(clusterJ.GetInnerVertex(), clusterI.GetInnerVertex().GetPosition()) == true
	      && this->IsAdjacent(clusterJ.GetOuterVertex(), clusterI.GetOuterVertex().GetPosition()) == true )
	    {
	        isShowerLike = true;
            }
	}

        if ( true == isVertexLike && false == isShowerLike )
	{
            fullClusterList.push_back(clusterI);
            possibleClusterList.push_back(clusterI);
	}
    }


    // Build lists of vertex clusters
    LArPointingClusterVertexList fullVertexList;
    LArPointingClusterVertexList possibleVertexList;

    for (LArPointingClusterList::const_iterator iter = fullClusterList.begin(), iterEnd = fullClusterList.end(); iter != iterEnd; ++iter)
    {
        const LArPointingCluster& cluster = *iter;
        fullVertexList.push_back(cluster.GetInnerVertex());
        fullVertexList.push_back(cluster.GetOuterVertex());
    }

    for (LArPointingClusterList::const_iterator iter = possibleClusterList.begin(), iterEnd = possibleClusterList.end(); iter != iterEnd; ++iter)
    {
        const LArPointingCluster& cluster = *iter;
        possibleVertexList.push_back(cluster.GetInnerVertex());
        possibleVertexList.push_back(cluster.GetOuterVertex());
    }

    for (LArPointingClusterList::const_iterator iter = cleanClusterList.begin(), iterEnd = cleanClusterList.end(); iter != iterEnd; ++iter)
    {
        const LArPointingCluster& cluster = *iter;
        cleanVertexList.push_back(cluster.GetInnerVertex());
        cleanVertexList.push_back(cluster.GetOuterVertex());
    }


    // Separate vertex-like and shower-like clusters
    for (LArPointingClusterVertexList::const_iterator iterI = possibleVertexList.begin(), iterEndI = possibleVertexList.end(); iterI != iterEndI; ++iterI)
    {
        const LArPointingCluster::Vertex& vertexI = *iterI;
        
        bool isVertexLike(false), isShowerLike(false);

        for (LArPointingClusterVertexList::const_iterator iterJ = fullVertexList.begin(), iterEndJ = fullVertexList.end(); iterJ != iterEndJ; ++iterJ)
        { 
            const LArPointingCluster::Vertex& vertexJ = *iterJ;
            
            if ( vertexI.GetCluster() == vertexJ.GetCluster() ) continue;
	   
            if ( this->IsAdjacent(vertexJ, vertexI.GetPosition()) == true
	      && this->IsDownstream(vertexJ, vertexI.GetPosition()) == true )
	    {
                isVertexLike = true;  
	    }

	    if ( this->IsNode(vertexJ.GetPosition(), vertexJ.GetDirection(), vertexI.GetPosition()) == true )
	    {
	        isVertexLike = true;  
	    }

            if ( this->IsAdjacent(vertexI, vertexJ.GetPosition()) == true
	      && this->IsAdjacent(vertexJ, vertexI.GetPosition()) == true 
	      && this->IsNode(vertexJ.GetPosition(), vertexJ.GetDirection(), vertexI.GetPosition()) == false
	      && this->IsNode(vertexI.GetPosition(), vertexI.GetDirection(), vertexJ.GetPosition()) == false )
	    {
	        isShowerLike = true;
	    }
        }

        if ( true == isVertexLike && false == isShowerLike )
        {
            cleanVertexList.push_back(vertexI);  
        }
    }

// ClusterList inputVertexList, outputVertexList;
// CollectClusters( fullVertexList, inputVertexList );
// CollectClusters( cleanVertexList, outputVertexList );
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&inputVertexList, "CLEAN CLUSTERS, FIRST PASS", GREEN);
// PandoraMonitoringApi::VisualizeClusters(&outputVertexList, "CLEAN CLUSTERS, SECOND PASS", BLUE);
// PandoraMonitoringApi::ViewEvent();

}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexFindingAlgorithm::GetListOfCandidateVertexPositions( const LArPointingClusterMap& pointingClusterMap, const LArPointingClusterVertexList& cleanClusterList, CartesianPointList& candidateVertexList )
{  
    // Put primary vertex positions into output list
    for (LArPointingClusterVertexList::const_iterator iter = cleanClusterList.begin(), iterEnd = cleanClusterList.end(); iter != iterEnd; ++iter)
        candidateVertexList.push_back((*iter).GetPosition());


    // Generate extra vertex positions by taking intersections between clusters
    for (LArPointingClusterVertexList::const_iterator iterI = cleanClusterList.begin(), iterEndI = cleanClusterList.end(); iterI != iterEndI; ++iterI)
    {
        const LArPointingCluster::Vertex& vertexI = *iterI;

        bool foundIntersect(false);
        float closestDisplacement(10.f);
        CartesianVector closestPosition(0.f,0.f,0.f);
        Cluster* closestCluster(NULL);

        for (LArPointingClusterVertexList::const_iterator iterJ = cleanClusterList.begin(), iterEndJ = cleanClusterList.end(); iterJ != iterEndJ; ++iterJ)
        {
            const LArPointingCluster::Vertex& vertexJ = *iterJ;

            if( vertexI.GetCluster() == vertexJ.GetCluster() ) continue;
   
            // Use closest end of each target cluster
            LArPointingClusterMap::const_iterator lookupJ = pointingClusterMap.find( vertexJ.GetCluster() );

            if ( lookupJ == pointingClusterMap.end() ) 
                throw pandora::StatusCodeException(STATUS_CODE_NOT_ALLOWED);

            const LArPointingCluster& clusterJ = lookupJ->second;
            const LArPointingCluster::Vertex& innerJ = clusterJ.GetInnerVertex();
            const LArPointingCluster::Vertex& outerJ = clusterJ.GetOuterVertex();

            const float innerDistanceSquared = (innerJ.GetPosition() - vertexI.GetPosition()).GetMagnitudeSquared();
            const float outerDistanceSquared = (outerJ.GetPosition() - vertexI.GetPosition()).GetMagnitudeSquared();

            if ( ( vertexJ.IsInner() == true  && innerDistanceSquared > outerDistanceSquared ) 
	      || ( vertexJ.IsInner() == false && innerDistanceSquared < outerDistanceSquared ) ) 
                continue;


            // Primary vertex must not be nodally associated with target cluster
            if ( this->IsAdjacent( vertexJ, vertexI.GetPosition() ) == true 
	      || this->IsNode( vertexJ.GetPosition(), vertexJ.GetDirection(), vertexI.GetPosition() ) == true )
	    {
	        foundIntersect = false;  break;
	    }


            // Check that vertex trajectory will pass through cluster
            if ( innerJ.GetDirection().GetDotProduct(vertexI.GetPosition()-innerJ.GetPosition()) < 0
	      || outerJ.GetDirection().GetDotProduct(vertexI.GetPosition()-outerJ.GetPosition()) < 0 ) 
                continue;


            // Calculate intersection between vertex trajectory and target cluster
            bool isPhysical(false);
            float intersectDisplacement(0.f);
            CartesianVector intersectPosition(0.f,0.f,0.f);

            this->GetIntersection( vertexI, vertexJ.GetCluster(), 
                                   intersectPosition, intersectDisplacement, isPhysical );

            if( false == isPhysical
             || intersectDisplacement > closestDisplacement ) 
                continue;
            

            // New vertex is closest to the original vertex
            if( intersectDisplacement < closestDisplacement )
	    {
                closestDisplacement = intersectDisplacement;
                foundIntersect = false;
	    }


            // New vertex must not be nodally associated with target cluster
            if( this->IsNode( vertexJ.GetPosition(), vertexJ.GetDirection(), intersectPosition ) == false 
             && this->IsNode( vertexI.GetPosition(), vertexI.GetDirection(), intersectPosition ) == false )
	    {
	        closestPosition = intersectPosition;  
                closestCluster = vertexJ.GetCluster();
                foundIntersect = true;
	    }
	}

        if( foundIntersect ) 
	{
            candidateVertexList.push_back( closestPosition );

// ClusterList tempListI, tempListJ;
// tempListI.insert( vertexI.GetCluster() );
// tempListJ.insert( closestCluster );
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&tempListI, "PRIMARY CLUSTER", GREEN);
// PandoraMonitoringApi::VisualizeClusters(&tempListJ, "TARGET CLUSTER", YELLOW);
// PandoraMonitoringApi::AddMarkerToVisualization(&vertexI.GetPosition(), "PRIMARY VERTEX", BLUE, 1.75);
// PandoraMonitoringApi::AddMarkerToVisualization(&closestPosition, "INTERSECTION", RED, 1.75);
// PandoraMonitoringApi::ViewEvent();

	}
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexFindingAlgorithm::GetListOfCleanVertexPositions( const LArPointingClusterMap& pointingClusterMap, const CartesianPointList& inputList, CartesianPointList& outputList )
{
    // Calculate sum of squared lengths for this view
    float totalLengthSquared(0.f);

    for (LArPointingClusterMap::const_iterator iter = pointingClusterMap.begin(), iterEnd = pointingClusterMap.end(); iter != iterEnd; ++iter)
    {
        const LArPointingCluster& cluster = iter->second;

        totalLengthSquared += LArClusterHelper::GetLengthSquared( cluster.GetCluster() );
    }


    // Select candidate vertex positions
    for ( CartesianPointList::const_iterator iter = inputList.begin(), iterEnd = inputList.end(); iter != iterEnd; ++iter )
    {
        const CartesianVector seedPosition = *iter;

        float thisLengthSquared(0.f);

        for (LArPointingClusterMap::const_iterator iter = pointingClusterMap.begin(), iterEnd = pointingClusterMap.end(); iter != iterEnd; ++iter)
        {
            const LArPointingCluster& cluster = iter->second;
            const LArPointingCluster::Vertex& innerVertex = cluster.GetInnerVertex();
            const LArPointingCluster::Vertex& outerVertex = cluster.GetOuterVertex();

            if ( this->IsNode( innerVertex.GetPosition(), innerVertex.GetDirection(), seedPosition )
	      || this->IsNode( outerVertex.GetPosition(), outerVertex.GetDirection(), seedPosition )
              || this->IsAdjacent( innerVertex, seedPosition ) || this->IsAdjacent( outerVertex, seedPosition ) )
	    {
	        thisLengthSquared += LArClusterHelper::GetLengthSquared( cluster.GetCluster() );
	    }
	}

        if ( thisLengthSquared > 100.f || thisLengthSquared > 0.1 * totalLengthSquared )
	{
            outputList.push_back(seedPosition);
	}
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexFindingAlgorithm::GetListOfMatchedVertexPositions( const LArPointingClusterMap& pointingClusterMapU, const LArPointingClusterMap& pointingClusterMapV, const LArPointingClusterMap& pointingClusterMapW, const CartesianPointList& inputListU, const CartesianPointList& inputListV, const CartesianPointList& inputListW, CartesianPointList& outputListU, CartesianPointList& outputListV, CartesianPointList& outputListW )
{

    CartesianPointList cleanListU, cleanListV, cleanListW;
    CartesianPointList matchedListU, matchedListV, matchedListW;
    CartesianPointList possibleListU, possibleListV, possibleListW;
        
    this->GetListOfCleanVertexPositions( pointingClusterMapU, inputListU, cleanListU );
    this->GetListOfCleanVertexPositions( pointingClusterMapV, inputListV, cleanListV );
    this->GetListOfCleanVertexPositions( pointingClusterMapW, inputListW, cleanListW );

    this->GetListOfMatchedVertexPositions2D( VIEW_U, VIEW_V, cleanListU, cleanListV, cleanListW, matchedListW, possibleListW );
    this->GetListOfMatchedVertexPositions2D( VIEW_V, VIEW_W, cleanListV, cleanListW, cleanListU, matchedListU, possibleListU );
    this->GetListOfMatchedVertexPositions2D( VIEW_W, VIEW_U, cleanListW, cleanListU, cleanListV, matchedListV, possibleListV );

    this->GetListOfCleanVertexPositions( pointingClusterMapU, possibleListU, matchedListU );
    this->GetListOfCleanVertexPositions( pointingClusterMapV, possibleListV, matchedListV );
    this->GetListOfCleanVertexPositions( pointingClusterMapW, possibleListW, matchedListW );

    this->GetListOfMatchedVertexPositions3D( matchedListU, matchedListV, matchedListW,
                                             outputListU, outputListV, outputListW );

    
/*

CartesianVector trueVertexU(130.5f, 0.f, 151.f);  // 128.2f, 0.f, 151.f 
CartesianVector trueVertexV(130.5f, 0.f, 151.f);  // 128.2f, 0.f, 151.f
CartesianVector trueVertexW(130.5f, 0.f, 100.f);  // 128.2f, 0.f, 100.f 


ClusterList clusterListU;
CollectClusters( pointingClusterMapU, clusterListU );

PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
PandoraMonitoringApi::VisualizeClusters(&clusterListU, "CLEAN CLUSTERS (U)", GREEN);
PandoraMonitoringApi::AddMarkerToVisualization(&trueVertexU, "vertexU", BLUE, 1.75);

for( CartesianPointList::const_iterator iter = matchedListU.begin(), iterEnd = matchedListU.end(); iter != iterEnd; ++iter ){
  const CartesianVector& outputVertex = *iter;
  PandoraMonitoringApi::AddMarkerToVisualization(&outputVertex, "vertexU", RED, 1.5);
}

PandoraMonitoringApi::ViewEvent();


ClusterList clusterListV;
CollectClusters( pointingClusterMapV, clusterListV );

PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
PandoraMonitoringApi::VisualizeClusters(&clusterListV, "CLEAN CLUSTERS (V)", GREEN);
PandoraMonitoringApi::AddMarkerToVisualization(&trueVertexV, "vertexV", BLUE, 1.75);

for( CartesianPointList::const_iterator iter = matchedListV.begin(), iterEnd = matchedListV.end(); iter != iterEnd; ++iter ){
  const CartesianVector& outputVertex = *iter;
  PandoraMonitoringApi::AddMarkerToVisualization(&outputVertex, "vertexV", RED, 1.5);
}

PandoraMonitoringApi::ViewEvent();


ClusterList clusterListW;
CollectClusters( pointingClusterMapW, clusterListW );

PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
PandoraMonitoringApi::VisualizeClusters(&clusterListW, "CLEAN CLUSTERS (W)", GREEN);
PandoraMonitoringApi::AddMarkerToVisualization(&trueVertexW, "vertexW", BLUE, 1.75);

for( CartesianPointList::const_iterator iter = matchedListW.begin(), iterEnd = matchedListW.end(); iter != iterEnd; ++iter ){
  const CartesianVector& outputVertex = *iter;
  PandoraMonitoringApi::AddMarkerToVisualization(&outputVertex, "vertexW", RED, 1.5);
}

PandoraMonitoringApi::ViewEvent();

*/

}


//------------------------------------------------------------------------------------------------------------------------------------------

void VertexFindingAlgorithm::GetListOfMatchedVertexPositions2D( const HitType view1, const HitType view2, const CartesianPointList& inputList1, const CartesianPointList& inputList2, const CartesianPointList& inputList3, CartesianPointList& matchedList3, CartesianPointList& projectedList3 )
{
    float m_maxSeparation = 2.5;
    float m_maxSeparationSquared = m_maxSeparation * m_maxSeparation;

    CartesianVector projectedPosition(0.f,0.f,0.f);

    float chi2(0.f);

    for (CartesianPointList::const_iterator iter1 = inputList1.begin(), iterEnd1 = inputList1.end(); iter1 != iterEnd1; ++iter1 )
    {
        const CartesianVector& position1 = *iter1;

        for (CartesianPointList::const_iterator iter2 = inputList2.begin(), iterEnd2 = inputList2.end(); iter2 != iterEnd2; ++iter2 )
	{
            const CartesianVector& position2 = *iter2;

            if ( std::fabs(position1.GetX()-position2.GetX()) > m_maxSeparation ) continue;
 
	    LArGeometryHelper::MergeTwoViews( view1, view2, position1, position2, projectedPosition, chi2 );
            

            bool foundMatch(false);
            float minSeparationSquared(m_maxSeparationSquared);
            CartesianVector matchedPosition(0.f,0.f,0.f);

            for (CartesianPointList::const_iterator iter3 = inputList3.begin(), iterEnd3 = inputList3.end(); iter3 != iterEnd3; ++iter3 )
	    {
                const CartesianVector& position3 = *iter3;

                float thisSeparationSquared((position3 - projectedPosition).GetMagnitudeSquared());

                if (thisSeparationSquared < minSeparationSquared)
		{
		    minSeparationSquared = thisSeparationSquared;
                    matchedPosition = position3;
                    foundMatch = true;
		}
	    }

            if ( foundMatch ) matchedList3.push_back(matchedPosition);
            else              projectedList3.push_back(projectedPosition);
	}
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexFindingAlgorithm::GetListOfMatchedVertexPositions3D( const CartesianPointList& inputListU, const CartesianPointList& inputListV, const CartesianPointList& inputListW, CartesianPointList& outputListU, CartesianPointList& outputListV, CartesianPointList& outputListW )
{
    float m_maxSeparation = 2.5;
    float m_maxSeparationSquared = m_maxSeparation * m_maxSeparation;

    CartesianVector projectedPositionU(0.f,0.f,0.f);
    CartesianVector projectedPositionV(0.f,0.f,0.f);
    CartesianVector projectedPositionW(0.f,0.f,0.f);

    float chi2(0.f);

    for (CartesianPointList::const_iterator iterU = inputListU.begin(), iterEndU = inputListU.end(); iterU != iterEndU; ++iterU )
    {
        const CartesianVector& positionU = *iterU;

        for (CartesianPointList::const_iterator iterV = inputListV.begin(), iterEndV = inputListV.end(); iterV != iterEndV; ++iterV )
	{
            const CartesianVector& positionV = *iterV;

            for (CartesianPointList::const_iterator iterW = inputListW.begin(), iterEndW = inputListW.end(); iterW != iterEndW; ++iterW )
	    {
                const CartesianVector& positionW = *iterW;

                if ( std::fabs(positionU.GetX()-positionV.GetX()) > m_maxSeparation 
                  || std::fabs(positionV.GetX()-positionW.GetX()) > m_maxSeparation 
                  || std::fabs(positionW.GetX()-positionU.GetX()) > m_maxSeparation ) 
                    continue;

		LArGeometryHelper::MergeThreeViews( positionU, positionV, positionW,
                                                    projectedPositionU, projectedPositionV, projectedPositionW, chi2 );

                if ( (positionU - projectedPositionU).GetMagnitudeSquared() > m_maxSeparationSquared
                  || (positionV - projectedPositionV).GetMagnitudeSquared() > m_maxSeparationSquared
		  || (positionW - projectedPositionW).GetMagnitudeSquared() > m_maxSeparationSquared )
		    continue;

                outputListU.push_back(positionU);
                outputListV.push_back(positionV);
                outputListW.push_back(positionW);
	    }
	}
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------


//
//
// NO MAN'S LAND...
//
//


//------------------------------------------------------------------------------------------------------------------------------------------

void VertexFindingAlgorithm::CollectMarkers( const LArPointingClusterVertexList& inputList, CartesianPointList& outputList )
{
    outputList.clear();

    for (LArPointingClusterVertexList::const_iterator iter = inputList.begin(), iterEnd = inputList.end(); iter != iterEnd; ++iter )
        outputList.push_back((*iter).GetPosition());
}
 
//------------------------------------------------------------------------------------------------------------------------------------------
 
void VertexFindingAlgorithm::CollectClusters( const LArPointingClusterVertexList& inputList, ClusterList& outputList )
{
    outputList.clear();

    for (LArPointingClusterVertexList::const_iterator iter = inputList.begin(), iterEnd = inputList.end(); iter != iterEnd; ++iter )
        outputList.insert((*iter).GetCluster());
}

//------------------------------------------------------------------------------------------------------------------------------------------
 
void VertexFindingAlgorithm::CollectClusters( const LArPointingClusterMap& inputList, ClusterList& outputList )
{
    outputList.clear();

    for (LArPointingClusterMap::const_iterator iter = inputList.begin(), iterEnd = inputList.end(); iter != iterEnd; ++iter )
        outputList.insert((iter->second).GetCluster());
}

//------------------------------------------------------------------------------------------------------------------------------------------


//
//
// NO MAN'S LAND...
//
//



//------------------------------------------------------------------------------------------------------------------------------------------

void VertexFindingAlgorithm::CleanUp( VertexFigureOfMeritMap& figureOfMeritMap )
{
    for( VertexFigureOfMeritMap::const_iterator iter = figureOfMeritMap.begin(), iterEnd = figureOfMeritMap.end(); iter != iterEnd; ++iter ){
        const LArVertexCandidate* thisVertex  = iter->first;
        if( thisVertex ) delete thisVertex;
    }
}




void VertexFindingAlgorithm::ProcessSingleView( const LArPointingClusterMap& pointingClusterMap, const LArPointingClusterVertexList& pointingVertexList, VertexFigureOfMeritMap& outputFigureOfMeritMap )
{

    // find connected vertices
    this->FindPossibleConnectedVertices( pointingClusterMap, pointingVertexList, outputFigureOfMeritMap );


   
}



void VertexFindingAlgorithm::ProcessSingleVertex( const ClusterList* const pClusterList, const VertexFigureOfMeritMap theFigureOfMeritMap, CartesianVector& bestVertex )
{


    

    // loop over figure of merit list
    float bestFigureOfMerit(0.f);
   
    for( VertexFigureOfMeritMap::const_iterator iter = theFigureOfMeritMap.begin(), iterEnd = theFigureOfMeritMap.end(); iter != iterEnd; ++iter ){
      const LArVertexCandidate* thisVertex  = iter->first;
      float              thisFigureOfMerit = iter->second;

      if( thisFigureOfMerit>bestFigureOfMerit ){
        bestFigureOfMerit = thisFigureOfMerit;
        bestVertex        = thisVertex->GetPosition();
      }
    }

}



void VertexFindingAlgorithm::FindPossibleConnectedVertices( const LArPointingClusterMap& pointingClusterMap, const LArPointingClusterVertexList& pointingClusterVertexCandidateList, VertexFigureOfMeritMap& outputFigureOfMeritMap )
{

    // Start by calculating the total energy 
    float totalEnergy(0.f);
    float totalEnergySquared(0.f); 
    
    for (LArPointingClusterVertexList::const_iterator iter = pointingClusterVertexCandidateList.begin(), iterEnd = pointingClusterVertexCandidateList.end(); iter != iterEnd; ++iter)
    {
        const LArPointingCluster::Vertex &thisCluster = *iter;
        const Cluster* pCluster = thisCluster.GetCluster();
        const float thisEnergy = LArClusterHelper::GetEnergyFromLength( pCluster );
      
        totalEnergy += 0.5 * thisEnergy;
        totalEnergySquared += 0.5 * thisEnergy * thisEnergy;
    }
     
    // Loop over clean clusters to identify possible vertex positions    
    for (LArPointingClusterVertexList::const_iterator iter = pointingClusterVertexCandidateList.begin(), iterEnd = pointingClusterVertexCandidateList.end(); iter != iterEnd; ++iter)
    {
        const LArPointingCluster::Vertex &clusterVertex = *iter;

        // Calculate energy of this cluster
        const Cluster* pCluster = clusterVertex.GetCluster();
        const float thisEnergy = LArClusterHelper::GetEnergyFromLength( pCluster );
        const float thisEnergySquared = thisEnergy * thisEnergy;
        
	// Find associated clusters
        float primaryEnergy(0.f);
        float associatedEnergy(0.f);
        float associatedMomentumModulus(0.f);
        unsigned int primaryClusters(0);
        CartesianVector associatedMomentum(0.f,0.f,0.f);

        CartesianVector testDirection(0.f,0.f,1.f);

        float energyFraction(0.f), momentumFraction(0.f);

        LArPointingClusterVertexList associatedClusterVertexList;
        RunFastReconstruction( pointingClusterMap, clusterVertex.GetPosition(), testDirection,//clusterVertex.GetDirection(), 
                               associatedClusterVertexList, primaryClusters, primaryEnergy, 
                               associatedEnergy, associatedMomentum );

        associatedMomentumModulus = associatedMomentum.GetMagnitude();

        if (totalEnergy > 0.f )
	{
            energyFraction = associatedEnergy/totalEnergy;
            momentumFraction = associatedMomentumModulus/totalEnergy;
	}

        // Some quality requirements on candidate vertex
        if ( associatedClusterVertexList.empty() )
	    continue;

        if ( m_runBeamMode && this->IsConsistentWithBeamDirection(associatedMomentum) == false ) 
            continue;

        if ( primaryEnergy / totalEnergy < 0.05 * ( 5.f - static_cast<float>(primaryClusters) )
	 && thisEnergySquared / totalEnergySquared < 0.25 ) 
            continue; 

        // Fill the output list
        float thisFigureOfMerit(0.f);

        if( totalEnergy > 0.f )
	    thisFigureOfMerit = energyFraction * momentumFraction;
	
        outputFigureOfMeritMap.insert( std::pair<const LArVertexCandidate*,float>
	    ( new LArVertexCandidate(clusterVertex.GetPosition(),clusterVertex.GetDirection(),energyFraction,momentumFraction),thisFigureOfMerit) );
       	      
// std::cout << "  thisVertex: energyFrac=" << thisEnergy/totalEnergy << " energySquaredFrac=" << thisEnergySquared/totalEnergySquared << " numPrimaryClusters=" << primaryClusters << " primaryEnergyFrac=" << primaryEnergy/totalEnergy << " associatedEnergyFrac=" << associatedEnergy/totalEnergy << " FigureOfMerit=" << (associatedEnergy*associatedMomentumModulus)/(totalEnergy*totalEnergy) << std::endl;
// CartesianVector thisVertex = clusterVertex.GetPosition();
// ClusterList clusterList;
// CollectClusters( associatedClusterVertexList, clusterList );
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&clusterList, "CLUSTERS", GREEN);
// PandoraMonitoringApi::AddMarkerToVisualization(&thisVertex, "VERTEX", RED, 1.);
// PandoraMonitoringApi::ViewEvent();

    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexFindingAlgorithm::RunFastReconstruction( const LArPointingClusterMap& inputMap, const CartesianVector& seedVertexPosition, const CartesianVector& seedVertexDirection, LArPointingClusterVertexList& outputList, unsigned int& primaryClusters, float& primaryEnergy, float& outputEnergy, CartesianVector& outputMomentum )
{
    outputList.clear();

    primaryClusters = 0;
    primaryEnergy = 0.f;

    outputEnergy = 0.f;
    outputMomentum.SetValues(0.f,0.f,0.f);

    LArPointingClusterVertexList fullList;
    LArPointingClusterVertexList associatedList;
    LArPointingClusterVertexList strongList;
    LArPointingClusterVertexList weakList;

    for (LArPointingClusterMap::const_iterator iter = inputMap.begin(), iterEnd = inputMap.end(); iter != iterEnd; ++iter)
    {
        const LArPointingCluster& pointingCluster = iter->second;
        fullList.push_back(pointingCluster.GetInnerVertex());
        fullList.push_back(pointingCluster.GetOuterVertex());
    }

    for (LArPointingClusterVertexList::const_iterator iter1 = fullList.begin(), iterEnd1 = fullList.end(); iter1 != iterEnd1; ++iter1)
    {
        const LArPointingCluster::Vertex& clusterVertex = *iter1;

        LArPointingClusterMap::const_iterator iter2 = inputMap.find(clusterVertex.GetCluster());

        if ( iter2 == inputMap.end() ) 
            throw pandora::StatusCodeException(STATUS_CODE_NOT_ALLOWED);

        const LArPointingCluster& cluster = iter2->second;

        const LArPointingCluster::Vertex& innerVertex = cluster.GetInnerVertex();
        const LArPointingCluster::Vertex& outerVertex = cluster.GetOuterVertex();

        const float innerDistanceSquared = (innerVertex.GetPosition() - seedVertexPosition ).GetMagnitudeSquared();
        const float outerDistanceSquared = (outerVertex.GetPosition() - seedVertexPosition ).GetMagnitudeSquared();

        if ( ( clusterVertex.IsInner()==true  && innerDistanceSquared < outerDistanceSquared ) 
	  || ( clusterVertex.IsInner()==false && innerDistanceSquared > outerDistanceSquared ) )
	    associatedList.push_back(clusterVertex);
    }

    for (LArPointingClusterVertexList::const_iterator iter0 = associatedList.begin(), iterEnd0 = associatedList.end(); iter0 != iterEnd0; ++iter0)
    {
        const LArPointingCluster::Vertex &clusterVertex = *iter0;

        const float thisEnergy( LArClusterHelper::GetEnergyFromLength( clusterVertex.GetCluster() ) );

        if ( this->IsPrimary( clusterVertex, seedVertexPosition, seedVertexDirection ) )
	{
	    strongList.push_back(clusterVertex);  
            ++primaryClusters;  primaryEnergy += thisEnergy;  
	}
        else
	{
            weakList.push_back(clusterVertex);
	}
    }

    for (LArPointingClusterVertexList::const_iterator iter1 = weakList.begin(), iterEnd1 = weakList.end(); iter1 != iterEnd1; ++iter1 )
    {
        const LArPointingCluster::Vertex &weakVertex = *iter1; 

        bool isAssociated(false);

        for (LArPointingClusterVertexList::const_iterator iter2 = strongList.begin(), iterEnd2 = strongList.end(); iter2 != iterEnd2; ++iter2 )
        {
	   const LArPointingCluster::Vertex &strongVertex = *iter2;

           if ( this->IsWeakAssociation( strongVertex, weakVertex ) )
	   {
	     isAssociated = true;   break;
	   }
	}

        if ( isAssociated ) outputList.push_back(weakVertex);
    }

    for (LArPointingClusterVertexList::const_iterator iter2 = strongList.begin(), iterEnd2 = strongList.end(); iter2 != iterEnd2; ++iter2 )
        outputList.push_back(*iter2);
    

    return this->GetEnergyAndMomentum( outputList, seedVertexPosition, outputEnergy, outputMomentum );
}
 
//------------------------------------------------------------------------------------------------------------------------------------------

void VertexFindingAlgorithm::GetEnergyAndMomentum( const LArPointingClusterVertexList& clusterList, const CartesianVector& thisVertex, float& outputEnergy, CartesianVector& outputMomentum )
{   
    outputEnergy = 0.f;
    outputMomentum.SetValues(0.f,0.f,0.f);

    for (LArPointingClusterVertexList::const_iterator iter = clusterList.begin(), iterEnd = clusterList.end(); iter != iterEnd; ++iter )
    {
        const LArPointingCluster::Vertex &pointingVertex = *iter;

	float thisEnergy(0.f);
        CartesianVector thisMomentum(0.f,0.f,0.f);

        this->GetEnergyAndMomentum( pointingVertex, thisVertex, thisEnergy, thisMomentum );

        outputEnergy   += thisEnergy;
	outputMomentum += thisMomentum;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexFindingAlgorithm::GetEnergyAndMomentum( const LArPointingCluster::Vertex& thisCluster, const CartesianVector& vertexPosition, float& outputEnergy, CartesianVector& outputMomentum )
{
    outputEnergy = 0.f;
    outputMomentum.SetValues(0.f,0.f,0.f);

    const Cluster* pCluster = thisCluster.GetCluster();
    
    const float clusterEnergy = LArClusterHelper::GetEnergyFromLength( thisCluster.GetCluster() );
    const float clusterLength = LArClusterHelper::GetLength( thisCluster.GetCluster() );

    const CartesianVector endPosition = thisCluster.GetPosition() + thisCluster.GetDirection() * clusterLength;
    const CartesianVector vertexDirection = (endPosition-vertexPosition).GetUnitVector();

    const float vertexLength = (endPosition-vertexPosition).GetMagnitude();
    const float cosTheta = thisCluster.GetDirection().GetDotProduct(vertexDirection);
    const float effectiveLength = cosTheta * std::min( vertexLength, clusterLength );

    if ( effectiveLength > 0.f )
    {
        outputEnergy = clusterEnergy * effectiveLength / clusterLength;
    }
    else 
    {
        outputEnergy = 0.f;
    }
  
    outputMomentum = thisCluster.GetDirection() * outputEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------
 
void VertexFindingAlgorithm::GetIntersection( const LArPointingCluster::Vertex& firstCluster, const LArPointingCluster::Vertex& secondCluster, CartesianVector& intersectPosition, CartesianVector& intersectDirection,  bool& isPhysical )
{
    isPhysical = false;

    const float cosRelativeAngle( (firstCluster.GetDirection()).GetDotProduct(secondCluster.GetDirection()) );

    const float maxCosRelativeAngle( 0.999 ); // 2.5 degrees (almost parallel)  
    const float maxDisplacement( 2.5 );

    float firstDisplacement(0.f), secondDisplacement(0.f);

    if ( std::fabs(cosRelativeAngle) < maxCosRelativeAngle ) // good pointing
    {
        LArPointingClusterHelper::GetIntersection( firstCluster, secondCluster, intersectPosition,
                                                   firstDisplacement, secondDisplacement );
	LArPointingClusterHelper::GetAverageDirection( firstCluster, secondCluster, intersectDirection );

        isPhysical = ( firstDisplacement < maxDisplacement && secondDisplacement < maxDisplacement );
    }

    /*
    // Possible treatment of anti-parallel clusters
    else if ( cosRelativeAngle < -maxCosRelativeAngle ) // anti-parallel (bad)
    {
        const float maxTransverseWidthAntiParallel( 1.f );
        const CartesianVector averagePosition = ( firstCluster.GetPosition() + secondCluster.GetPosition() ) * 0.5;
        const CartesianVector averageDirection = firstCluster.GetDirection(); // don't care, pick a safe value
        const float transverseWidth = averageDirection.GetCrossProduct(secondCluster.GetPosition() - firstCluster.GetPosition()).GetMagnitude(); 

        intersectPosition = averagePosition;
        intersectDirection = averageDirection;

        firstDisplacement = firstCluster.GetDirection().GetDotProduct(intersectPosition - firstCluster.GetPosition());
        secondDisplacement = secondCluster.GetDirection().GetDotProduct(intersectPosition - secondCluster.GetPosition());
        isPhysical = ( transverseWidth < maxTransverseWidthAntiParallel );
    }

    // Possible treatment of parallel clusters
    else if ( cosRelativeAngle > +maxCosRelativeAngle ) // parallel (good)
    {
        const float maxTransverseWidthParallel( 10.f );
        const CartesianVector averagePosition = (firstCluster.GetPosition() + secondCluster.GetPosition()) * 0.5;
        const CartesianVector averageDirection = (firstCluster.GetDirection() + secondCluster.GetDirection()).GetUnitVector();
        const float transverseWidth = averageDirection.GetCrossProduct(secondCluster.GetPosition() - firstCluster.GetPosition()).GetMagnitude(); 
        const float transverseLength = std::max( averageDirection.GetDotProduct(averagePosition - firstCluster.GetPosition()),
                                                 averageDirection.GetDotProduct(averagePosition - secondCluster.GetPosition()) );
        intersectPosition = averagePosition - averageDirection * (20.f * transverseWidth + transverseLength + 2.5 );
        intersectDirection = averageDirection;

        firstDisplacement = firstCluster.GetDirection().GetDotProduct(intersectPosition - firstCluster.GetPosition());
        secondDisplacement = secondCluster.GetDirection().GetDotProduct(intersectPosition - secondCluster.GetPosition());
        isPhysical = ( transverseWidth < maxTransverseWidthParallel );
    }
    */
}
 
//------------------------------------------------------------------------------------------------------------------------------------------

void VertexFindingAlgorithm::GetIntersection( const LArPointingCluster::Vertex& vertexCluster, const Cluster* pTargetCluster, CartesianVector& intersectPosition, float& intersectDisplacement, bool& isPhysical )
{            
    isPhysical = false;

    float displacementL(0.f), displacementT(0.f);

    LArPointingClusterHelper::GetIntersection( vertexCluster, pTargetCluster, intersectPosition, displacementL, displacementT );

    if (displacementT < +1.f && displacementL > -1.f) isPhysical = true;

    intersectDisplacement = displacementL;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexFindingAlgorithm::IsStrongAssociation( const LArPointingCluster::Vertex& parentCluster, const LArPointingCluster::Vertex& daughterCluster )
{
    return ( this->IsNode(parentCluster,daughterCluster) || this->IsEmission(parentCluster,daughterCluster) );
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexFindingAlgorithm::IsWeakAssociation( const LArPointingCluster::Vertex& parentCluster, const LArPointingCluster::Vertex& daughterCluster )
{
    return ( this->IsStrongSecondary(parentCluster,daughterCluster) || this->IsWeakSecondary(parentCluster,daughterCluster) );
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexFindingAlgorithm::IsNode( const LArPointingCluster::Vertex& parentCluster, const LArPointingCluster::Vertex& daughterCluster )
{
    return this->IsNode( parentCluster.GetPosition(), parentCluster.GetDirection(), daughterCluster.GetPosition(), daughterCluster.GetDirection() );
}
  
//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexFindingAlgorithm::IsEmission( const LArPointingCluster::Vertex& parentCluster, const LArPointingCluster::Vertex& daughterCluster )
{
    return this->IsEmission( parentCluster.GetPosition(), parentCluster.GetDirection(), daughterCluster.GetPosition(), daughterCluster.GetDirection() );
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexFindingAlgorithm::IsWeakSecondary( const LArPointingCluster::Vertex& parentCluster, const LArPointingCluster::Vertex& daughterCluster )
{
    if ( !this->IsDaughterDiverging(parentCluster,daughterCluster) ) return false;

    const float cosRelativeAngle( (parentCluster.GetDirection()).GetDotProduct(daughterCluster.GetDirection()) );

    if ( cosRelativeAngle < 0.866 ) return false;

    float rL(0.f), rT(0.f);

    LArVertexHelper::GetImpactParameters( parentCluster.GetPosition(), parentCluster.GetDirection() * -1.f, daughterCluster.GetPosition(), rL, rT );

    if ( rL < 2.5 ) return false;

    static const float tanTheta( std::tan( M_PI * 10.0 / 180.0 ) );

    if ( rT/rL < tanTheta )  return true;  else return false;
}
 
//------------------------------------------------------------------------------------------------------------------------------------------
   
bool VertexFindingAlgorithm::IsStrongSecondary( const LArPointingCluster::Vertex& parentCluster, const LArPointingCluster::Vertex& daughterCluster )
{ 
    const float cosRelativeAngle( (parentCluster.GetDirection()).GetDotProduct(daughterCluster.GetDirection()) );

    if ( cosRelativeAngle < 0.707 ) return false;

    float rL(0.f), rT(0.f);

    LArVertexHelper::GetImpactParameters( parentCluster.GetPosition(), parentCluster.GetDirection() * -1.f, daughterCluster.GetPosition(), rL, rT );

    if ( rL < 2.5 ) return false;

    return this->IsAdjacent( parentCluster, daughterCluster.GetPosition() );
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexFindingAlgorithm::IsDaughterDiverging( const LArPointingCluster::Vertex& parentCluster, const LArPointingCluster::Vertex& daughterCluster )
{
   const float cosRelativeAngle( (parentCluster.GetDirection()).GetDotProduct(daughterCluster.GetDirection()) );

   if ( cosRelativeAngle > 0.966 ) return true;

   if ( cosRelativeAngle < 0.25 ) return false;

   const CartesianVector& parentPosition = parentCluster.GetPosition();
   const CartesianVector& parentDirection = parentCluster.GetDirection();

   const CartesianVector& daughterPosition = daughterCluster.GetPosition();
   const CartesianVector& daughterDirection = daughterCluster.GetDirection();

   const float daughterLength = LArClusterHelper::GetLength( daughterCluster.GetCluster() );

   const float initialSeparation = parentDirection.GetCrossProduct(daughterPosition-parentPosition).GetMagnitudeSquared();
   const float finalSeparation = parentDirection.GetCrossProduct(daughterPosition+daughterDirection*daughterLength-parentPosition).GetMagnitudeSquared();

   return ( finalSeparation > initialSeparation );
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexFindingAlgorithm::IsPrimary( const LArPointingCluster::Vertex& clusterTrajectory, const CartesianVector& targetPosition, const CartesianVector& targetDirection )
{
  return ( this->IsNode( targetPosition, targetDirection, clusterTrajectory.GetPosition(), clusterTrajectory.GetDirection() )
	|| this->IsEmission( targetPosition, targetDirection, clusterTrajectory.GetPosition(), clusterTrajectory.GetDirection() ) );
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexFindingAlgorithm::IsAdjacent( const LArPointingCluster::Vertex& clusterTrajectory, const CartesianVector& targetPosition )
{
   const Cluster* thisCluster = clusterTrajectory.GetCluster();
  
   if( LArClusterHelper::GetClosestDistance( targetPosition, thisCluster ) < 2.5 ) return true;  else return false;
  
   return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexFindingAlgorithm::IsDownstream( const LArPointingCluster::Vertex& clusterTrajectory, const CartesianVector& targetPosition )
{
    if( this->IsNode( clusterTrajectory.GetPosition(), clusterTrajectory.GetDirection(), targetPosition ) ) return true;
    
    if( (clusterTrajectory.GetPosition()-targetPosition).GetUnitVector().GetDotProduct(clusterTrajectory.GetDirection()) > 0.707 ) return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexFindingAlgorithm::IsNode( const CartesianVector& parentPosition, const CartesianVector& parentDirection, const CartesianVector& daughterPosition, const CartesianVector& daughterDirection ) const
{    
    const float cosRelativeAngle( parentDirection.GetDotProduct(daughterDirection) );

    if ( cosRelativeAngle < -0.966 ) return false;

    return ( this->IsNode( parentPosition, parentDirection, daughterPosition ) 
	  || this->IsNode( daughterPosition, daughterDirection, parentPosition ) );
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexFindingAlgorithm::IsEmission( const CartesianVector& parentPosition, const CartesianVector& parentDirection, const CartesianVector& daughterPosition, const CartesianVector& daughterDirection ) const
{   
    const float cosRelativeAngle( parentDirection.GetDotProduct(daughterDirection) );

    if ( cosRelativeAngle < 0.25 ) return false;

    return this->IsEmission( daughterPosition, daughterDirection, parentPosition );
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexFindingAlgorithm::IsNode( const CartesianVector& clusterVertex, const CartesianVector& clusterDirection, const CartesianVector& targetPosition ) const
{
    float rL(0.f), rT(0.f);

    LArVertexHelper::GetImpactParameters( clusterVertex, clusterDirection, targetPosition, rL, rT );

    if ( rL > -2.5 && rL < 2.5 && rT < 2.5 ) return true;  else return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexFindingAlgorithm::IsEmission( const CartesianVector& clusterVertex, const CartesianVector& clusterDirection, const CartesianVector& targetPosition ) const
{
    float rL(0.f), rT(0.f);

    static const float tanSqTheta( pow( tan( M_PI * 2.0 / 180.0 ), 2.0 ) );

    LArVertexHelper::GetImpactParameters( clusterVertex, clusterDirection, targetPosition, rL, rT );

    if ( rL > 0.f && rL < 75.f && rT*rT < 2.5*2.5 + rL*rL*tanSqTheta ) return true;  else return false;
}
 
//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexFindingAlgorithm::IsConsistentWithBeamDirection( const CartesianVector& thisMomentum ) const
{
    if( thisMomentum.GetUnitVector().GetZ() > -0.25 ) return true;  else return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexFindingAlgorithm::SetTrueVertex()
{
    CartesianVector trueVertexU(130.5f, 0.f, 151.f);  // 128.2f,0.f,151.f 
    CartesianVector trueVertexV(130.5f, 0.f, 151.f);  // 128.2f,0.f,151.f
    CartesianVector trueVertexW(130.5f, 0.f, 100.f);  // 128.2f,0.f,100.f 

    LArVertexHelper::AddVertex(m_vertexNameU, trueVertexU);
    LArVertexHelper::AddVertex(m_vertexNameV, trueVertexV);
    LArVertexHelper::AddVertex(m_vertexNameW, trueVertexW);
    
    LArVertexHelper::SetCurrentVertex(m_vertexNameW);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexFindingAlgorithm::SetVertex(const CartesianVector& eventVertex, std::string vertexName)
{
    LArVertexHelper::AddVertex(vertexName, eventVertex);
    LArVertexHelper::SetCurrentVertex(vertexName);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexFindingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_vertexNameU = "";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "VertexNameU", m_vertexNameU));

    m_vertexNameV = "";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "VertexNameV", m_vertexNameV));

    m_vertexNameW = "";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "VertexNameW", m_vertexNameW));

    m_clusterListNameU = "";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "ClusterListNameU", m_clusterListNameU));

    m_clusterListNameV = "";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "ClusterListNameV", m_clusterListNameV));

    m_clusterListNameW = "";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "ClusterListNameW", m_clusterListNameW));

    m_useTrueVertex = false;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "UseTrueVertex", m_useTrueVertex));

    m_runBeamMode = true; 
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "RunBeamMode", m_runBeamMode));


    return STATUS_CODE_SUCCESS;
}


} // namespace lar
