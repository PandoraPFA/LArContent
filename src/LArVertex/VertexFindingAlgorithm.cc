/**
 *  @file   LArContent/src/LArVertex/VertexFindingAlgorithm.cc
 * 
 *  @brief  Implementation of the vertex finding algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"
#include "LArHelpers/LArPointingClusterHelper.h"
#include "LArHelpers/LArVertexHelper.h"

#include "LArVertex/VertexFindingAlgorithm.h"

#include <fstream>
#include <cmath>

using namespace pandora;

namespace lar
{

StatusCode VertexFindingAlgorithm::Run()
{ 
    // Build true vertex (TODO: Replace these hard-coded numbers!)
    CartesianVector trueVertexU(130.5f, 0.f, 151.f);  // 128.2f, 0.f, 151.f 
    CartesianVector trueVertexV(130.5f, 0.f, 151.f);  // 128.2f, 0.f, 151.f
    CartesianVector trueVertexW(130.5f, 0.f, 100.f);  // 128.2f, 0.f, 100.f 
    CartesianVector trueVertex3D(130.5f, 0.f, 100.f); // 128.2f, 0.f, 100.f 

    LArVertexHelper::AddVertex("TrueVertexU", trueVertexU);
    LArVertexHelper::AddVertex("TrueVertexV", trueVertexV);
    LArVertexHelper::AddVertex("TrueVertexW", trueVertexW);
    LArVertexHelper::AddVertex("TrueVertex3D", trueVertex3D);

    // Cheat the vertex
    if ( m_useTrueVertex )
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, SetVertex( trueVertexU, m_vertexNameU ));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, SetVertex( trueVertexV, m_vertexNameV ));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, SetVertex( trueVertexW, m_vertexNameW ));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, SetVertex( trueVertex3D, m_vertexName3D ));

        return STATUS_CODE_SUCCESS;
    }


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


    // Use 3D information to form a final set of candidate vertex positions 
    CartesianPointList matchedVertexListU, matchedVertexListV,  matchedVertexListW; 

    this->GetListOfMatchedVertexPositions( pointingClusterMapU,  pointingClusterMapV,  pointingClusterMapW, 
                                           candidateVertexListU, candidateVertexListV, candidateVertexListW,
                                           matchedVertexListU,   matchedVertexListV,   matchedVertexListW );




    // Process individual views
    VertexFigureOfMeritMap theFigureOfMeritMapU;
    VertexFigureOfMeritMap theFigureOfMeritMapV;
    VertexFigureOfMeritMap theFigureOfMeritMapW;

    ProcessView1D( pointingClusterMapU, matchedVertexListU, theFigureOfMeritMapU ); 
    ProcessView1D( pointingClusterMapV, matchedVertexListV, theFigureOfMeritMapV ); 
    ProcessView1D( pointingClusterMapW, matchedVertexListW, theFigureOfMeritMapW ); 




    // Process vertices
    CartesianVector recoVertexU(0.f,0.f,0.f);
    CartesianVector recoVertexV(0.f,0.f,0.f);
    CartesianVector recoVertexW(0.f,0.f,0.f);

    /*
    ProcessVertex1D( theFigureOfMeritMapU, recoVertexU ); 
    ProcessVertex1D( theFigureOfMeritMapV, recoVertexV ); 
    ProcessVertex1D( theFigureOfMeritMapW, recoVertexW ); 
    */

    ProcessVertex3D( theFigureOfMeritMapU, theFigureOfMeritMapV, theFigureOfMeritMapW, 
                     recoVertexU, recoVertexV, recoVertexW ); 



    // Calculate 3D Vertex
    float chiSquared(0.f);
    CartesianVector recoVertex3D(0.f,0.f,0.f);
    
    LArGeometryHelper::MergeThreePositions3D(VIEW_U, VIEW_V, VIEW_W,
                                             recoVertexU, recoVertexV, recoVertexW,
                                             recoVertex3D, chiSquared);
                          


    // Clean up
    this->CleanUp( theFigureOfMeritMapU ); 
    this->CleanUp( theFigureOfMeritMapV ); 
    this->CleanUp( theFigureOfMeritMapW ); 

   
    // Set vertices
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, SetVertex( recoVertexU, m_vertexNameU ));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, SetVertex( recoVertexV, m_vertexNameV ));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, SetVertex( recoVertexW, m_vertexNameW ));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, SetVertex( recoVertex3D, m_vertexName3D ));



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
// Cluster* closestCluster(NULL); // (For event display below)

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

            if ( ( vertexJ.IsInnerVertex() == true  && innerDistanceSquared > outerDistanceSquared ) 
                || ( vertexJ.IsInnerVertex() == false && innerDistanceSquared < outerDistanceSquared ) ) 
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
// closestCluster = vertexJ.GetCluster(); // (For event display below)
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

// ClusterList clusterList;
// CollectClusters( pointingClusterMap, clusterList );
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&clusterList, "CLEAN CLUSTERS", GREEN);
// for( CartesianPointList::const_iterator iter = outputList.begin(), iterEnd = outputList.end(); iter != iterEnd; ++iter ){
//   const CartesianVector& outputVertex = *iter;
//   PandoraMonitoringApi::AddMarkerToVisualization(&outputVertex, "vertex", RED, 1.5);
// }
// PandoraMonitoringApi::ViewEvent();

}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexFindingAlgorithm::GetListOfMatchedVertexPositions( const LArPointingClusterMap& pointingClusterMapU, const LArPointingClusterMap& pointingClusterMapV, const LArPointingClusterMap& pointingClusterMapW, const CartesianPointList& inputListU, const CartesianPointList& inputListV, const CartesianPointList& inputListW, CartesianPointList& outputListU, CartesianPointList& outputListV, CartesianPointList& outputListW )
{

    CartesianPointList cleanListU, cleanListV, cleanListW;
    CartesianPointList possibleListU, possibleListV, possibleListW;
        
    this->GetListOfCleanVertexPositions( pointingClusterMapU, inputListU, cleanListU );
    this->GetListOfCleanVertexPositions( pointingClusterMapV, inputListV, cleanListV );
    this->GetListOfCleanVertexPositions( pointingClusterMapW, inputListW, cleanListW );

    this->GetListOfMatchedVertexPositions2D( VIEW_U, VIEW_V, 
                                             cleanListU, cleanListV, possibleListW );
    this->GetListOfMatchedVertexPositions2D( VIEW_V, VIEW_W, 
                                             cleanListV, cleanListW, possibleListU );
    this->GetListOfMatchedVertexPositions2D( VIEW_W, VIEW_U, 
                                             cleanListW, cleanListU, possibleListV );

    this->GetListOfCleanVertexPositions( pointingClusterMapU, possibleListU, cleanListU );
    this->GetListOfCleanVertexPositions( pointingClusterMapV, possibleListV, cleanListV );
    this->GetListOfCleanVertexPositions( pointingClusterMapW, possibleListW, cleanListW );

    this->GetListOfMatchedVertexPositions3D( cleanListU, cleanListV, cleanListW,
                                             outputListU, outputListV, outputListW );

    

// CartesianVector trueVertexU(130.5f, 0.f, 151.f);  // 128.2f, 0.f, 151.f 
// CartesianVector trueVertexV(130.5f, 0.f, 151.f);  // 128.2f, 0.f, 151.f
// CartesianVector trueVertexW(130.5f, 0.f, 100.f);  // 128.2f, 0.f, 100.f 

// ClusterList clusterListU;
// CollectClusters( pointingClusterMapU, clusterListU );
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&clusterListU, "CLEAN CLUSTERS (U)", GREEN);
// PandoraMonitoringApi::AddMarkerToVisualization(&trueVertexU, "vertexU", BLUE, 1.75);
// for( CartesianPointList::const_iterator iter = outputListU.begin(), iterEnd = outputListU.end(); iter != iterEnd; ++iter ){
//   const CartesianVector& outputVertex = *iter;
//   PandoraMonitoringApi::AddMarkerToVisualization(&outputVertex, "vertexU", RED, 1.5);
// }
// PandoraMonitoringApi::ViewEvent();

// ClusterList clusterListV;
// CollectClusters( pointingClusterMapV, clusterListV );
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&clusterListV, "OUTPUT CLUSTERS (V)", GREEN);
// PandoraMonitoringApi::AddMarkerToVisualization(&trueVertexV, "vertexV", BLUE, 1.75);
// for( CartesianPointList::const_iterator iter = outputListV.begin(), iterEnd = outputListV.end(); iter != iterEnd; ++iter ){
//   const CartesianVector& outputVertex = *iter;
//   PandoraMonitoringApi::AddMarkerToVisualization(&outputVertex, "vertexV", RED, 1.5);
// }
// PandoraMonitoringApi::ViewEvent();

// ClusterList clusterListW;
// CollectClusters( pointingClusterMapW, clusterListW );
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&clusterListW, "OUTPUT CLUSTERS (W)", GREEN);
// PandoraMonitoringApi::AddMarkerToVisualization(&trueVertexW, "vertexW", BLUE, 1.75);
// for( CartesianPointList::const_iterator iter = outputListW.begin(), iterEnd = outputListW.end(); iter != iterEnd; ++iter ){
//   const CartesianVector& outputVertex = *iter;
//   PandoraMonitoringApi::AddMarkerToVisualization(&outputVertex, "vertexW", RED, 1.5);
// }
// PandoraMonitoringApi::ViewEvent();

}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexFindingAlgorithm::GetListOfMatchedVertexPositions2D( const HitType view1, const HitType view2, const CartesianPointList& inputList1, const CartesianPointList& inputList2, CartesianPointList& projectedList3 )
{
    // Form matches between two views
    float m_maxSeparation = 2.5;
    
    CartesianVector position3(0.f,0.f,0.f);

    float chi2(0.f);

    for (CartesianPointList::const_iterator iter1 = inputList1.begin(), iterEnd1 = inputList1.end(); iter1 != iterEnd1; ++iter1 )
    {
        const CartesianVector& position1 = *iter1;

        for (CartesianPointList::const_iterator iter2 = inputList2.begin(), iterEnd2 = inputList2.end(); iter2 != iterEnd2; ++iter2 )
        {
            const CartesianVector& position2 = *iter2;

            if ( std::fabs(position1.GetX()-position2.GetX()) > m_maxSeparation ) continue;
 
            LArGeometryHelper::MergeTwoPositions( view1, view2, position1, position2, position3, chi2 );

            projectedList3.push_back(position3);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexFindingAlgorithm::GetListOfMatchedVertexPositions3D( const CartesianPointList& inputListU, const CartesianPointList& inputListV, const CartesianPointList& inputListW, CartesianPointList& outputListU, CartesianPointList& outputListV, CartesianPointList& outputListW )
{
    //
    // TODO: Match the directions as well as positions...
    //

    // Configurable parameters
    float m_maxSeparation = 2.5;
    float m_maxSeparationSquared = m_maxSeparation * m_maxSeparation;


    // Form matches between three views 
    CartesianPointList matchedListU, matchedListV, matchedListW;

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

                LArGeometryHelper::MergeThreePositions( positionU, positionV, positionW,
                                                        projectedPositionU, projectedPositionV, projectedPositionW, chi2 );

                if ( (positionU - projectedPositionU).GetMagnitudeSquared() > m_maxSeparationSquared
                    || (positionV - projectedPositionV).GetMagnitudeSquared() > m_maxSeparationSquared
                    || (positionW - projectedPositionW).GetMagnitudeSquared() > m_maxSeparationSquared )
                    continue;

                matchedListU.push_back(positionU);
                matchedListV.push_back(positionV);
                matchedListW.push_back(positionW);
            }
        }
    }

    this->GetReducedListOfVertexPositions( matchedListU, outputListU );
    this->GetReducedListOfVertexPositions( matchedListV, outputListV );
    this->GetReducedListOfVertexPositions( matchedListW, outputListW );
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexFindingAlgorithm::GetReducedListOfVertexPositions( const CartesianPointList& inputList, CartesianPointList& outputList )
{
    // Generate a more compact list
    float m_minSeparation = 0.5;
    float m_minSeparationSquared = m_minSeparation * m_minSeparation;

    for (unsigned int i=0; i<inputList.size(); ++i )
    {
        const CartesianVector& positionI = inputList.at(i);

        bool vetoThisEntry(false);

        for (unsigned int j=i+1; j<inputList.size(); ++j )
        {
            const CartesianVector& positionJ = inputList.at(j);

            if ( (positionI - positionJ).GetMagnitudeSquared() < m_minSeparationSquared )
            {
                vetoThisEntry = true;  break;
            }
        }

        if ( false == vetoThisEntry )
        {
            outputList.push_back(positionI);
        }
    }
}

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

void VertexFindingAlgorithm::CleanUp( VertexFigureOfMeritMap& figureOfMeritMap )
{
    for( VertexFigureOfMeritMap::const_iterator iter = figureOfMeritMap.begin(), iterEnd = figureOfMeritMap.end(); iter != iterEnd; ++iter ){
        const LArVertexCandidate* thisVertex  = iter->first;
        if( thisVertex ) delete thisVertex;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexFindingAlgorithm::ProcessVertex1D( const VertexFigureOfMeritMap theFigureOfMeritMap, CartesianVector& bestVertex )
{

    // loop over figure of merit list
    float bestFigureOfMerit(0.f);

    bool foundVertex(false);

    for( VertexFigureOfMeritMap::const_iterator iter = theFigureOfMeritMap.begin(), iterEnd = theFigureOfMeritMap.end(); iter != iterEnd; ++iter ){
      const LArVertexCandidate* thisVertex  = iter->first;
      float              thisFigureOfMerit = iter->second;

      if ( m_runBeamMode 
        && this->IsConsistentWithBeamDirectionW(thisVertex->GetDirection()) == false ) 
          continue;

      if( thisFigureOfMerit>bestFigureOfMerit ){
        foundVertex       = true;
        bestFigureOfMerit = thisFigureOfMerit;
        bestVertex        = thisVertex->GetPosition();
      }
    }

    if( false == foundVertex )
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexFindingAlgorithm::ProcessVertex3D( const VertexFigureOfMeritMap theFigureOfMeritMapU, const VertexFigureOfMeritMap theFigureOfMeritMapV, const VertexFigureOfMeritMap theFigureOfMeritMapW, CartesianVector& bestVertexU, CartesianVector& bestVertexV, CartesianVector& bestVertexW )
{  
    // Configurable parameters
    float m_maxSeparation = 2.5;
    float m_maxSeparationSquared = m_maxSeparation * m_maxSeparation;

    // Big Nested Loop 
    CartesianVector mergedVertexU(0.f,0.f,0.f);
    CartesianVector mergedVertexV(0.f,0.f,0.f);
    CartesianVector mergedVertexW(0.f,0.f,0.f);

    float chiSquared(0.f);
    float mergedFigureOfMerit(0.f);
    float bestFigureOfMerit(-99999.f);

    float theMagicNumber(0.1);

    bool foundVertex(false);

    for( VertexFigureOfMeritMap::const_iterator iterW = theFigureOfMeritMapW.begin(), iterEndW = theFigureOfMeritMapW.end(); iterW != iterEndW; ++iterW )
    {
        const LArVertexCandidate* candidateW = iterW->first;
        const CartesianVector& vertexW    = candidateW->GetPosition();
        const CartesianVector& directionW = candidateW->GetDirection();
        float              figureOfMeritW = iterW->second;

        if ( m_runBeamMode 
          && this->IsConsistentWithBeamDirectionW(directionW) == false ) 
            continue;

        for( VertexFigureOfMeritMap::const_iterator iterV = theFigureOfMeritMapV.begin(), iterEndV = theFigureOfMeritMapV.end(); iterV != iterEndV; ++iterV )
        {
            const LArVertexCandidate* candidateV = iterV->first;
            const CartesianVector& vertexV    = candidateV->GetPosition(); 
            const CartesianVector& directionV = candidateV->GetDirection();
            float              figureOfMeritV = iterV->second;

            for( VertexFigureOfMeritMap::const_iterator iterU = theFigureOfMeritMapU.begin(), iterEndU = theFigureOfMeritMapU.end(); iterU != iterEndU; ++iterU )
            {
                const LArVertexCandidate* candidateU = iterU->first;
                const CartesianVector& vertexU    = candidateU->GetPosition();
                const CartesianVector& directionU = candidateU->GetDirection();
                float              figureOfMeritU = iterU->second;

            if ( m_runBeamMode 
                && this->IsConsistentWithBeamDirectionUV(directionU,directionV) == false ) 
                    continue;

                //
                // TODO: Put in a requirement of consistent U,V,W directions
                //

                if ( std::fabs(vertexU.GetX()-vertexV.GetX()) > m_maxSeparation 
                  || std::fabs(vertexV.GetX()-vertexW.GetX()) > m_maxSeparation 
                  || std::fabs(vertexW.GetX()-vertexU.GetX()) > m_maxSeparation ) 
                    continue;

                LArGeometryHelper::MergeThreePositions( vertexU, vertexV, vertexW,
                                                        mergedVertexU, mergedVertexV, mergedVertexW,
                                                        chiSquared );

                if ( (vertexU - mergedVertexU).GetMagnitudeSquared() > m_maxSeparationSquared
                    || (vertexV - mergedVertexV).GetMagnitudeSquared() > m_maxSeparationSquared
                    || (vertexW - mergedVertexW).GetMagnitudeSquared() > m_maxSeparationSquared )
                    continue;

// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::AddMarkerToVisualization(&vertexU, "testVertexU", RED, 3.5);
// PandoraMonitoringApi::AddMarkerToVisualization(&vertexV, "testVertexV", BLUE, 3.5);
// PandoraMonitoringApi::AddMarkerToVisualization(&vertexW, "testVertexW", BLACK, 3.5);
// PandoraMonitoringApi::ViewEvent();

                mergedFigureOfMerit = figureOfMeritU + figureOfMeritV + figureOfMeritW - theMagicNumber * chiSquared;

                if( mergedFigureOfMerit>bestFigureOfMerit )
                {
                    foundVertex = true;
                    bestVertexU = mergedVertexU;
                    bestVertexV = mergedVertexV;
                    bestVertexW = mergedVertexW;
                    bestFigureOfMerit = mergedFigureOfMerit;
                }
            }
        }
    }

    if( false == foundVertex )
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexFindingAlgorithm::ProcessView1D( const LArPointingClusterMap& pointingClusterMap, const CartesianPointList& inputVertexList, VertexFigureOfMeritMap& outputFigureOfMeritMap )
{

    float totalEnergy(0.f);

    for (LArPointingClusterMap::const_iterator iter = pointingClusterMap.begin(), iterEnd = pointingClusterMap.end(); iter != iterEnd; ++iter)
    {
        const LArPointingCluster& cluster = iter->second;

        totalEnergy += LArClusterHelper::GetEnergyFromLength( cluster.GetCluster() );
    }


    for (CartesianPointList::const_iterator iter = inputVertexList.begin(), iterEnd = inputVertexList.end(); iter != iterEnd; ++iter )
    {
        const CartesianVector& seedVertexPosition = *iter;

         
        // 
        // TODO: Calculate a seed direction for this vertex
        //

        const CartesianVector seedVertexDirection(0.f,0.f,1.f); // beam direction

  
        
        // Find associated clusters
        float associatedEnergy(0.f);
        float associatedMomentumModulus(0.f);
        CartesianVector associatedMomentum(0.f,0.f,0.f);

        float energyFraction(0.f), momentumFraction(0.f);

        LArPointingClusterVertexList associatedClusterVertexList;
        RunFastReconstruction( pointingClusterMap, seedVertexPosition, seedVertexDirection,
                               associatedClusterVertexList, associatedEnergy, associatedMomentum );

        associatedMomentumModulus = associatedMomentum.GetMagnitude();

        if (totalEnergy > 0.f)
        {
            energyFraction = associatedEnergy/totalEnergy;
            momentumFraction = associatedMomentumModulus/totalEnergy;
        }

        // Some quality requirements on candidate vertex
        if ( associatedClusterVertexList.empty() )
            continue;

        // Fill the output list
        float thisFigureOfMerit(0.f);

        if( totalEnergy > 0.f )
            thisFigureOfMerit = energyFraction * momentumFraction;

        outputFigureOfMeritMap.insert( std::pair<const LArVertexCandidate*,float>
            ( new LArVertexCandidate(seedVertexPosition,associatedMomentum.GetUnitVector(),
              energyFraction,momentumFraction),thisFigureOfMerit) );
    }

// ClusterList clusterList;
// CollectClusters( pointingClusterMap, clusterList );
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&clusterList, "CLEAN CLUSTERS", GREEN); 
// for( VertexFigureOfMeritMap::const_iterator iter = outputFigureOfMeritMap.begin(), iterEnd = outputFigureOfMeritMap.end(); iter != iterEnd; ++iter )
// {
//   const CartesianVector thisVertex = (iter->first)->GetPosition();
//   PandoraMonitoringApi::AddMarkerToVisualization(&thisVertex, "vertex", RED, 1.5);
// }
// PandoraMonitoringApi::ViewEvent();

}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexFindingAlgorithm::RunFastReconstruction( const LArPointingClusterMap& inputMap, const CartesianVector& seedVertexPosition, const CartesianVector& seedVertexDirection, LArPointingClusterVertexList& outputList, float& outputEnergy, CartesianVector& outputMomentum )
{
    outputList.clear();

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

        if ( ( clusterVertex.IsInnerVertex()==true  && innerDistanceSquared < outerDistanceSquared ) 
            || ( clusterVertex.IsInnerVertex()==false && innerDistanceSquared > outerDistanceSquared ) )
            associatedList.push_back(clusterVertex);
    }

    for (LArPointingClusterVertexList::const_iterator iter0 = associatedList.begin(), iterEnd0 = associatedList.end(); iter0 != iterEnd0; ++iter0)
    {
        const LArPointingCluster::Vertex &clusterVertex = *iter0;

        if ( this->IsPrimary( clusterVertex, seedVertexPosition, seedVertexDirection ) )
        {
            strongList.push_back(clusterVertex);  
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

        //
        // TODO: Deal with split clusters! Here is a temporary fix...
        //

        if ( this->IsStrongSecondary( weakVertex, seedVertexPosition, seedVertexDirection ) )
            isAssociated = true;


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

    return this->IsWeakSecondary( parentCluster, daughterCluster.GetPosition(), daughterCluster.GetDirection() ); 
}
 
//------------------------------------------------------------------------------------------------------------------------------------------
   
bool VertexFindingAlgorithm::IsStrongSecondary( const LArPointingCluster::Vertex& parentCluster, const LArPointingCluster::Vertex& daughterCluster )
{ 
    return this->IsStrongSecondary( parentCluster, daughterCluster.GetPosition(), daughterCluster.GetDirection() );
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

bool VertexFindingAlgorithm::IsStrongSecondary( const LArPointingCluster::Vertex& clusterTrajectory, const CartesianVector& targetPosition, const CartesianVector& targetDirection )
{
    const float cosRelativeAngle( (clusterTrajectory.GetDirection()).GetDotProduct(targetDirection) );

    if ( cosRelativeAngle < 0.707 ) return false;

    float rL(0.f), rT(0.f);

    LArVertexHelper::GetImpactParameters( clusterTrajectory.GetPosition(), clusterTrajectory.GetDirection() * -1.f, targetPosition, rL, rT );

    if ( rL < 2.5 ) return false;

    return this->IsAdjacent( clusterTrajectory, targetPosition ); 
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexFindingAlgorithm::IsWeakSecondary( const LArPointingCluster::Vertex& clusterTrajectory, const CartesianVector& targetPosition, const CartesianVector& targetDirection )
{
    const float cosRelativeAngle( (clusterTrajectory.GetDirection()).GetDotProduct(targetDirection) );

    if ( cosRelativeAngle < 0.866 ) return false;

    float rL(0.f), rT(0.f);

    LArVertexHelper::GetImpactParameters( clusterTrajectory.GetPosition(), clusterTrajectory.GetDirection() * -1.f, targetPosition, rL, rT );

    if ( rL < 2.5 ) return false;

    static const float tanTheta( std::tan( M_PI * 10.0 / 180.0 ) );

    if ( rT/rL < tanTheta )  return true;  else return false;
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

bool VertexFindingAlgorithm::IsConsistentWithBeamDirectionW( const CartesianVector& momentumW ) const
{
    CartesianVector directionW( momentumW.GetUnitVector() );

    if ( directionW.GetZ() > -0.25 ) return true;  else return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexFindingAlgorithm::IsConsistentWithBeamDirectionUV( const CartesianVector& momentumU, const CartesianVector& momentumV ) const
{
    CartesianVector directionU( momentumU.GetUnitVector() );
    CartesianVector directionV( momentumV.GetUnitVector() );

    if ( directionU.GetX() * directionV.GetX() < 0 )
    {
        directionU.SetValues( -directionU.GetX(), 0.f, directionU.GetZ() );
    }

    CartesianVector directionW( LArGeometryHelper::MergeTwoDirections(VIEW_U,VIEW_V,directionU,directionV) );

    return this->IsConsistentWithBeamDirectionW( directionW );
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

    m_vertexName3D = "";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "VertexName3D", m_vertexName3D));

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
