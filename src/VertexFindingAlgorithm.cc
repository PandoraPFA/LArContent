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
#include "LArParticleId.h"

#include <fstream>
#include <cmath>

using namespace pandora;

namespace lar
{

StatusCode VertexFindingAlgorithm::Run()
{
    // Cheat the vertex
    if ( m_useTrueVertex ) return SetTrueVertex();


    const ClusterList *pClusterListU = NULL;    
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetClusterList(*this, m_clusterListNameU, pClusterListU));

    const ClusterList *pClusterListV = NULL;    
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetClusterList(*this, m_clusterListNameV, pClusterListV));

    const ClusterList *pClusterListW = NULL;    
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetClusterList(*this, m_clusterListNameW, pClusterListW));


    // Process individual views
    VertexFigureOfMeritList theFigureOfMeritListU;
    VertexFigureOfMeritList theFigureOfMeritListV;
    VertexFigureOfMeritList theFigureOfMeritListW;

    ProcessSingleView( pClusterListU, theFigureOfMeritListU ); 

    ProcessSingleView( pClusterListV, theFigureOfMeritListV ); 

    ProcessSingleView( pClusterListW, theFigureOfMeritListW ); 




    // Process individual vertices
    CartesianVector recoVertexU(0.f,0.f,0.f);
    CartesianVector recoVertexV(0.f,0.f,0.f);
    CartesianVector recoVertexW(0.f,0.f,0.f);

    ProcessSingleVertex( pClusterListU, theFigureOfMeritListU, recoVertexU ); 

    ProcessSingleVertex( pClusterListV, theFigureOfMeritListV, recoVertexV ); 

    ProcessSingleVertex( pClusterListW, theFigureOfMeritListW, recoVertexW ); 




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

    for( VertexFigureOfMeritList::const_iterator iterU = theFigureOfMeritListU.begin(), iterEndU = theFigureOfMeritListU.end(); iterU != iterEndU; ++iterU )
    {
        const CartesianVector vertexU = (iterU->first)->GetPosition();
        float          figureOfMeritU = iterU->second;

        for( VertexFigureOfMeritList::const_iterator iterV = theFigureOfMeritListV.begin(), iterEndV = theFigureOfMeritListV.end(); iterV != iterEndV; ++iterV )
        {
            const CartesianVector vertexV = (iterV->first)->GetPosition();
            float          figureOfMeritV = iterV->second;

            for( VertexFigureOfMeritList::const_iterator iterW = theFigureOfMeritListW.begin(), iterEndW = theFigureOfMeritListW.end(); iterW != iterEndW; ++iterW )
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
    this->CleanUp( theFigureOfMeritListU ); 
    
    this->CleanUp( theFigureOfMeritListV ); 

    this->CleanUp( theFigureOfMeritListW ); 

   
    // Set vertices
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, SetVertex( bestMergedVertexU, m_vertexNameU ));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, SetVertex( bestMergedVertexV, m_vertexNameV ));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, SetVertex( bestMergedVertexW, m_vertexNameW ));



    return STATUS_CODE_SUCCESS;
}

void VertexFindingAlgorithm::CleanUp( VertexFigureOfMeritList& figureOfMeritList )
{
    for( VertexFigureOfMeritList::const_iterator iter = figureOfMeritList.begin(), iterEnd = figureOfMeritList.end(); iter != iterEnd; ++iter ){
        const LArPointingVertex* thisVertex  = iter->first;
        if( thisVertex ) delete thisVertex;
    }
}

void VertexFindingAlgorithm::ProcessSingleView( const pandora::ClusterList* const pClusterList, VertexFigureOfMeritList& outputFigureOfMeritList )
{
 
    // Select a set of clean clusters
    ClusterVector clusterVector;
    this->GetListOfCleanClusters(pClusterList, clusterVector);
    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNOccupiedLayers);


    // Generate a list of clean pointing clusters
    LArPointingClusterMap pointingClusterMap;
    LArPointingClusterList pointingClusterList;
    LArPointingClusterVertexList pointingClusterVertexList;
    LArPointingClusterVertexList pointingClusterVertexCandidateList;

    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
    {
        pointingClusterMap.insert( std::pair<Cluster*,LArPointingCluster>(*iter,LArPointingCluster(*iter)) );
    }

    for (LArPointingClusterMap::const_iterator iter = pointingClusterMap.begin(), iterEnd = pointingClusterMap.end(); iter != iterEnd; ++iter)
    {
        const LArPointingCluster& pointingCluster = iter->second;
        pointingClusterList.push_back(pointingCluster);
        pointingClusterVertexList.push_back(pointingCluster.GetInnerVertex());
        pointingClusterVertexList.push_back(pointingCluster.GetOuterVertex());
    }

    this->GetListOfCleanVertexClusters( pointingClusterVertexList, pointingClusterVertexCandidateList );


    //ClusterList inputClusterVertexList, outputClusterVertexList;
    //CollectClusters( pointingClusterVertexList, inputClusterVertexList );
    //CollectClusters( pointingClusterVertexCandidateList, outputClusterVertexList );
    //PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
    //PandoraMonitoringApi::VisualizeClusters(&inputClusterVertexList, "CLEAN CLUSTERS, FIRST PASS", GREEN);
    //PandoraMonitoringApi::VisualizeClusters(&outputClusterVertexList, "CLEAN CLUSTERS, SECOND PASS", BLUE);
    //PandoraMonitoringApi::ViewEvent();

 
    // find connected vertices
    this->FindPossibleConnectedVertices( pointingClusterMap, pointingClusterVertexCandidateList, outputFigureOfMeritList );

    // find displaced vertices
    this->FindPossibleDisplacedVertices( pointingClusterMap, pointingClusterVertexCandidateList, outputFigureOfMeritList );

}

void VertexFindingAlgorithm::ProcessSingleVertex( const ClusterList* const pClusterList, const VertexFigureOfMeritList theFigureOfMeritList, pandora::CartesianVector& bestVertex )
{


    // set first pass vertex    
    GetFirstPassVertex(pClusterList,bestVertex); 

    // loop over figure of merit list
    float bestFigureOfMerit(0.f);
   
    for( VertexFigureOfMeritList::const_iterator iter = theFigureOfMeritList.begin(), iterEnd = theFigureOfMeritList.end(); iter != iterEnd; ++iter ){
      const LArPointingVertex* thisVertex  = iter->first;
      float              thisFigureOfMerit = iter->second;

      if( thisFigureOfMerit>bestFigureOfMerit ){
        bestFigureOfMerit = thisFigureOfMerit;
        bestVertex        = thisVertex->GetPosition();
      }
    }

}

void VertexFindingAlgorithm::GetFirstPassVertex( const ClusterList* const pClusterList, CartesianVector& firstVertex )
{
    if ( pClusterList->empty() ) return;
        

    ClusterVector clusterVector;
    for ( ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter ) 
        clusterVector.push_back(*iter);
    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNOccupiedLayers);
    
    const Cluster* pCluster = clusterVector.at(0);

    if ( m_runBeamMode )
    {
        firstVertex = pCluster->GetCentroid( pCluster->GetInnerPseudoLayer() );
    }
    else
    {
        firstVertex = ( pCluster->GetCentroid(pCluster->GetInnerPseudoLayer())
		      + pCluster->GetCentroid(pCluster->GetOuterPseudoLayer()) ) * 0.5;
    }

    
}

void VertexFindingAlgorithm::FindPossibleConnectedVertices( const LArPointingClusterMap& pointingClusterMap, const LArPointingClusterVertexList& pointingClusterVertexCandidateList, VertexFigureOfMeritList& outputFigureOfMeritList )
{

    // Start by calculating the total energy 
    float totalEnergy(0.f);
    float totalEnergySquared(0.f); 
    
    for (LArPointingClusterVertexList::const_iterator iter = pointingClusterVertexCandidateList.begin(), iterEnd = pointingClusterVertexCandidateList.end(); iter != iterEnd; ++iter)
    {
        const LArPointingCluster::Vertex &thisCluster = *iter;
        const Cluster* pCluster = thisCluster.GetCluster();
        const float thisEnergy = GetEnergy( pCluster );
      
        totalEnergy += 0.5 * thisEnergy;
        totalEnergySquared += 0.5 * thisEnergy * thisEnergy;
    }
     
    // Loop over clean clusters to identify possible vertex positions    
    for (LArPointingClusterVertexList::const_iterator iter = pointingClusterVertexCandidateList.begin(), iterEnd = pointingClusterVertexCandidateList.end(); iter != iterEnd; ++iter)
    {
        const LArPointingCluster::Vertex &clusterVertex = *iter;

        // Calculate energy of this cluster
        const Cluster* pCluster = clusterVertex.GetCluster();
        const float thisEnergy = GetEnergy( pCluster );
        const float thisEnergySquared = thisEnergy * thisEnergy;
        
	// Find associated clusters
        float primaryEnergy(0.f);
        float associatedEnergy(0.f);
        float associatedMomentumModulus(0.f);
        unsigned int primaryClusters(0);
        CartesianVector associatedMomentum(0.f,0.f,0.f);

        LArPointingClusterVertexList associatedClusterVertexList;
        RunFastReconstruction( clusterVertex.GetPosition(), clusterVertex.GetDirection(), 
                               pointingClusterMap, pointingClusterVertexCandidateList, associatedClusterVertexList, 
                               primaryClusters, primaryEnergy, associatedEnergy, associatedMomentum );
        associatedMomentumModulus = associatedMomentum.GetMagnitude();

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
            thisFigureOfMerit = (associatedEnergy*associatedMomentumModulus)/(totalEnergy*totalEnergy);
	
        outputFigureOfMeritList.insert( std::pair<const LArPointingVertex*,float>
	    ( new LArPointingVertex(clusterVertex.GetPosition(),clusterVertex.GetDirection()),thisFigureOfMerit) );
       	      
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
  
void VertexFindingAlgorithm::FindPossibleDisplacedVertices( const LArPointingClusterMap& pointingClusterMap, const LArPointingClusterVertexList& pointingClusterVertexCandidateList, VertexFigureOfMeritList& outputFigureOfMeritList )
{
    
    // Start by calculating the total energy 
    float totalEnergy(0.f);
    float totalEnergySquared(0.f);
    
    for (LArPointingClusterVertexList::const_iterator iter = pointingClusterVertexCandidateList.begin(), iterEnd = pointingClusterVertexCandidateList.end(); iter != iterEnd; ++iter)
    {
        const LArPointingCluster::Vertex &thisCluster = *iter;
        const Cluster* pCluster = thisCluster.GetCluster();
        const float thisEnergy = GetEnergy( pCluster );
        const float thisEnergySquared = thisEnergy * thisEnergy;

        totalEnergy += 0.5 * thisEnergy;
        totalEnergySquared += 0.5 * thisEnergy * thisEnergy;
    }
   
    // Loop over pairs of clean clusters to identify possible displaced vertex positions   
    for (LArPointingClusterVertexList::const_iterator iterI = pointingClusterVertexCandidateList.begin(), iterEndI = pointingClusterVertexCandidateList.end(); iterI != iterEndI; ++iterI)
    {
        const LArPointingCluster::Vertex& clusterI = *iterI;
        const Cluster* pClusterI = clusterI.GetCluster();
        const float clusterEnergyI = this->GetEnergy( pClusterI );

        for (LArPointingClusterVertexList::const_iterator iterJ = iterI, iterEndJ = pointingClusterVertexCandidateList.end(); iterJ != iterEndJ; ++iterJ)
        {
            const LArPointingCluster::Vertex& clusterJ = *iterJ;
            const Cluster* pClusterJ = clusterJ.GetCluster();
            const float clusterEnergyJ = this->GetEnergy( pClusterJ );

            if ( pClusterI == pClusterJ ) continue;

            if ( clusterEnergyI < 0.05 * totalEnergy || clusterEnergyJ < 0.05 * totalEnergy ) continue;

            bool isPhysical(false);
            CartesianVector intersectPosition(0.f,0.f,0.f);
            CartesianVector intersectDirection(0.f,0.f,0.f);

            this->GetIntersection( clusterI, clusterJ, intersectPosition, intersectDirection, isPhysical );

	    if( false==isPhysical ) continue;

// std::cout << " DISTANCE1=" << (clusterI.GetPosition()-intersectPosition).GetMagnitude() << " EMISSION1=" << IsEmission(clusterI.GetPosition(),clusterI.GetDirection(),intersectPosition) << " DISTANCE2=" << (clusterJ.GetPosition()-intersectPosition).GetMagnitude() << " EMISSION2=" << IsEmission(clusterJ.GetPosition(),clusterJ.GetDirection(),intersectPosition) << std::endl;
// ClusterList clusterList1;
// Cluster* pClusterII = (Cluster*)(pClusterI);
// Cluster* pClusterJJ = (Cluster*)(pClusterJ);
// clusterList1.insert(pClusterII);
// clusterList1.insert(pClusterJJ);
// CartesianVector vertexI = clusterI.GetPosition();
// CartesianVector vertexJ = clusterJ.GetPosition();
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&clusterList1, "CLUSTERS", GREEN);
// PandoraMonitoringApi::AddMarkerToVisualization(&vertexI,"VERTEX",BLUE, 1.);
// PandoraMonitoringApi::AddMarkerToVisualization(&vertexJ,"VERTEX",BLUE, 1.);
// PandoraMonitoringApi::AddMarkerToVisualization(&intersectPosition, "INTERSECTION", RED, 1.);
// PandoraMonitoringApi::ViewEvent();

            // Find associated clusters
            float primaryEnergy(0.f);
            float associatedEnergy(0.f);
            float associatedMomentumModulus(0.f);
            unsigned int primaryClusters(0);
            CartesianVector associatedMomentum(0.f,0.f,0.f);

            LArPointingClusterVertexList associatedClusterVertexList;
            RunFastReconstruction( intersectPosition, intersectDirection, 
                                   pointingClusterMap, pointingClusterVertexCandidateList, associatedClusterVertexList, 
                                   primaryClusters, primaryEnergy, associatedEnergy, associatedMomentum );
            associatedMomentumModulus = associatedMomentum.GetMagnitude();


            // Some quality requirements on candidate vertex
            if ( associatedClusterVertexList.empty() )
	        continue;

            if ( m_runBeamMode && this->IsConsistentWithBeamDirection(associatedMomentum) == false ) 
                continue;

            if ( primaryEnergy / totalEnergy < 0.05 * ( 5.f - static_cast<float>(primaryClusters) ) ) 
                continue; 

            // Fill the output list
            float thisFigureOfMerit = 0.f;

            static const float epsilon(0.05);

            if( totalEnergy > 0.f )
	      thisFigureOfMerit = (associatedEnergy*associatedMomentumModulus)/((1.f+epsilon)*totalEnergy*totalEnergy);

            outputFigureOfMeritList.insert( std::pair<const LArPointingVertex*,float>
	       ( new LArPointingVertex(intersectPosition,intersectDirection),thisFigureOfMerit) );
           
// std::cout << "  thisVertex: numPrimaryClusters=" << primaryClusters << " primaryEnergyFrac=" << primaryEnergy/totalEnergy << " associatedEnergyFrac=" << associatedEnergy/totalEnergy << " FigureOfMerit=" << (associatedEnergy*associatedMomentumModulus)/(totalEnergy*totalEnergy) << std::endl;
// ClusterList clusterList2;
// CollectClusters( associatedClusterVertexList, clusterList2 );
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&clusterList2, "CLUSTERS", GREEN);
// PandoraMonitoringApi::AddMarkerToVisualization(&intersectPosition, "VERTEX", RED, 1.);
// PandoraMonitoringApi::ViewEvent();

	}
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexFindingAlgorithm::CollectMarkers( const LArPointingClusterVertexList& inputList, pandora::CartesianPointList& outputList )
{
    outputList.clear();

    for (LArPointingClusterVertexList::const_iterator iter = inputList.begin(), iterEnd = inputList.end(); iter != iterEnd; ++iter )
        outputList.push_back((*iter).GetPosition());
}
 
//------------------------------------------------------------------------------------------------------------------------------------------
 
void VertexFindingAlgorithm::CollectClusters( const LArPointingClusterVertexList& inputList, pandora::ClusterList& outputList )
{
    outputList.clear();

    for (LArPointingClusterVertexList::const_iterator iter = inputList.begin(), iterEnd = inputList.end(); iter != iterEnd; ++iter )
        outputList.insert((*iter).GetCluster());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexFindingAlgorithm::GetListOfCleanVertexClusters( const LArPointingClusterVertexList& inputList, LArPointingClusterVertexList& outputList )
{
    outputList.clear();

    float totalEnergy(0.f);

    for (LArPointingClusterVertexList::const_iterator iter = inputList.begin(), iterEnd = inputList.end(); iter != iterEnd; ++iter)
    {
        const LArPointingCluster::Vertex& cluster = *iter;
        const Cluster* pCluster = cluster.GetCluster();

        totalEnergy += 0.5 * this->GetEnergy( pCluster );
    }

    for (LArPointingClusterVertexList::const_iterator iterI = inputList.begin(), iterEndI = inputList.end(); iterI != iterEndI; ++iterI)
    {
        const LArPointingCluster::Vertex& clusterI = *iterI;
        const Cluster* pClusterI = clusterI.GetCluster();

        if ( this->GetLength( pClusterI ) > 10.f 
	  || this->GetEnergy( pClusterI ) > 0.1 * totalEnergy )
        {
	    outputList.push_back(clusterI);  continue;
        }

        bool isAssociated(false);

        for (LArPointingClusterVertexList::const_iterator iterJ = inputList.begin(), iterEndJ = inputList.end(); iterJ != iterEndJ; ++iterJ)
        { 
            const LArPointingCluster::Vertex& clusterJ = *iterJ;
            const Cluster* pClusterJ = clusterJ.GetCluster();

            if ( pClusterI == pClusterJ ) continue;
	   
            if ( this->IsNode( clusterI.GetPosition(), clusterI.GetDirection(), clusterJ.GetPosition() ) 
	      || this->IsNode( clusterJ.GetPosition(), clusterJ.GetDirection(), clusterI.GetPosition() ) )
	    {
	        isAssociated = true;  break;
	    }

            if ( this->IsAdjacent( clusterI, clusterJ.GetPosition() )
	      || this->IsAdjacent( clusterJ, clusterI.GetPosition() ) )
	    {
                isAssociated = true;  break;
	    }
        }

        if ( isAssociated )
        {
            outputList.push_back(clusterI);  
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexFindingAlgorithm::RunFastReconstruction( const CartesianVector& seedVertexPosition, const CartesianVector& seedVertexDirection, const LArPointingClusterMap& inputMap, const LArPointingClusterVertexList& inputList, LArPointingClusterVertexList& outputList, unsigned int& primaryClusters, float& primaryEnergy, float& outputEnergy, CartesianVector& outputMomentum )
{
    outputList.clear();

    primaryClusters = 0;
    primaryEnergy = 0.f;
    outputEnergy = 0.f;
    outputMomentum.SetValues(0.f,0.f,0.f);

    LArPointingClusterVertexList candidateList;
    LArPointingClusterVertexList strongList;
    LArPointingClusterVertexList weakList;

    for (LArPointingClusterVertexList::const_iterator iter1 = inputList.begin(), iterEnd1 = inputList.end(); iter1 != iterEnd1; ++iter1)
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
	    candidateList.push_back(clusterVertex);
    }

    for (LArPointingClusterVertexList::const_iterator iter0 = candidateList.begin(), iterEnd0 = candidateList.end(); iter0 != iterEnd0; ++iter0)
    {
        const LArPointingCluster::Vertex &clusterVertex = *iter0;

        const float thisEnergy( this->GetEnergy( clusterVertex.GetCluster() ) );

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
    
    const float clusterEnergy = this->GetEnergy( thisCluster.GetCluster() );
    const float clusterLength = this->GetLength( thisCluster.GetCluster() );

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
 
void VertexFindingAlgorithm::GetIntersection( const LArPointingCluster::Vertex& firstCluster, const LArPointingCluster::Vertex& secondCluster, CartesianVector& intersectPosition, CartesianVector& intersectDirection, bool& isPhysical )
{
    const float cosRelativeAngle( (firstCluster.GetDirection()).GetDotProduct(secondCluster.GetDirection()) );

    const float cosAlmostParallelAngle( 0.999 ); // 2.5 degrees

    if ( std::fabs(cosRelativeAngle) < cosAlmostParallelAngle ) // good pointing
    {
        LArPointingCluster::Vertex::GetIntersection( firstCluster, secondCluster, intersectPosition, isPhysical );
	LArPointingCluster::Vertex::GetAverageDirection( firstCluster, secondCluster, intersectDirection );
    }
    else if ( cosRelativeAngle < -cosAlmostParallelAngle ) // anti-parallel (bad)
    {
        intersectPosition = ( firstCluster.GetPosition() + secondCluster.GetPosition() ) * 0.5;
        intersectDirection = firstCluster.GetDirection(); // don't care, pick a safe value
        isPhysical = false;
    }
    else if ( cosRelativeAngle > +cosAlmostParallelAngle ) // parallel (good)
    {
        const CartesianVector averagePosition = (firstCluster.GetPosition() + secondCluster.GetPosition()) * 0.5;
        const CartesianVector averageDirection = (firstCluster.GetDirection() + secondCluster.GetDirection()).GetUnitVector();
        const float Width = averageDirection.GetCrossProduct(secondCluster.GetPosition() - firstCluster.GetPosition()).GetMagnitude(); 
        const float deltaLength = std::max( averageDirection.GetDotProduct(averagePosition - firstCluster.GetPosition()),
                                            averageDirection.GetDotProduct(averagePosition - secondCluster.GetPosition()) );
        intersectPosition = averagePosition - averageDirection * (20.f * Width + deltaLength + 2.5);
        intersectDirection = averageDirection;
        isPhysical = true;

// std::cout << "ALMOST PARALLEL" << std::endl;
// ClusterList clusterList;
// Cluster* pCluster1 = firstCluster.GetCluster();
// Cluster* pCluster2 = secondCluster.GetCluster();
// clusterList.insert(pCluster1);
// clusterList.insert(pCluster2);
// CartesianVector vertex1 = firstCluster.GetPosition();
// CartesianVector vertex2 = secondCluster.GetPosition();
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&clusterList, "CLUSTERS", GREEN);
// PandoraMonitoringApi::AddMarkerToVisualization(&vertex1,"VERTEX",BLUE, 1.);
// PandoraMonitoringApi::AddMarkerToVisualization(&vertex2,"VERTEX",BLUE, 1.);
// PandoraMonitoringApi::AddMarkerToVisualization(&intersectPosition, "INTERSECTION", RED, 1.);
// PandoraMonitoringApi::ViewEvent();

    }
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

    static const float tanTheta( tan( M_PI * 10.0 / 180.0 ) );

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

   const float daughterLength = this->GetLength( daughterCluster.GetCluster() );

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

    if ( rL > 0.f && rL < 125.f && rT*rT < 2.5*2.5 + rL*rL*tanSqTheta ) return true;  else return false;
}
 
//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexFindingAlgorithm::IsConsistentWithBeamDirection( const CartesianVector& thisMomentum ) const
{
    if( thisMomentum.GetUnitVector().GetZ() > -0.25 ) return true;  else return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexFindingAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector)
{
    for ( ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter ) 
    {
        Cluster *pCluster = *iter;

        if ( 1 + pCluster->GetOuterPseudoLayer() - pCluster->GetInnerPseudoLayer() < 10 ) continue;

        const CartesianVector innerVertex = pCluster->GetCentroid( pCluster->GetInnerPseudoLayer() );
        const CartesianVector outerVertex = pCluster->GetCentroid( pCluster->GetOuterPseudoLayer() );

        if ( (outerVertex-innerVertex).GetMagnitudeSquared() < 15.f ) continue;

        if ( LArParticleId::LArLayerOccupancy(pCluster) < 0.75 ) continue;

        clusterVector.push_back(pCluster);
    }
}
//------------------------------------------------------------------------------------------------------------------------------------------

float VertexFindingAlgorithm::GetEnergy( const Cluster* const pCluster ) const
{
    static const float dEdX( 0.002 ); // approximately 2 MeV/cm

    return dEdX * GetLength( pCluster );
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VertexFindingAlgorithm::GetLength( const Cluster* const pCluster ) const
{
    const CartesianVector innerCentroid(pCluster->GetCentroid(pCluster->GetInnerPseudoLayer()));
    const CartesianVector outerCentroid(pCluster->GetCentroid(pCluster->GetOuterPseudoLayer()));

    return (outerCentroid-innerCentroid).GetMagnitude();
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
