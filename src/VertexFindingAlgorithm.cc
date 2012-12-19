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

    const ClusterList* pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentClusterList(*this, pClusterList));

    CartesianVector recoVertex(0.f,0.f,0.f);

    // Set a first pass vertex
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, GetFirstPassVertex(pClusterList,recoVertex)); 


// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(pClusterList, "CLUSTERS", GREEN);
// PandoraMonitoringApi::AddMarkerToVisualization(&recoVertex, "RECO", RED, 1.);
// PandoraMonitoringApi::ViewEvent();

    
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


// ClusterList inputClusterVertexList, outputClusterVertexList;
// CollectClusters( pointingClusterVertexList, inputClusterVertexList );
// CollectClusters( pointingClusterVertexCandidateList, outputClusterVertexList );
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&inputClusterVertexList, "CLEAN CLUSTERS, FIRST PASS", GREEN);
// PandoraMonitoringApi::VisualizeClusters(&outputClusterVertexList, "CLEAN CLUSTERS, SECOND PASS", BLUE);
// PandoraMonitoringApi::ViewEvent();


    // The figure of merit
    float figureOfMerit(0.f);
    
    static const float epsilon(0.05);

    // Find the best connected vertex
    float bestConnectedEnergy(0.f);
    float bestConnectedFigureOfMerit(0.f);
    CartesianVector bestConnectedVertex(0.f,0.f,0.f);
    CartesianVector bestConnectedMomentum(0.f,0.f,0.f);
    LArPointingClusterVertexList bestConnectedClusterVertexList;
 
    this->FindBestConnectedVertex( pointingClusterMap, pointingClusterVertexCandidateList, bestConnectedClusterVertexList,
                                   bestConnectedVertex, bestConnectedEnergy, bestConnectedMomentum, bestConnectedFigureOfMerit );

    if ( bestConnectedFigureOfMerit > figureOfMerit )
    {
        figureOfMerit = bestConnectedFigureOfMerit;
        recoVertex = bestConnectedVertex;   
    }


    // Find the best displaced vertex
    float bestDisplacedEnergy(0.f);
    float bestDisplacedFigureOfMerit(0.f);
    CartesianVector bestDisplacedVertex(0.f,0.f,0.f);
    CartesianVector bestDisplacedMomentum(0.f,0.f,0.f);
    LArPointingClusterVertexList bestDisplacedClusterVertexList;

    this->FindBestDisplacedVertex( pointingClusterMap, pointingClusterVertexCandidateList, bestDisplacedClusterVertexList,
                                   bestDisplacedVertex, bestDisplacedEnergy, bestDisplacedMomentum, bestDisplacedFigureOfMerit );

    if ( bestDisplacedFigureOfMerit > (1.f + epsilon) * figureOfMerit )
    {
        figureOfMerit = bestDisplacedFigureOfMerit;
        recoVertex = bestDisplacedVertex;   
    }


    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, SetVertex(recoVertex));


// ClusterList cleanClusterList;
// ClusterList cleanVertexClusterList;
// for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
// {
// Cluster* pCluster = *iter;
// if ( LArVertexHelper::IsConnectedToCurrentVertex(pCluster) ) cleanVertexClusterList.insert(pCluster);
// else cleanClusterList.insert(pCluster);
// }
// const CartesianVector &theVertex(LArVertexHelper::GetCurrentVertex());
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&cleanClusterList, "CLUSTERS", GREEN);
// PandoraMonitoringApi::VisualizeClusters(&cleanVertexClusterList, "VERTEX CLUSTERS", BLUE);
// PandoraMonitoringApi::AddMarkerToVisualization(&theVertex, "RECO", RED, 1.);
// PandoraMonitoringApi::AddMarkerToVisualization(&trueVertex, "TRUE", BLUE, 1.);
// PandoraMonitoringApi::ViewEvent();


    return STATUS_CODE_SUCCESS;
}

StatusCode VertexFindingAlgorithm::GetFirstPassVertex( const ClusterList* const pClusterList, CartesianVector& firstVertex )
{
    if ( pClusterList->empty() )
        return STATUS_CODE_NOT_INITIALIZED;

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

    return STATUS_CODE_SUCCESS;
}

void VertexFindingAlgorithm::FindBestConnectedVertex(const LArPointingClusterMap& pointingClusterMap, const LArPointingClusterVertexList& pointingClusterVertexCandidateList, LArPointingClusterVertexList& outputClusterVertexList, CartesianVector& outputVertex, float& outputEnergy, CartesianVector& outputMomentum, float& outputFigureOfMerit )
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
     
    // Loop over clean clusters to identify the best candidate vertex position    
    float bestEnergy(0.f);
    float bestMomentumModulus(0.f);
    CartesianVector bestMomentum(0.f,0.f,0.f);
    CartesianVector bestVertex(0.f,0.f,0.f);
    LArPointingClusterVertexList bestClusterVertexList;

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

       
        // Best candidate vertex
        if ( associatedEnergy*associatedMomentumModulus > bestEnergy*bestMomentumModulus )
	{
            bestVertex = clusterVertex.GetPosition();
            bestEnergy = associatedEnergy;
            bestMomentum = associatedMomentum;
            bestMomentumModulus = associatedMomentumModulus;
            bestClusterVertexList = associatedClusterVertexList;
	}
	      
// std::cout << "  thisVertex: energyFrac=" << thisEnergy/totalEnergy << " energySquaredFrac=" << thisEnergySquared/totalEnergySquared << " numPrimaryClusters=" << primaryClusters << " primaryEnergyFrac=" << primaryEnergy/totalEnergy << " associatedEnergyFrac=" << associatedEnergy/totalEnergy << " FigureOfMerit=" << (associatedEnergy*associatedMomentumModulus)/(totalEnergy*totalEnergy) << std::endl;
// CartesianVector thisVertex = clusterVertex.GetPosition();
// ClusterList clusterList;
// CollectClusters( associatedClusterVertexList, clusterList );
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&clusterList, "CLUSTERS", GREEN);
// PandoraMonitoringApi::AddMarkerToVisualization(&thisVertex, "VERTEX", RED, 1.);
// PandoraMonitoringApi::ViewEvent();

    }

    outputVertex = bestVertex;
    outputEnergy = bestEnergy;
    outputMomentum = bestMomentum;
    outputClusterVertexList = bestClusterVertexList;

    if( totalEnergy > 0.f )
        outputFigureOfMerit = (bestEnergy*bestMomentumModulus)/(totalEnergy*totalEnergy);
    else
        outputFigureOfMerit = 0.f;

// std::cout << " bestVertex: bestMomentumFrac=" << bestMomentumModulus/totalEnergy << " bestEnergyFrac=" << bestEnergy/totalEnergy << " FigureOfMerit=" << (bestMomentumModulus*bestEnergy) / (totalEnergy*totalEnergy) << std::endl;
// ClusterList clusterList;
// CartesianPointList pointList;
// CollectClusters( bestClusterVertexList, clusterList );
// CollectMarkers( bestClusterVertexList, pointList );
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&clusterList, "CLUSTERS", GREEN);
// for( unsigned int n=0; n<pointList.size(); ++n ){
//    CartesianVector myPoint = pointList.at(n);
//    PandoraMonitoringApi::AddMarkerToVisualization(&myPoint, "CLUSTERS", BLUE, 1.);
// }
// PandoraMonitoringApi::AddMarkerToVisualization(&bestVertex, "VERTEX", RED, 1.);
// PandoraMonitoringApi::ViewEvent();

}

//------------------------------------------------------------------------------------------------------------------------------------------
  
void VertexFindingAlgorithm::FindBestDisplacedVertex( const LArPointingClusterMap& pointingClusterMap, const LArPointingClusterVertexList& pointingClusterVertexCandidateList, LArPointingClusterVertexList& outputClusterVertexList, CartesianVector& outputVertex, float& outputEnergy, CartesianVector& outputMomentum, float& outputFigureOfMerit )
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
   
    // Loop over pairs of clean clusters to identify the best displaced vertex position    
    float bestEnergy(0.f);
    float bestMomentumModulus(0.f);
    CartesianVector bestMomentum(0.f,0.f,0.f);
    CartesianVector bestVertex(0.f,0.f,0.f);
    LArPointingClusterVertexList bestClusterVertexList;

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

       
            // Best candidate vertex
            if ( associatedEnergy*associatedMomentumModulus > bestEnergy*bestMomentumModulus )
	    {
                bestVertex = intersectPosition;
                bestEnergy = associatedEnergy;
                bestMomentum = associatedMomentum;
                bestMomentumModulus = associatedMomentumModulus;
                bestClusterVertexList = associatedClusterVertexList;
	    }

// std::cout << "  thisVertex: numPrimaryClusters=" << primaryClusters << " primaryEnergyFrac=" << primaryEnergy/totalEnergy << " associatedEnergyFrac=" << associatedEnergy/totalEnergy << " FigureOfMerit=" << (associatedEnergy*associatedMomentumModulus)/(totalEnergy*totalEnergy) << std::endl;
// ClusterList clusterList2;
// CollectClusters( associatedClusterVertexList, clusterList2 );
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&clusterList2, "CLUSTERS", GREEN);
// PandoraMonitoringApi::AddMarkerToVisualization(&intersectPosition, "VERTEX", RED, 1.);
// PandoraMonitoringApi::ViewEvent();

	}
    }

  
    outputVertex = bestVertex;
    outputEnergy = bestEnergy;
    outputMomentum = bestMomentum;
    outputClusterVertexList = bestClusterVertexList;

    if( totalEnergy > 0.f )
        outputFigureOfMerit = (bestEnergy*bestMomentumModulus)/(totalEnergy*totalEnergy);
    else
        outputFigureOfMerit = 0.f;

// std::cout << " bestVertex: bestMomentumFrac=" << bestMomentumModulus/totalEnergy << " bestEnergyFrac=" << bestEnergy/totalEnergy << " FigureOfMerit=" << (bestMomentumModulus*bestEnergy) / (totalEnergy*totalEnergy) << std::endl;
// ClusterList clusterList;
// CartesianPointList pointList;
// CollectClusters( bestClusterVertexList, clusterList );
// CollectMarkers( bestClusterVertexList, pointList );
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&clusterList, "CLUSTERS", GREEN);
// for( unsigned int n=0; n<pointList.size(); ++n ){
//    CartesianVector myPoint = pointList.at(n);
//    PandoraMonitoringApi::AddMarkerToVisualization(&myPoint, "CLUSTERS", BLUE, 1.);
// }
// PandoraMonitoringApi::AddMarkerToVisualization(&bestVertex, "VERTEX", RED, 1.);
// PandoraMonitoringApi::ViewEvent();   

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
    CartesianVector trueVertex(130.5f, 0.f, 100.f); // 128.2f,0.f,100.f 

    LArVertexHelper::AddVertex(m_vertexName, trueVertex);
    LArVertexHelper::SetCurrentVertex(m_vertexName);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexFindingAlgorithm::SetVertex(const CartesianVector& eventVertex)
{
    LArVertexHelper::AddVertex(m_vertexName, eventVertex);
    LArVertexHelper::SetCurrentVertex(m_vertexName);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexFindingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "VertexName", m_vertexName));

    m_useTrueVertex = false;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "UseTrueVertex", m_useTrueVertex));

    m_runBeamMode = true; 
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "RunBeamMode", m_runBeamMode));


    return STATUS_CODE_SUCCESS;
}


///////////////////////////////////////////////////////////////////////////////////
// 

StatusCode VertexFindingAlgorithm::RunOldMethod(const ClusterList *const pClusterList)
{
  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PrepareData(pClusterList));
  
  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, FindVertex());

  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexFindingAlgorithm::PrepareData(const ClusterList *const pClusterList)
{
  nClusters = 0;
  nClusterHits = 0;

  ClusterVector clusterVector;
  GetListOfCleanClusters(pClusterList, clusterVector);
  std::sort(clusterVector.begin(), clusterVector.end(), Cluster::SortByInnerLayer);

  unsigned int ctr0 = 0;

  for( ClusterVector::const_iterator iter = clusterVector.begin(); iter != clusterVector.end(); ++iter ) {
    Cluster* pCluster = *iter;
    OrderedCaloHitList clusterHitList = pCluster->GetOrderedCaloHitList();

    PseudoLayer innerLayer = pCluster->GetInnerPseudoLayer();
    PseudoLayer outerLayer = pCluster->GetOuterPseudoLayer();

    ClusterHelper::ClusterFitResult innerLayerFit;
    ClusterHelper::ClusterFitResult outerLayerFit;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, ClusterHelper::FitStart(pCluster, 50, innerLayerFit));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, ClusterHelper::FitEnd(pCluster, 50, outerLayerFit));

    MyCluster* newCluster = NewCluster();
    (*newCluster).icluster = ctr0;
    (*newCluster).vx       = pCluster->GetCentroid(innerLayer).GetX();
    (*newCluster).vz       = pCluster->GetCentroid(innerLayer).GetZ();
    (*newCluster).ex       = pCluster->GetCentroid(outerLayer).GetX();
    (*newCluster).ez       = pCluster->GetCentroid(outerLayer).GetZ();
    (*newCluster).px       = innerLayerFit.GetDirection().GetX();
    (*newCluster).pz       = innerLayerFit.GetDirection().GetZ();
    (*newCluster).qx       = outerLayerFit.GetDirection().GetX();
    (*newCluster).qz       = outerLayerFit.GetDirection().GetZ();
    (*newCluster).n        = 1+outerLayer-innerLayer;

    for ( OrderedCaloHitList::const_iterator clusterHitIter = clusterHitList.begin(); clusterHitIter != clusterHitList.end(); ++clusterHitIter ) {
      PseudoLayer  layer         = clusterHitIter->first;
      CaloHitList* layerHitList  = clusterHitIter->second;

      unsigned int ctr1 = 0;

      for( CaloHitList::const_iterator layerHitIter = layerHitList->begin(); layerHitIter != layerHitList->end(); ++layerHitIter ) {
        CaloHit* pCaloHit = *layerHitIter;

        MyClusterHit* newHit = NewClusterHit();
        (*newHit).icluster = ctr0;
        (*newHit).ilayer   = layer;
        (*newHit).ihit     = ctr1;
        (*newHit).x        = pCaloHit->GetPositionVector().GetX();
        (*newHit).z        = pCaloHit->GetPositionVector().GetZ();
        (*newHit).n1       = 1+outerLayer-layer; 
        (*newHit).n2       = 1+layer-innerLayer; 
      }
      ++ctr1;
    }
    ++ctr0;
  }

  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexFindingAlgorithm::FindVertex()
{
  CartesianVector eventVertex(CalcNearestHitToVertex());

  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, SetVertex(eventVertex));

  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector VertexFindingAlgorithm::CalcNearestHitToVertex()
{
  // Fit Parameters
  double vtx_x = 0.0;
  double vtx_z = 0.0;

  double lnl_best = 0.0;
  double lnl = 0.0;

  double x = 0.0;
  double z = 0.0;

  // Loop over all hits
  for( unsigned int i=0; i<nClusterHits; ++i ){
    MyClusterHit* hit = (MyClusterHit*)(vClusterHits.at(i));

    x = (*hit).x;
    z = (*hit).z;
    lnl = CalcLikelihoodFunction( x, z );

    if( lnl<lnl_best ){
      vtx_x = x;
      vtx_z = z;
      lnl_best = lnl;
    }
  }

  return CartesianVector(vtx_x, 0., vtx_z);
}

//------------------------------------------------------------------------------------------------------------------------------------------

double VertexFindingAlgorithm::CalcLikelihoodFunction( double vtxX, double vtxZ )
{
  return CalcLikelihoodFunctionMethod1( vtxX, vtxZ );
}

//------------------------------------------------------------------------------------------------------------------------------------------

double VertexFindingAlgorithm::CalcLikelihoodFunctionMethod1( double vtxX, double vtxZ )
{
  double f  = 0.0;

  double vx  = 0.0;
  double vz  = 0.0;
  double ex  = 0.0;
  double ez  = 0.0;
  double px = 0.0;
  double pz = 0.0;
  double qx = 0.0;
  double qz = 0.0;
  double n  = 0;

  double dx = 0.0;
  double dz = 0.0;
  double dr = 0.0;

  double costh = 0.0;
  double pcosth = 0.0;
  double qcosth = 0.0;

  double dz1 = 0.0;
  double dz2 = 0.0;

  for( unsigned i=0; i<nClusters; ++i ){
    MyCluster* cluster = (MyCluster*)(vClusters.at(i));
  
    vx = (*cluster).vx;
    vz = (*cluster).vz;
    ex = (*cluster).ex;
    ez = (*cluster).ez;
    px = (*cluster).px;
    pz = (*cluster).pz;
    qx = (*cluster).qx;
    qz = (*cluster).qz;
    n  = (*cluster).n;

    dx = (vx-vtxX);
    dz = (vz-vtxZ);
    dr = sqrt(dx*dx+dz*dz);
    if( dr>0.0 ) pcosth = +(dx*px+dz*pz)/dr;
    else         pcosth = 1.0;  
      
    dx = (ex-vtxX);
    dz = (ez-vtxZ);
    dr = sqrt(dx*dx+dz*dz);
    if( dr>0.0 ) qcosth = -(dx*px+dz*pz)/dr;
    else         qcosth = 1.0;  

    costh = std::max(pcosth,qcosth);

    dz1 = ez-vz;
    dz2 = ez-vz;
    if( vtxZ>vz && vtxZ<ez ){
      dz1 = +(vtxZ-vz);
      dz2 = -(vtxZ-ez);
    }

    f += dz1*dz2*costh;
  }

  return -f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

double VertexFindingAlgorithm::CalcLikelihoodFunctionMethod2( double vtxX, double vtxZ )
{
  double f  = 0.0;

  double vx = 0.0;
  double vz = 0.0;
  double et = 0.0;
  double ex = 0.0;
  double ez = 0.0;
  double px = 0.0;
  double pz = 0.0;
  double qx = 0.0;
  double qz = 0.0;
  double n  = 0;

  double dx = 0.0;
  double dz = 0.0;

  double vT = 0.0;
  double vL = 0.0;
  double eT = 0.0;
  double eL = 0.0;
  double dT = 0.0;
  double dL = 0.0;

  double L2 = 0.0;
  double X2 = 0.0;
  double T2 = 2.0*2.0;
  
  double k = 2.0;
  double x = 0.0;

  for( unsigned i=0; i<nClusters; ++i ){
    MyCluster* cluster = (MyCluster*)(vClusters.at(i));
  
    vx = (*cluster).vx;
    vz = (*cluster).vz;
    ex = (*cluster).ex;
    ez = (*cluster).ez;
    px = (*cluster).px;
    pz = (*cluster).pz;
    qx = (*cluster).qx;
    qz = (*cluster).qz;
    n  = (*cluster).n;

    vL = (vx-vtxX)*px + (vz-vtxZ)*pz;
    vT = (vx-vtxX)*pz - (vz-vtxZ)*px;

    eL = (ex-vtxX)*qx + (ez-vtxZ)*qz;
    eT = (ex-vtxX)*qz - (ez-vtxZ)*qx;

    L2 = (ex-vx)*(ex-vx) + (ez-vz)*(ez-vz); 

    // project back from start of cluster
         if( vtxZ<vz ) x = 0.0;
    else if( vtxZ<ez ) x = (vtxZ-vz)/(ez-vz);
    else               x = 1.0;

    X2 = T2*(1.0+k*vL*vL/L2);

    if( 1.0-(vT*vT)/X2>0.0 ){
      f += 0.573*L2*( 1.0/(1.0+7.0*x*x)-0.125 )
                   *( 1.0-(vT*vT)/X2 )/sqrt(X2);
    }

    // project forward from end of cluster
         if( vtxZ>ez ) x = 0.0;
    else if( vtxZ>vz ) x = (ez-vtxZ)/(ez-vz);
    else               x = 1.0;

    X2 = T2*(1.0+k*eL*eL/L2);

    if( 1.0-(eT*eT)/X2>0.0 ){
      f += 0.573*L2*( 1.0/(1.0+7.0*(x-1.0)*(x-1.0))-0.125 )
                   *( 1.0-(eT*eT)/X2 )/sqrt(X2);
    }

  }

  return -f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

VertexFindingAlgorithm::MyCluster* VertexFindingAlgorithm::NewCluster()
{
  MyCluster* cluster = NULL;

  if( nClusters<vClusters.size() ){
    cluster = (MyCluster*)(vClusters.at(nClusters));                 
  }
  else{
    cluster = new MyCluster();
    vClusters.push_back(cluster);
  }

  (*cluster).icluster = -1;
  (*cluster).vx       = 0.0;
  (*cluster).vz       = 0.0;
  (*cluster).ex       = 0.0;
  (*cluster).ez       = 0.0;
  (*cluster).px       = 0.0;
  (*cluster).pz       = 0.0;
  (*cluster).qx       = 0.0;
  (*cluster).qz       = 0.0;
  (*cluster).n        = 0;

  ++nClusters;

  return cluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
VertexFindingAlgorithm::MyClusterHit* VertexFindingAlgorithm::NewClusterHit()
{
  MyClusterHit* hit = NULL;

  if( nClusterHits<vClusterHits.size() ){
    hit = (MyClusterHit*)(vClusterHits.at(nClusterHits));                 
  }
  else{
    hit = new MyClusterHit();
    vClusterHits.push_back(hit);
  }

  (*hit).icluster = -1;
  (*hit).ilayer   = -1;
  (*hit).ihit     = -1;
  (*hit).x        = 0.0;
  (*hit).z        = 0.0;
  (*hit).n1       = 0; 
  (*hit).n2       = 0; 

  ++nClusterHits;

  return hit;
}

//
///////////////////////////////////////////////////////////////////////////////////


} // namespace lar
