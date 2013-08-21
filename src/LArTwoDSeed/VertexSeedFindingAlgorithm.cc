/**
 *  @file   LArContent/src/LArTwoDSeed/VertexSeedFindingAlgorithm.cc
 * 
 *  @brief  Implementation of the vertex seed finding algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArPointingClusterHelper.h"
#include "LArHelpers/LArVertexHelper.h"

#include "LArTwoDSeed/VertexSeedFindingAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode VertexSeedFindingAlgorithm::Run()
{
    // Get current cluster list
    const ClusterList *pInputClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentClusterList(*this, pInputClusterList));


    // Get list of clean clusters
    ClusterVector clusterVector;
    this->GetListOfCleanClusters(pInputClusterList, clusterVector);
    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNOccupiedLayers);


    // Identify clusters associated with vertex 
    ClusterList vertexSeedClusterList;
    this->GetListOfVertexClusters(clusterVector,vertexSeedClusterList);
    

    // Cluster list management
    ClusterList nonSeedClusterList(*pInputClusterList);

    for (ClusterList::const_iterator iter = vertexSeedClusterList.begin(), iterEnd = vertexSeedClusterList.end(); iter != iterEnd; ++iter)
    {
        nonSeedClusterList.erase(*iter);
    }

    if (!vertexSeedClusterList.empty())
    {
        if (!nonSeedClusterList.empty())
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveClusterList(*this, m_nonSeedClusterListName, nonSeedClusterList));
        }

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveClusterList(*this, m_seedClusterListName, vertexSeedClusterList));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentClusterList(*this, m_seedClusterListName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSeedFindingAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (LArClusterHelper::GetLayerSpan(pCluster) < m_minClusterLayers)
	    continue;

        if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLengthSquared)
	    continue;

        clusterVector.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSeedFindingAlgorithm::GetListOfVertexClusters(const ClusterVector &clusterVector, ClusterList& seedClusterList) const
{
    // Get the current vertex
    const CartesianVector eventVertex(LArVertexHelper::GetCurrentVertex());


    // Generate a map of pointing clusters 
    LArPointingClusterMap pointingClusterMap;

    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
        pointingClusterMap.insert( std::pair<Cluster*,LArPointingCluster>(*iter,LArPointingCluster(*iter)) );


    // Identify nodes and emissions
    LArPointingClusterVertexList emissions, associations, seeds;

    for (LArPointingClusterMap::const_iterator iter = pointingClusterMap.begin(), iterEnd = pointingClusterMap.end(); iter != iterEnd; ++iter)
    {
        const LArPointingCluster& pointingCluster = iter->second;
        const LArPointingCluster::Vertex& innerVertex = pointingCluster.GetInnerVertex();
        const LArPointingCluster::Vertex& outerVertex = pointingCluster.GetOuterVertex();

        const bool possibleNode( pointingCluster.GetCluster()->GetNCaloHits() >= m_minClusterHitsNode );
        const bool possibleEmission( pointingCluster.GetCluster()->GetNCaloHits() >= m_minClusterHitsEmission );

        const float innerDisplacementSquared( (eventVertex - innerVertex.GetPosition()).GetMagnitudeSquared() );
        const float outerDisplacementSquared( (eventVertex - outerVertex.GetPosition()).GetMagnitudeSquared() );

        if (innerDisplacementSquared < outerDisplacementSquared)
	{
	    if (LArPointingClusterHelper::IsNode(eventVertex,innerVertex))
	    {
	        if (possibleNode) seeds.push_back(innerVertex);
	    }
            else if(LArPointingClusterHelper::IsEmission(eventVertex,innerVertex))
	    {
	        if (possibleEmission) emissions.push_back(innerVertex);
	    }
	}
        else
	{
            if (LArPointingClusterHelper::IsNode(eventVertex,outerVertex))
	    {
	        if (possibleNode) seeds.push_back(outerVertex);
	    }
            else if(LArPointingClusterHelper::IsEmission(eventVertex,outerVertex))
	    {
	        if (possibleEmission) emissions.push_back(outerVertex);
	    }
	} 
    }

    associations.insert(associations.end(), seeds.begin(), seeds.end());
    associations.insert(associations.end(), emissions.begin(), emissions.end());

    // Sort through the list of emissions
    for (LArPointingClusterVertexList::const_iterator iterI = emissions.begin(), iterEndI = emissions.end(); iterI != iterEndI; ++iterI)
    {
        const LArPointingCluster::Vertex& vertexI = *iterI;

        bool isSeed(true);

        for (LArPointingClusterVertexList::const_iterator iterJ = associations.begin(), iterEndJ = associations.end(); iterJ != iterEndJ; ++iterJ)
        {
            const LArPointingCluster::Vertex& vertexJ = *iterJ;

            if (vertexI.GetCluster() == vertexJ.GetCluster())
	        continue;

            // Check proximity between vertexI and vertexJ
            if (2 * vertexI.GetCluster()->GetOrderedCaloHitList().size() < vertexJ.GetCluster()->GetOrderedCaloHitList().size() &&
                (vertexI.GetPosition() - vertexJ.GetPosition()).GetMagnitudeSquared() > 5.f * 5.f &&
                LArClusterHelper::GetClosestDistance(vertexI.GetPosition(),vertexJ.GetCluster()) < 2.5 )
	    {
                isSeed = false;
                break;
	    }

            // Check proximity between vertexI and endJ
            LArPointingClusterMap::const_iterator lookupJ = pointingClusterMap.find(vertexJ.GetCluster());

            if ( lookupJ == pointingClusterMap.end() ) 
                throw pandora::StatusCodeException(STATUS_CODE_FAILURE);

            const LArPointingCluster& clusterJ = lookupJ->second;
            const LArPointingCluster::Vertex& endJ = vertexJ.IsInner() ? clusterJ.GetOuterVertex() : clusterJ.GetInnerVertex();
        
            if ( LArPointingClusterHelper::IsNode(endJ.GetPosition(),vertexI) ||
                (LArPointingClusterHelper::IsEmission(endJ.GetPosition(),vertexI) &&
                 endJ.GetDirection().GetDotProduct(vertexI.GetDirection()) < -0.707 ) )
	    {
	        isSeed = false;
                break;
	    }
	}

        if (isSeed) 
            seeds.push_back(vertexI);
    }

    // Populate list of vertex clusters
    for (LArPointingClusterVertexList::const_iterator iter = seeds.begin(), iterEnd = seeds.end(); iter != iterEnd; ++iter)
    {
        const LArPointingCluster::Vertex& vertexCluster = *iter;
        seedClusterList.insert(vertexCluster.GetCluster());
    }

// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&seedClusterList, "SeedClusters", AUTOITER);
// PandoraMonitoringApi::AddMarkerToVisualization(&eventVertex, "Vertex", BLACK, 2.5);
// PandoraMonitoringApi::ViewEvent();

}

//------------------------------------------------------------------------------------------------------------------------------------------

/*
// KEEP OLD METHOD FOR NOW
void VertexSeedFindingAlgorithm::GetListOfVertexClusters(const ClusterVector &clusterVector, ClusterList &vertexSeedClusterList) const
{ 
    // Get the current vertex
    const CartesianVector eventVertex(LArVertexHelper::GetCurrentVertex());


    // Create list of clusters associated with the vertex
    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
    {
        LArPointingCluster pointingCluster(*iter);

        if (LArPointingClusterHelper::IsNode(eventVertex, pointingCluster.GetInnerVertex().GetPosition()) ||
            LArPointingClusterHelper::IsNode(eventVertex, pointingCluster.GetOuterVertex().GetPosition()) ||
            LArPointingClusterHelper::IsEmitted(eventVertex, pointingCluster.GetInnerVertex()) ||
            LArPointingClusterHelper::IsEmitted(eventVertex, pointingCluster.GetOuterVertex()))
        {
            vertexSeedClusterList.insert(*iter);
        }
    }


    // Merge or delete clusters
    for (ClusterList::iterator iterI = vertexSeedClusterList.begin(); iterI != vertexSeedClusterList.end(); ++iterI)
    {
        Cluster *pClusterI = *iterI;

        for (ClusterList::iterator iterJ = vertexSeedClusterList.begin(); iterJ != vertexSeedClusterList.end(); )
        {
            Cluster *pClusterJ = *iterJ;

            if (pClusterI == pClusterJ)
            {
                ++iterJ;
                continue;
            }

            // TODO, improve efficiency - currently need to recalculate every time
            LArPointingCluster pointingClusterI(*iterI);
            const bool isNodeInnerI(LArPointingClusterHelper::IsNode(eventVertex, pointingClusterI.GetInnerVertex().GetPosition()));
            const bool isNodeOuterI(LArPointingClusterHelper::IsNode(eventVertex, pointingClusterI.GetOuterVertex().GetPosition()));

            LArPointingCluster pointingClusterJ(*iterJ);
            const bool isNodeInnerJ(LArPointingClusterHelper::IsNode(eventVertex, pointingClusterJ.GetInnerVertex().GetPosition()));
            const bool isNodeOuterJ(LArPointingClusterHelper::IsNode(eventVertex, pointingClusterJ.GetOuterVertex().GetPosition()));
            const bool isEmittedInnerJ(LArPointingClusterHelper::IsEmitted(eventVertex, pointingClusterJ.GetInnerVertex()));
            const bool isEmittedOuterJ(LArPointingClusterHelper::IsEmitted(eventVertex, pointingClusterJ.GetOuterVertex()));

            // Whether to remove cluster J from list of vertex seeds
            bool removeDaughter(false);

            if (isNodeInnerI && !isNodeInnerJ && isEmittedInnerJ &&
                (2 * pClusterJ->GetOrderedCaloHitList().size() < pClusterI->GetOrderedCaloHitList().size()) &&
                ((pointingClusterJ.GetInnerVertex().GetPosition() - pointingClusterI.GetInnerVertex().GetPosition()).GetMagnitudeSquared() > 5.f * 5.f) &&
                (LArClusterHelper::GetClosestDistance(pointingClusterJ.GetInnerVertex().GetPosition(), pClusterI) < 2.f))
            {
                removeDaughter = true;
            }
            else if (isNodeOuterI && !isNodeOuterJ && isEmittedOuterJ &&
                (2 * pClusterJ->GetOrderedCaloHitList().size() < pClusterI->GetOrderedCaloHitList().size()) &&
                ((pointingClusterJ.GetOuterVertex().GetPosition() - pointingClusterI.GetOuterVertex().GetPosition()).GetMagnitudeSquared() > 5.f * 5.f) &&
                (LArClusterHelper::GetClosestDistance(pointingClusterJ.GetOuterVertex().GetPosition(), pClusterI) < 2.f))
            {
                removeDaughter = true;
            }

            // Whether to merge clusters I and J
            bool mergeDaughter(false);

            if (!isNodeOuterI && !isNodeInnerJ && isEmittedInnerJ &&
                ((pointingClusterI.GetOuterVertex().GetPosition() - pointingClusterJ.GetInnerVertex().GetPosition()).GetMagnitudeSquared() < 2.f * 2.f))
            {
                mergeDaughter = true;
            }
            else if (!isNodeInnerI && !isNodeOuterJ && isEmittedOuterJ &&
                ((pointingClusterI.GetInnerVertex().GetPosition() - pointingClusterJ.GetOuterVertex().GetPosition()).GetMagnitudeSquared() < 2.f * 2.f))
            {
                mergeDaughter = true;
            }

// if (removeDaughter || mergeDaughter)
// {
//    if (removeDaughter) std::cout << "Remove daughter " << std::endl;
//    if (mergeDaughter) std::cout << "Merge daughter " << std::endl;
//    if (removeDaughter && mergeDaughter) std::cout << "Both merge and remove daughter - PROBLEM! " << std::endl;
//    ClusterList parent, daughter; parent.insert(pClusterI); daughter.insert(pClusterJ);
//    PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
//    PandoraMonitoringApi::VisualizeClusters(&vertexSeedClusterList, "vertexSeedClusterList", BLUE);
//    PandoraMonitoringApi::VisualizeClusters(&parent, "parent", RED);
//    PandoraMonitoringApi::VisualizeClusters(&daughter, "daughter", GREEN);
//    PandoraMonitoringApi::ViewEvent();
// }
            if (removeDaughter)
            {
                vertexSeedClusterList.erase(iterJ);
                iterJ = vertexSeedClusterList.begin();
            }
            else if (mergeDaughter)
            {
                vertexSeedClusterList.erase(iterJ);
                iterJ = vertexSeedClusterList.begin();
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pClusterI, pClusterJ));
            }
            else
            {
                ++iterJ;
            }
        }
    }
}
*/

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexSeedFindingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "SeedClusterListName", m_seedClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NonSeedClusterListName", m_nonSeedClusterListName));

    float minClusterLength = 3.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "MinClusterLength", minClusterLength));
    m_minClusterLengthSquared = minClusterLength * minClusterLength;

    m_minClusterLayers = 5;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "MinClusterLayers", m_minClusterLayers));

    m_minClusterHitsNode = 5;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "MinClusterHitsNode", m_minClusterHitsNode));

    m_minClusterHitsEmission = 15;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "MinClusterHitsEmission", m_minClusterHitsEmission));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
