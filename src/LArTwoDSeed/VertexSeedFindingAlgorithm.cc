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

void VertexSeedFindingAlgorithm::GetSeedClusterList(const ClusterVector &candidateClusters, ClusterList &seedClusterList) const
{
    LArPointingClusterMap pointingClusterMap;

    for (ClusterVector::const_iterator iter = candidateClusters.begin(), iterEnd = candidateClusters.end(); iter != iterEnd; ++iter)
    {
        pointingClusterMap.insert(LArPointingClusterMap::value_type(*iter, LArPointingCluster(*iter)));
    }

    // Identify nodes and emissions
    LArPointingClusterVertexList emissions, associations, seeds;
    const CartesianVector &eventVertex(LArVertexHelper::GetCurrentVertex());

    for (LArPointingClusterMap::const_iterator iter = pointingClusterMap.begin(), iterEnd = pointingClusterMap.end(); iter != iterEnd; ++iter)
    {
        const LArPointingCluster &pointingCluster(iter->second);
        const LArPointingCluster::Vertex &innerVertex(pointingCluster.GetInnerVertex());
        const LArPointingCluster::Vertex &outerVertex(pointingCluster.GetOuterVertex());

        const bool possibleNode(pointingCluster.GetCluster()->GetNCaloHits() >= m_minClusterHitsNode);
        const bool possibleEmission(pointingCluster.GetCluster()->GetNCaloHits() >= m_minClusterHitsEmission);

        const float innerDisplacementSquared((eventVertex - innerVertex.GetPosition()).GetMagnitudeSquared());
        const float outerDisplacementSquared((eventVertex - outerVertex.GetPosition()).GetMagnitudeSquared());

        if (innerDisplacementSquared < outerDisplacementSquared)
        {
            if (LArPointingClusterHelper::IsNode(eventVertex, innerVertex))
            {
                if (possibleNode)
                    seeds.push_back(innerVertex);
            }
            else if(LArPointingClusterHelper::IsEmission(eventVertex, innerVertex))
            {
                if (possibleEmission)
                    emissions.push_back(innerVertex);
            }
        }
        else
        {
            if (LArPointingClusterHelper::IsNode(eventVertex, outerVertex))
            {
                if (possibleNode)
                    seeds.push_back(outerVertex);
            }
            else if(LArPointingClusterHelper::IsEmission(eventVertex, outerVertex))
            {
                if (possibleEmission)
                    emissions.push_back(outerVertex);
            }
        }
    }

    associations.insert(associations.end(), seeds.begin(), seeds.end());
    associations.insert(associations.end(), emissions.begin(), emissions.end());

    // Sort through the list of emissions
    for (LArPointingClusterVertexList::const_iterator iterI = emissions.begin(), iterEndI = emissions.end(); iterI != iterEndI; ++iterI)
    {
        const LArPointingCluster::Vertex &vertexI(*iterI);
        bool isSeed(true);

        for (LArPointingClusterVertexList::const_iterator iterJ = associations.begin(), iterEndJ = associations.end(); iterJ != iterEndJ; ++iterJ)
        {
            const LArPointingCluster::Vertex &vertexJ(*iterJ);

            if (vertexI.GetCluster() == vertexJ.GetCluster())
                continue;

            // Check proximity between vertexI and vertexJ
            if ((2 * vertexI.GetCluster()->GetOrderedCaloHitList().size() < vertexJ.GetCluster()->GetOrderedCaloHitList().size()) &&
                ((vertexI.GetPosition() - vertexJ.GetPosition()).GetMagnitudeSquared() > 5.f * 5.f) &&
                (LArClusterHelper::GetClosestDistance(vertexI.GetPosition(),vertexJ.GetCluster()) < 2.5f))
            {
                isSeed = false;
                break;
            }

            // Check proximity between vertexI and endJ
            LArPointingClusterMap::const_iterator lookupJ = pointingClusterMap.find(vertexJ.GetCluster());

            if (pointingClusterMap.end() == lookupJ)
                throw pandora::StatusCodeException(STATUS_CODE_FAILURE);

            const LArPointingCluster &clusterJ(lookupJ->second);
            const LArPointingCluster::Vertex &endJ(vertexJ.IsInnerVertex() ? clusterJ.GetOuterVertex() : clusterJ.GetInnerVertex());

            if (LArPointingClusterHelper::IsNode(endJ.GetPosition(), vertexI) ||
                (LArPointingClusterHelper::IsEmission(endJ.GetPosition(), vertexI) &&
                (endJ.GetDirection().GetDotProduct(vertexI.GetDirection()) < -0.707f)) )
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
        seedClusterList.insert(iter->GetCluster());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexSeedFindingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_minClusterHitsNode = 5;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "MinClusterHitsNode", m_minClusterHitsNode));

    m_minClusterHitsEmission = 15;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "MinClusterHitsEmission", m_minClusterHitsEmission));

    return SeedFindingBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
