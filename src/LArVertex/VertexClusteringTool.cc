/**
 *  @file   ExampleContent/src/VertexClusteringTools/VertexClusteringTool.cc
 * 
 *  @brief  Implementation of the example algorithm tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArVertex/VertexClusteringTool.h"
#include "LArHelpers/LArMCParticleHelper.h"
#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"

using namespace pandora;

namespace lar_content
{

VertexClusteringTool::VertexClusteringTool() :
    m_clusteringRadius(2.5f),
    m_maxVertexToCentroidDistance(5.f),
    m_minClusterSize(5)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexClusteringTool::SortVerticesByZ(const pandora::Vertex *const pLhs, const pandora::Vertex *const pRhs)
{
    const CartesianVector deltaPosition(pRhs->GetPosition() - pLhs->GetPosition());

    if (std::fabs(deltaPosition.GetZ()) > std::numeric_limits<float>::epsilon())
        return (deltaPosition.GetZ() > std::numeric_limits<float>::epsilon());

    if (std::fabs(deltaPosition.GetX()) > std::numeric_limits<float>::epsilon())
        return (deltaPosition.GetX() > std::numeric_limits<float>::epsilon());

    // ATTN No way to distinguish between vertices if still have a tie in y coordinate
    return (deltaPosition.GetY() > std::numeric_limits<float>::epsilon());
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<const VertexList*> VertexClusteringTool::ClusterVertices(const VertexList* pVertexList)
{
    if (this->GetPandora().GetSettings()->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << this->GetType() << std::endl;

    std::vector<const Vertex*> sortedVertexVector;
    
    for (const Vertex *const pVertex : (*pVertexList))
        sortedVertexVector.push_back(pVertex);

    std::sort(sortedVertexVector.begin(), sortedVertexVector.end(), SortVerticesByZ);

    VertexClusterList vertexClusterList;
    VertexList usedVertices;
    
    VertexCluster *const pVertexClusterSeed = new VertexCluster();
    CartesianVector currentClusterCentroid(0., 0., 0.);
    
    for (const Vertex *const pVertex : sortedVertexVector)
    {
        if (usedVertices.count(pVertex) != 0)
            continue;
        
        if (pVertexClusterSeed->GetVertexList().size() == 0)
        {
            pVertexClusterSeed->AddVertex(pVertex);
            usedVertices.insert(pVertex);
        }

        for (const Vertex *const pVertexIterator : sortedVertexVector)
        {
            if (usedVertices.count(pVertexIterator) != 0)
                continue;
            
            currentClusterCentroid = pVertexClusterSeed->GetCentroidPosition();
            
            if ((pVertexClusterSeed->GetVertexList().size() <= m_minClusterSize) && (CheckVertexToClusterDistance(pVertexIterator, pVertexClusterSeed)))
            {
                pVertexClusterSeed->AddVertex(pVertexIterator);
                usedVertices.insert(pVertexIterator);
            }
            else if ((pVertexClusterSeed->GetVertexList().size() > m_minClusterSize) && ((((currentClusterCentroid - (pVertexIterator->GetPosition())).GetMagnitude() <= m_maxVertexToCentroidDistance)) && (CheckVertexToClusterDistance(pVertexIterator, pVertexClusterSeed))))
            {
                pVertexClusterSeed->AddVertex(pVertexIterator);
                usedVertices.insert(pVertexIterator);
            }
        }
        
        VertexCluster *const pVertexCluster = new VertexCluster(*pVertexClusterSeed);
        
        vertexClusterList.push_back(pVertexCluster);
        pVertexClusterSeed->ClearVertexCluster();
        
        //pVertexClusterSeed->AddVertex(pVertex);
        //usedVertices.insert(pVertex);
    }
    
    if (pVertexClusterSeed->GetVertexList().size() != 0)
    {
        VertexCluster *const pLastVertexCluster = new VertexCluster(*pVertexClusterSeed);
        vertexClusterList.push_back(pLastVertexCluster);
    }
    
    std::vector<const VertexList*> outputVertexListVector;
    
    for (VertexCluster* pVertexCluster : vertexClusterList)
        outputVertexListVector.push_back(&(pVertexCluster->GetVertexList())); //all for now: later m_nSelectedVerticesPerCluster

    return outputVertexListVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::CartesianVector VertexClusteringTool::VertexCluster::GetCentroidPosition() const
{
    if (this->GetVertexList().empty())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

    pandora::CartesianVector centroid(0.f, 0.f, 0.f);

    for (const pandora::Vertex *const pVertex : this->GetVertexList())
        centroid += pVertex->GetPosition();

    centroid *= static_cast<float>(1.f / this->GetVertexList().size());
    return centroid;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexClusteringTool::CheckVertexToClusterDistance(const Vertex *const pVertex, VertexCluster *const pVertexCluster) const
{
    for (const Vertex *const pClusterVertex : (pVertexCluster->GetVertexList()))
    {
        if (((pVertex->GetPosition() - pClusterVertex->GetPosition()).GetMagnitude()) <= m_clusteringRadius)
            return true;
    }
    
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexClusteringTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    // Read settings from xml file here
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusteringRadius", m_clusteringRadius));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxVertexToCentroidDistance", m_maxVertexToCentroidDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterSize", m_minClusterSize));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
