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
    m_maxVertexToCentroidDistance(5.f),
    m_removeSmallClusters(false),
    m_minClusterSize(5),
    m_monteCarloClusterCheck(false),
    m_recoEndPointCheck(false)
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

std::vector<const VertexList*> VertexClusteringTool::ClusterVertices(const Algorithm *const pAlgorithm, const VertexList* pVertexList)
{
    if (this->GetPandora().GetSettings()->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << this->GetType() << std::endl;

    //-------------------------------------------------------------------------------------------------

    const MCParticleList *pMCParticleList = NULL;
    PandoraContentApi::GetList(*pAlgorithm, "MCParticleList3D", pMCParticleList);
    MCParticleVector mcPrimaryVector;
    LArMCParticleHelper::GetPrimaryMCParticleList(pMCParticleList, mcPrimaryVector);

    std::vector<CartesianVector> endpointVector;

    for (const MCParticle* mcParticle : (*pMCParticleList))
    {
        endpointVector.push_back(mcParticle->GetVertex());   
        endpointVector.push_back(mcParticle->GetEndpoint()); 
    }

   //-------------------------------------------------------------------------------------------------

//    const ClusterList *pClusterListU(NULL);
//    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*pAlgorithm, "ClustersU", pClusterListU));

//    const ClusterList *pClusterListV(NULL);
//    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*pAlgorithm, "ClustersV", pClusterListV));

    const ClusterList *pClusterListW(NULL);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*pAlgorithm, "ClustersW", pClusterListW));

    std::vector<CartesianVector> reconstructedEndpointVectorW;

    for (const Cluster* pCluster : (*pClusterListW))
    {
        OrderedCaloHitList orderedCaloHitList(pCluster->GetOrderedCaloHitList());
        CartesianVector innerCoordinate(0.f, 0.f, 0.f);
        CartesianVector outerCoordinate(0.f, 0.f, 0.f);

        LArClusterHelper::GetExtremalCoordinates(orderedCaloHitList, innerCoordinate, outerCoordinate);
        reconstructedEndpointVectorW.push_back(innerCoordinate);   
        reconstructedEndpointVectorW.push_back(outerCoordinate); 
    }

    //-------------------------------------------------------------------------------------------------

    std::vector<const Vertex*> sortedVertexVector;
    for (const Vertex *const pVertex : (*pVertexList))
        sortedVertexVector.push_back(pVertex);

    std::sort(sortedVertexVector.begin(), sortedVertexVector.end(), SortVerticesByZ);

    VertexClusterList vertexClusterList;
    VertexList usedVertices;
    
    VertexCluster *const pVertexClusterSeed = new VertexCluster();
    const Vertex* pPreviousVertex(*(sortedVertexVector.begin()));

    for (const Vertex *const pVertex : sortedVertexVector)
    {
        if (usedVertices.count(pVertex) == 1)
            continue;
        
        CartesianVector currentClusterCentroid(0., 0., 0.);

        if (pVertexClusterSeed->GetVertexList().size() == 0)
        {
            pVertexClusterSeed->AddVertex(pVertex);
            usedVertices.insert(pVertex);
            continue;
        }
        else
            currentClusterCentroid = pVertexClusterSeed->GetCentroidPosition();
        
        if ((pVertexClusterSeed->GetVertexList().size() <= m_minClusterSize) && (((pVertex->GetPosition() - pPreviousVertex->GetPosition()).GetMagnitude()) <= 2.5))
        {
            pVertexClusterSeed->AddVertex(pVertex);
            usedVertices.insert(pVertex);
        }
        else if ((pVertexClusterSeed->GetVertexList().size() > m_minClusterSize) && ((((currentClusterCentroid - (pVertex->GetPosition())).GetMagnitude() <= m_maxVertexToCentroidDistance)) && (((pVertex->GetPosition() - pPreviousVertex->GetPosition()).GetMagnitude()) <= 2.5)))
        {
            pVertexClusterSeed->AddVertex(pVertex);
            usedVertices.insert(pVertex);
        }
        else
        {
            VertexCluster *const pNewVertexCluster = new VertexCluster(*pVertexClusterSeed);

            if (m_monteCarloClusterCheck)
            {
                bool MCCheckFulfilled(false);
                for (const Vertex *const pAnotherVertex : (pNewVertexCluster->GetVertexList()))
                {
                    if (MCCheckFulfilled)
                        break;

                    for (CartesianVector &endPoint : endpointVector)
                    {
                        if (((endPoint - (pAnotherVertex->GetPosition())).GetMagnitude()) < m_maxVertexToCentroidDistance)
                        {
                            //std::cout << "True cluster." << std::endl;
                            vertexClusterList.push_back(pNewVertexCluster);
                            MCCheckFulfilled = true;
                            break;  
                        }
                    }
                }
            }
            else if (m_recoEndPointCheck)
            {
                bool recoCheckFulfilled(false);
                for (const Vertex *const pYetAnotherVertex : (pNewVertexCluster->GetVertexList()))
                {
                    if (recoCheckFulfilled)
                        break;

                    const CartesianVector vertexProjectionW(lar_content::LArGeometryHelper::ProjectPosition(this->GetPandora(), pYetAnotherVertex->GetPosition(), TPC_VIEW_W));

                    for (CartesianVector &endPoint : reconstructedEndpointVectorW)
                    {
                        if (((endPoint - vertexProjectionW).GetMagnitude()) < m_maxVertexToCentroidDistance)
                        {
                            //std::cout << "True cluster." << std::endl;
                            vertexClusterList.push_back(pNewVertexCluster);
                            recoCheckFulfilled = true;
                            break;  
                        }
                    }
                }
            }
            else
                vertexClusterList.push_back(pNewVertexCluster);

            pVertexClusterSeed->ClearVertexCluster();
            pVertexClusterSeed->AddVertex(pVertex);
            usedVertices.insert(pVertex);
        }

        pPreviousVertex = pVertex;     
    }

    VertexCluster *const pNewVertexCluster = new VertexCluster(*pVertexClusterSeed);

    if (m_monteCarloClusterCheck)
    {
        for (const Vertex *const pVertex : (pNewVertexCluster->GetVertexList()))
        {
            for (CartesianVector &endPoint : endpointVector)
            {
                if (((endPoint - (pVertex->GetPosition())).GetMagnitude()) < m_maxVertexToCentroidDistance)
                {
                    //std::cout << "True cluster." << std::endl;
                    vertexClusterList.push_back(pNewVertexCluster);
                    break;   
                }
            }
        }
    }
    else if (m_recoEndPointCheck)
    {
        for (const Vertex *const pYetAnotherVertex : (pNewVertexCluster->GetVertexList()))
        {
            const CartesianVector vertexProjectionW(lar_content::LArGeometryHelper::ProjectPosition(this->GetPandora(), pYetAnotherVertex->GetPosition(), TPC_VIEW_W));

            for (CartesianVector &endPoint : reconstructedEndpointVectorW)
            {
                if (((endPoint - vertexProjectionW).GetMagnitude()) < m_maxVertexToCentroidDistance)
                {
                    //std::cout << "True cluster." << std::endl;
                    vertexClusterList.push_back(pNewVertexCluster);
                    break;  
                }
            }
        }
    }
    else
        vertexClusterList.push_back(pNewVertexCluster);

    std::vector<const VertexList*> outputVertexListVector;
    
    if (m_removeSmallClusters == true)
        this->RemoveSmallClusters(vertexClusterList);
    
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

void VertexClusteringTool::RemoveSmallClusters(VertexClusterList &vertexClusterList)
{
    for (VertexClusterList::const_iterator iter = vertexClusterList.begin(); iter != vertexClusterList.end(); )
    {
        const VertexCluster* pVertexCluster(*iter);
        
        if (pVertexCluster->GetVertexList().size() == 1)
            iter = vertexClusterList.erase(iter); //calls destructor
        else
            ++iter;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexClusteringTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    // Read settings from xml file here

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxVertexToCentroidDistance", m_maxVertexToCentroidDistance));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "RemoveSmallClusters", m_removeSmallClusters));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterSize", m_minClusterSize));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MonteCarloClusterCheck", m_monteCarloClusterCheck));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "RecoEndPointCheck", m_recoEndPointCheck));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
