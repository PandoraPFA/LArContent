/**
 *  @file   larpandoracontent/MyDlVtxFeatureTool.cc
 *
 *  @brief  Implementation of the myDlVtx feature tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArVertex/MyDlVtxFeatureTool.h"

using namespace pandora;

namespace lar_content
{

MyDlVtxFeatureTool::MyDlVtxFeatureTool() :
    m_myDlVtxConstant(17),
    m_inputVertexListName(),
    m_numClusterCaloHitsPar(5)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MyDlVtxFeatureTool::Run(LArMvaHelper::MvaFeatureVector &featureVector, const VertexSelectionBaseAlgorithm *const pAlgorithm, const Vertex * const pVertex,
    const VertexSelectionBaseAlgorithm::SlidingFitDataListMap &/*slidingFitDataListMap*/, const VertexSelectionBaseAlgorithm::ClusterListMap &clusterListMap,
    const VertexSelectionBaseAlgorithm::KDTreeMap &, const VertexSelectionBaseAlgorithm::ShowerClusterListMap &, const float, float &)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    const VertexList *pInputVertexList(NULL);
    std::string vertexListName(m_inputVertexListName);
    const StatusCode statusCode(PandoraContentApi::GetList(*pAlgorithm, vertexListName.c_str(), pInputVertexList));

    if(STATUS_CODE_SUCCESS != statusCode)
        featureVector.push_back(0.f);
    else
    {
        const Vertex * const pDlVertex(pInputVertexList->front());
        float score((pVertex->GetPosition()-pDlVertex->GetPosition()).GetMagnitude());

        ClusterList clustersW(clusterListMap.at(TPC_VIEW_W));
        ClusterList sortedLongClusters;
        for (const Cluster *const pCluster : clustersW)
            if (pCluster->GetNCaloHits() >= m_numClusterCaloHitsPar) 
                sortedLongClusters.push_back(pCluster);

        if(sortedLongClusters.empty())
            featureVector.push_back(0.f);
        else
        {
            CartesianVector inner(0.f, 0.f, 0.f), outer(0.f, 0.f, 0.f);
            LArClusterHelper::GetExtremalCoordinates(sortedLongClusters, inner, outer);
            double rdist = std::max(fabs(outer.GetZ()-inner.GetZ()), fabs(outer.GetX()-inner.GetX()));
            rdist = std::max(rdist,(1e-9));

            float MyDlVtxFeature(-m_myDlVtxConstant*score/rdist);

            featureVector.push_back(MyDlVtxFeature);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MyDlVtxFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MyDlVtxConstant", m_myDlVtxConstant));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "InputVertexListName", m_inputVertexListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NumClusterCaloHitsPar", m_numClusterCaloHitsPar));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
