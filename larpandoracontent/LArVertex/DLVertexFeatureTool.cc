/**
 *  @file   larpandoracontent/LArVertex/DLVertexFeatureTool.cc
 *
 *  @brief  Implementation of the DLVertex feature tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArVertex/DLVertexFeatureTool.h"

using namespace pandora;

namespace lar_content
{

DLVertexFeatureTool::DLVertexFeatureTool() :
    m_DLVertexFeatureConstant(17),
    m_numClusterCaloHitsPar(5),
    m_minRDist(1e-9)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLVertexFeatureTool::Run(LArMvaHelper::MvaFeatureVector &featureVector, const VertexSelectionBaseAlgorithm *const pAlgorithm, const Vertex * const pVertex,
    const VertexSelectionBaseAlgorithm::SlidingFitDataListMap &, const VertexSelectionBaseAlgorithm::ClusterListMap &clusterListMap,
    const VertexSelectionBaseAlgorithm::KDTreeMap &, const VertexSelectionBaseAlgorithm::ShowerClusterListMap &, const float, float &)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    const VertexList *pInputVertexList(nullptr);

    if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*pAlgorithm, m_inputVertexListName.c_str(), pInputVertexList))
    {
        featureVector.push_back(0.f);
    }
    else
    {
        const Vertex *const pDLVertex(pInputVertexList->front());
        const float score((pVertex->GetPosition()-pDLVertex->GetPosition()).GetMagnitude());

        const ClusterList &clustersW(clusterListMap.at(TPC_VIEW_W));
        ClusterList longClusters;
        for (const Cluster *const pCluster : clustersW)
            if (pCluster->GetNCaloHits() >= m_numClusterCaloHitsPar) 
                longClusters.push_back(pCluster);

        if (longClusters.empty())
        {
            featureVector.push_back(0.f);
        }
        else
        {
            CartesianVector inner(0.f, 0.f, 0.f), outer(0.f, 0.f, 0.f);
            LArClusterHelper::GetExtremalCoordinates(longClusters, inner, outer);
            double rDist = std::max(std::fabs(outer.GetZ()-inner.GetZ()), std::fabs(outer.GetX()-inner.GetX()));
            rDist = std::max(rDist,m_minRDist);

            const float DLVertexFeature(-m_DLVertexFeatureConstant*score/rDist);

            featureVector.push_back(DLVertexFeature);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLVertexFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DLVertexFeatureConstant", m_DLVertexFeatureConstant));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "InputVertexListName", m_inputVertexListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NumClusterCaloHitsPar", m_numClusterCaloHitsPar));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinRDist", m_minRDist));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
