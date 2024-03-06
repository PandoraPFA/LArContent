/**
 *  @file   larpandoracontent/LArVertex/VertexRefinementAlgorithm.cc
 *
 *  @brief  Implementation of the vertex refinement algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"

#include "larpandoracontent/LArVertex/VertexRefinementAlgorithm.h"

#include <Eigen/Dense>

using namespace pandora;

namespace lar_content
{

VertexRefinementAlgorithm::VertexRefinementAlgorithm() :
    m_chiSquaredCut(2.f), m_distanceCut(5.f), m_minimumHitsCut(5), m_twoDDistanceCut(10.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexRefinementAlgorithm::Run()
{
    const VertexList *pInputVertexList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pInputVertexList));

    if (!pInputVertexList || pInputVertexList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "VertexRefinementAlgorithm: unable to find current vertex list " << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    const VertexList *pOutputVertexList(NULL);
    std::string temporaryListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pOutputVertexList, temporaryListName));

    ClusterList clusterListU, clusterListV, clusterListW;
    this->GetClusterLists(m_inputClusterListNames, clusterListU, clusterListV, clusterListW);

    this->RefineVertices(pInputVertexList, clusterListU, clusterListV, clusterListW);

    if (!pOutputVertexList->empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_outputVertexListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_outputVertexListName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexRefinementAlgorithm::GetClusterLists(
    const StringVector &inputClusterListNames, ClusterList &clusterListU, ClusterList &clusterListV, ClusterList &clusterListW) const
{
    for (const std::string &clusterListName : inputClusterListNames)
    {
        const ClusterList *pClusterList(nullptr);
        PANDORA_THROW_RESULT_IF_AND_IF(
            STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));

        if (!pClusterList || pClusterList->empty())
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "VertexRefinementAlgorithm: unable to find cluster list " << clusterListName << std::endl;

            continue;
        }

        const HitType hitType(LArClusterHelper::GetClusterHitType(pClusterList->front()));

        if ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        ClusterList &outputClusterList((TPC_VIEW_U == hitType) ? clusterListU : (TPC_VIEW_V == hitType) ? clusterListV : clusterListW);
        outputClusterList.insert(outputClusterList.end(), pClusterList->begin(), pClusterList->end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexRefinementAlgorithm::RefineVertices(const VertexList *const pVertexList, const ClusterList &clusterListU,
    const ClusterList &clusterListV, const ClusterList &clusterListW) const
{
    for (const Vertex *const pVertex : *pVertexList)
    {
        const CartesianVector originalPosition(pVertex->GetPosition());

        const CartesianVector vtxU(
            this->RefineVertexTwoD(clusterListU, LArGeometryHelper::ProjectPosition(this->GetPandora(), originalPosition, TPC_VIEW_U)));
        const CartesianVector vtxV(
            this->RefineVertexTwoD(clusterListV, LArGeometryHelper::ProjectPosition(this->GetPandora(), originalPosition, TPC_VIEW_V)));
        const CartesianVector vtxW(
            this->RefineVertexTwoD(clusterListW, LArGeometryHelper::ProjectPosition(this->GetPandora(), originalPosition, TPC_VIEW_W)));

        CartesianVector vtxUV(0.f, 0.f, 0.f), vtxUW(0.f, 0.f, 0.f), vtxVW(0.f, 0.f, 0.f), vtx3D(0.f, 0.f, 0.f), position3D(0.f, 0.f, 0.f);
        float chi2UV(0.f), chi2UW(0.f), chi2VW(0.f), chi23D(0.f), chi2(0.f);

        LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, vtxU, vtxV, vtxUV, chi2UV);
        LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_W, vtxU, vtxW, vtxUW, chi2UW);
        LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), TPC_VIEW_V, TPC_VIEW_W, vtxV, vtxW, vtxVW, chi2VW);
        LArGeometryHelper::MergeThreePositions3D(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, vtxU, vtxV, vtxW, vtx3D, chi23D);

        if (chi2UV < chi2UW && chi2UV < chi2VW && chi2UV < chi23D)
        {
            position3D = vtxUV;
            chi2 = chi2UV;
        }
        else if (chi2UW < chi2VW && chi2UW < chi23D)
        {
            position3D = vtxUW;
            chi2 = chi2UW;
        }
        else if (chi2VW < chi23D)
        {
            position3D = vtxVW;
            chi2 = chi2VW;
        }
        else
        {
            position3D = vtx3D;
            chi2 = chi23D;
        }

        if (chi2 > m_chiSquaredCut)
            position3D = originalPosition;

        if ((position3D - originalPosition).GetMagnitude() > m_distanceCut)
            position3D = originalPosition;

        PandoraContentApi::Vertex::Parameters parameters;
        parameters.m_position = position3D;
        parameters.m_vertexLabel = VERTEX_INTERACTION;
        parameters.m_vertexType = VERTEX_3D;

        const Vertex *pNewVertex(NULL);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pNewVertex));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector VertexRefinementAlgorithm::RefineVertexTwoD(const ClusterList &clusterList, const CartesianVector &originalVtxPos) const
{
    CartesianPointVector intercepts, directions;
    FloatVector weights;

    for (const Cluster *const pCluster : clusterList)
    {
        if (LArClusterHelper::GetClosestDistance(originalVtxPos, pCluster) > 10)
            continue;

        if (pCluster->GetNCaloHits() < m_minimumHitsCut)
            continue;

        CartesianVector centroid(0.f, 0.f, 0.f);
        LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
        LArPcaHelper::EigenVectors eigenVectors;
        CartesianPointVector pointVector;

        LArClusterHelper::GetCoordinateVector(pCluster, pointVector);
        LArPcaHelper::RunPca(pointVector, centroid, eigenValues, eigenVectors);

        intercepts.push_back(LArClusterHelper::GetClosestPosition(originalVtxPos, pCluster));
        directions.push_back(eigenVectors.at(0).GetUnitVector());
        weights.push_back(1.f / ((LArClusterHelper::GetClosestPosition(originalVtxPos, pCluster) - originalVtxPos).GetMagnitudeSquared() + 1));
    }

    CartesianVector newVtxPos(originalVtxPos);
    if (intercepts.size() > 1)
        GetBestFitPoint(intercepts, directions, weights, newVtxPos);

    if ((newVtxPos - originalVtxPos).GetMagnitude() > m_twoDDistanceCut)
        return originalVtxPos;

    return newVtxPos;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexRefinementAlgorithm::GetBestFitPoint(const CartesianPointVector &intercepts, const CartesianPointVector &directions,
    const FloatVector &weights, CartesianVector &bestFitPoint) const
{
    const int n(intercepts.size());

    Eigen::VectorXd d = Eigen::VectorXd::Zero(3 * n);
    Eigen::MatrixXd G = Eigen::MatrixXd::Zero(3 * n, n + 3);

    for (int i = 0; i < n; ++i)
    {
        d(3 * i) = intercepts[i].GetX() * weights[i];
        d(3 * i + 1) = intercepts[i].GetY() * weights[i];
        d(3 * i + 2) = intercepts[i].GetZ() * weights[i];

        G(3 * i, 0) = weights[i];
        G(3 * i + 1, 1) = weights[i];
        G(3 * i + 2, 2) = weights[i];

        G(3 * i, i + 3) = -directions[i].GetX() * weights[i];
        G(3 * i + 1, i + 3) = -directions[i].GetY() * weights[i];
        G(3 * i + 2, i + 3) = -directions[i].GetZ() * weights[i];
    }

    if ((G.transpose() * G).determinant() < std::numeric_limits<float>::epsilon())
    {
        bestFitPoint = CartesianVector(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
        return;
    }

    Eigen::VectorXd m = (G.transpose() * G).inverse() * G.transpose() * d;

    bestFitPoint = CartesianVector(m[0], m[1], m[2]);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexRefinementAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputClusterListNames", m_inputClusterListNames));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputVertexListName", m_inputVertexListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputVertexListName", m_outputVertexListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ChiSquaredCut", m_chiSquaredCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "DistanceCut", m_distanceCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinimumHitsCut", m_minimumHitsCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TwoDDistanceCut", m_twoDDistanceCut));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
