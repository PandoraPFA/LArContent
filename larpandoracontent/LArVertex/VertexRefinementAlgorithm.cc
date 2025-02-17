/**
 *  @file   larpandoracontent/LArVertex/VertexRefinementAlgorithm.cc
 *
 *  @brief  Implementation of the vertex refinement algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArEigenHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArVertex/VertexRefinementAlgorithm.h"

#define _USE_MATH_DEFINES
#include <cmath>

using namespace pandora;

namespace lar_content
{

VertexRefinementAlgorithm::VertexRefinementAlgorithm() :
    m_caloHitListName{"CaloHitList2D"},
    m_hitRadii(10.f),
    m_vetoPrimaryRegion{true}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexRefinementAlgorithm::Run()
{
    const VertexList *pOutputVertexList{nullptr};
    std::string originalListName, temporaryListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Vertex>(*this, originalListName));

    const VertexList *pInputVertexList{nullptr};
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_inputVertexListName, pInputVertexList));

    if (!pInputVertexList || pInputVertexList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "VertexRefinementAlgorithm: unable to find input vertex list " << std::endl;

        return STATUS_CODE_SUCCESS;
    }
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    if (!pCaloHitList || pCaloHitList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "VertexRefinementAlgorithm: unable to find calo hit list \'" << m_caloHitListName << "\'" << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pOutputVertexList, temporaryListName));

    this->RefineVertices(*pInputVertexList, *pCaloHitList);

    if (!pOutputVertexList->empty())
    {

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_outputVertexListName));
    }
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, originalListName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexRefinementAlgorithm::RefineVertices(const VertexList &vertexList, const CaloHitList &caloHitList) const
{
    const VertexList *pPrimaryVertexList{nullptr};
    PandoraContentApi::GetList(*this, m_primaryVertexListName, pPrimaryVertexList);

    if (m_vetoPrimaryRegion && (!pPrimaryVertexList || pPrimaryVertexList->empty()))
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "VertexRefinementAlgorithm: unable to find primary vertex list " << std::endl;

        return;
    }
    const CartesianVector &primaryVertex{pPrimaryVertexList->front()->GetPosition()};

    CaloHitVector caloHitVectorU, caloHitVectorV, caloHitVectorW;
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        switch (pCaloHit->GetHitType())
        {
            case TPC_VIEW_U:
                caloHitVectorU.emplace_back(pCaloHit);
                break;
            case TPC_VIEW_V:
                caloHitVectorV.emplace_back(pCaloHit);
                break;
            case TPC_VIEW_W:
                caloHitVectorW.emplace_back(pCaloHit);
                break;
            default:
                break;
        }
    }

    for (const Vertex *const pVertex : vertexList)
    {
        const CartesianVector originalPosition(pVertex->GetPosition());
        if (m_vetoPrimaryRegion && ((originalPosition - primaryVertex).GetMagnitude() < m_hitRadii))
            continue;
        const CartesianVector &originalVtxU{LArGeometryHelper::ProjectPosition(this->GetPandora(), originalPosition, TPC_VIEW_U)};
        const CartesianVector &originalVtxV{LArGeometryHelper::ProjectPosition(this->GetPandora(), originalPosition, TPC_VIEW_V)};
        const CartesianVector &originalVtxW{LArGeometryHelper::ProjectPosition(this->GetPandora(), originalPosition, TPC_VIEW_W)};

        CaloHitList nearbyHitListU;
        this->GetNearbyHits(caloHitVectorU, originalVtxU, nearbyHitListU);
        const CartesianVector vtxU(this->RefineVertexTwoD(nearbyHitListU, originalVtxU));
        CaloHitList nearbyHitListV;
        this->GetNearbyHits(caloHitVectorV, originalVtxV, nearbyHitListV);
        const CartesianVector vtxV(this->RefineVertexTwoD(nearbyHitListV, originalVtxV));
        CaloHitList nearbyHitListW;
        this->GetNearbyHits(caloHitVectorW, originalVtxW, nearbyHitListW);
        const CartesianVector vtxW(this->RefineVertexTwoD(nearbyHitListW, originalVtxW));

        CartesianVector vtxUV(0.f, 0.f, 0.f), vtxUW(0.f, 0.f, 0.f), vtxVW(0.f, 0.f, 0.f), vtx3D(0.f, 0.f, 0.f), position3D(0.f, 0.f, 0.f);
        float chi2UV(0.f), chi2UW(0.f), chi2VW(0.f), chi23D(0.f);

        LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, vtxU, vtxV, vtxUV, chi2UV);
        LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_W, vtxU, vtxW, vtxUW, chi2UW);
        LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), TPC_VIEW_V, TPC_VIEW_W, vtxV, vtxW, vtxVW, chi2VW);
        LArGeometryHelper::MergeThreePositions3D(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, vtxU, vtxV, vtxW, vtx3D, chi23D);

        if (chi2UV < chi2UW && chi2UV < chi2VW && chi2UV < chi23D)
            position3D = vtxUV;
        else if (chi2UW < chi2VW && chi2UW < chi23D)
            position3D = vtxUW;
        else if (chi2VW < chi23D)
            position3D = vtxVW;
        else
            position3D = vtx3D;

        PandoraContentApi::Vertex::Parameters parameters;
        parameters.m_position = position3D;
        parameters.m_vertexLabel = VERTEX_INTERACTION;
        parameters.m_vertexType = VERTEX_3D;

        const Vertex *pNewVertex(NULL);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pNewVertex));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexRefinementAlgorithm::GetNearbyHits(const CaloHitVector &hitVector, const CartesianVector &centroid, CaloHitList &nearbyHitList) const
{
    Eigen::MatrixXf hitMatrix(hitVector.size(), 2);
    LArEigenHelper::Vectorize(hitVector, hitMatrix);
    Eigen::RowVectorXf vertex(2);
    vertex << centroid.GetX(), centroid.GetZ();
    Eigen::MatrixXf norms((hitMatrix.rowwise() - vertex).array().pow(2).rowwise().sum());
    for (int r = 0; r < hitMatrix.rows(); ++r)
    {
        if (norms(r, 0) < m_hitRadii * m_hitRadii)
            nearbyHitList.emplace_back(hitVector.at(r));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector VertexRefinementAlgorithm::RefineVertexTwoD(const CaloHitList &caloHitList, const CartesianVector &seedVertex) const
{
    if (caloHitList.empty())
        return seedVertex;
    const int nBins{72};
    const float pi{static_cast<float>(M_PI)};
    Eigen::MatrixXf hitMatrix(caloHitList.size(), 2);
    Eigen::MatrixXf loMatrix(caloHitList.size(), 2);
    Eigen::MatrixXf hiMatrix(caloHitList.size(), 2);
    LArEigenHelper::Vectorize(caloHitList, hitMatrix, loMatrix, hiMatrix);
    Eigen::RowVectorXf results{Eigen::RowVectorXf::Zero(hitMatrix.rows())};
    float best{0.f};
    for (int r = 0; r < hitMatrix.rows(); ++r)
    {
        Eigen::RowVectorXf row(2);
        row << hitMatrix(r, 0), hitMatrix(r, 1);
        // Compute dx, dz and angle between each hit and the candidate vertex hit
        Eigen::RowVectorXf loPhis(loMatrix.rows());
        Eigen::RowVectorXf hiPhis(hiMatrix.rows());
        LArEigenHelper::GetAngles(loMatrix, row, loPhis);
        LArEigenHelper::GetAngles(hiMatrix, row, hiPhis);
        loPhis = loPhis.array() * nBins / (2 * pi);
        hiPhis = hiPhis.array() * nBins / (2 * pi);

        Eigen::RowVectorXf counts{Eigen::RowVectorXf::Zero(nBins)};
        Eigen::RowVectorXf counts2{Eigen::RowVectorXf::Zero(nBins)};
        for (int i = 0; i < hitMatrix.rows(); ++i)
        {
            if (i == r)
                continue;
            const float lo{loPhis(i) < hiPhis(i) ? loPhis(i) : hiPhis(i)};
            const float hi{loPhis(i) < hiPhis(i) ? hiPhis(i) : loPhis(i)};
            const float arc{hi - lo};
            const int bin0{loPhis(i) < hiPhis(i) ? static_cast<int>(loPhis(i)) : static_cast<int>(hiPhis(i))};
            const int bin1{loPhis(i) < hiPhis(i) ? static_cast<int>(hiPhis(i)) : static_cast<int>(loPhis(i))};
            for (int b = bin0; b <= bin1; ++b)
            {
                const float theta0{std::max(lo, static_cast<float>(b))}, theta1{std::min(hi, static_cast<float>(b + 1))};
                const float frac{(theta1 - theta0) / arc};
                counts(b) += frac;                                                       // * rWeightVector(i);
                counts2(b > (nBins >> 1) ? b - (nBins >> 1) : b + (nBins >> 1)) -= frac; // * rWeightVector(i);
            }
        }
        const CartesianVector centroid(hitMatrix(r, 0), 0, hitMatrix(r, 1));
        const float shift{(centroid - seedVertex).GetMagnitude()};
        if (shift > 3.f)
            results(r) = (counts + counts2).array().pow(2).sum() / (shift / 3.f);
        else
            results(r) = (counts + counts2).array().pow(2).sum();
        if (results(r) > best)
        {
            best = results(r);
            /*            std::cout << "Total: " << results(r) << std::endl;

            PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitList, "near", BLACK));
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &centroid, "vtx", RED, 3));
            for (int i = 0; i < nBins; ++i)
            {
                const float theta{2.f * i * pi / nBins};
                const CartesianVector b{CartesianVector(std::cos(theta), 0, std::sin(theta)) * 2 * m_hitRadii + centroid};
                PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &centroid, &b, "bin " + std::to_string(i), GRAY, 1, 1));
            }
            PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));*/
        }
    }
    Eigen::Index index;
    results.maxCoeff(&index);

    const CartesianVector centroid(hitMatrix(index, 0), 0, hitMatrix(index, 1));
    /*PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &seedVertex, "seed", MAGENTA, 1));
    PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &centroid, "vtx", RED, 1));
    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitList, "near", BLUE));
    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));*/

    return centroid;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexRefinementAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputVertexListName", m_inputVertexListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputVertexListName", m_outputVertexListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PrimaryVertexListName", m_primaryVertexListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "HitRadii", m_hitRadii));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VetoPrimaryRegion", m_vetoPrimaryRegion));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
