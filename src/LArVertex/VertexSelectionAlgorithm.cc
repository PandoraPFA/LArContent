/**
 *  @file   LArContent/src/LArVertex/VertexSelectionAlgorithm.cc
 * 
 *  @brief  Implementation of the vertex selection algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"

#include "LArVertex/VertexSelectionAlgorithm.h"

using namespace pandora;

namespace lar_content
{

VertexSelectionAlgorithm::VertexSelectionAlgorithm() :
    m_replaceCurrentVertexList(true),
    m_histogramNPhiBins(200),
    m_histogramPhiMin(-1.1f * M_PI),
    m_histogramPhiMax(+1.1f * M_PI),
    m_maxHitVertexDisplacement(std::numeric_limits<float>::max()),
    m_hitDeweightingPower(-0.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexSelectionAlgorithm::Run()
{
    const VertexList *pInputVertexList(NULL);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pInputVertexList));

    Vertex *pBestVertex(NULL);
    float bestFigureOfMerit(0.f);

    for (VertexList::const_iterator iter = pInputVertexList->begin(), iterEnd = pInputVertexList->end(); iter != iterEnd; ++iter)
    {
        Vertex *const pVertex(*iter);
        Histogram histogramU(m_histogramNPhiBins, m_histogramPhiMin, m_histogramPhiMax);
        Histogram histogramV(m_histogramNPhiBins, m_histogramPhiMin, m_histogramPhiMax);
        Histogram histogramW(m_histogramNPhiBins, m_histogramPhiMin, m_histogramPhiMax);

        this->FillHistogram(pVertex, TPC_VIEW_U, m_inputClusterListNameU, histogramU);
        this->FillHistogram(pVertex, TPC_VIEW_V, m_inputClusterListNameV, histogramV);
        this->FillHistogram(pVertex, TPC_VIEW_W, m_inputClusterListNameW, histogramW);

        const float figureOfMerit(this->GetFigureOfMerit(histogramU, histogramV, histogramW));

        if (figureOfMerit > bestFigureOfMerit)
        {
            pBestVertex = pVertex;
            bestFigureOfMerit = figureOfMerit;
        }
    }

    if (NULL != pBestVertex)
    {
        VertexList selectedVertexList;
        selectedVertexList.insert(pBestVertex);

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_outputVertexListName, selectedVertexList));

        if (m_replaceCurrentVertexList)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_outputVertexListName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionAlgorithm::FillHistogram(const Vertex *const pVertex, const HitType hitType, const std::string &clusterListName,
    Histogram &histogram) const
{
    const ClusterList *pClusterList = NULL;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));

    for (ClusterList::const_iterator cIter = pClusterList->begin(), cIterEnd = pClusterList->end(); cIter != cIterEnd; ++cIter)
    {
        const Cluster *const pCluster(*cIter);
        const HitType clusterHitType(LArClusterHelper::GetClusterHitType(pCluster));

        if (hitType != clusterHitType)
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), hitType));
        this->FillHistogram(vertexPosition2D, pCluster, histogram);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionAlgorithm::FillHistogram(const CartesianVector &vertexPosition2D, const Cluster *const pCluster, Histogram &histogram) const
{
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hIter = iter->second->begin(), hIterEnd = iter->second->end(); hIter != hIterEnd; ++hIter)
        {
            const CaloHit *const pCaloHit(*hIter);

            const CartesianVector displacement(pCaloHit->GetPositionVector() - vertexPosition2D);
            const float magnitude(displacement.GetMagnitude());

            if (magnitude > m_maxHitVertexDisplacement)
                continue;

            const float phi(std::atan2(displacement.GetZ(), displacement.GetX()));
            histogram.Fill(phi, std::pow(magnitude, m_hitDeweightingPower));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VertexSelectionAlgorithm::GetFigureOfMerit(const Histogram &histogramU, const Histogram &histogramV, const Histogram &histogramW) const
{
    float figureOfMerit(0.f);
    figureOfMerit += this->GetFigureOfMerit(histogramU);
    figureOfMerit += this->GetFigureOfMerit(histogramV);
    figureOfMerit += this->GetFigureOfMerit(histogramW);

    return figureOfMerit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VertexSelectionAlgorithm::GetFigureOfMerit(const Histogram &histogram) const
{
    float sumSquaredEntries(0.f);

    for (int xBin = 0; xBin < histogram.GetNBinsX(); ++xBin)
    {
        const float binContent(histogram.GetBinContent(xBin));
        sumSquaredEntries += binContent * binContent;
    }

    return sumSquaredEntries;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexSelectionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameU", m_inputClusterListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameV", m_inputClusterListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameW", m_inputClusterListNameW));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputVertexListName", m_outputVertexListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ReplaceCurrentVertexList", m_replaceCurrentVertexList));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "HistogramNPhiBins", m_histogramNPhiBins));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "HistogramPhiMin", m_histogramPhiMin));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "HistogramPhiMax", m_histogramPhiMax));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxHitVertexDisplacement", m_maxHitVertexDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "HitDeweightingPower", m_hitDeweightingPower));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
