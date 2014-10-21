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
    m_maxOnHitDisplacement(1.f),
    m_hitDeweightingPower(-0.5f),
    m_maxTopScoreCandidates(50),
    m_maxTopScoreSelections(3),
    m_minCandidateDisplacement(2.f),
    m_minCandidateScoreFraction(0.4f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexSelectionAlgorithm::Run()
{
    const VertexList *pInputVertexList(NULL);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pInputVertexList));

    VertexScoreList vertexScoreList;

    for (VertexList::const_iterator iter = pInputVertexList->begin(), iterEnd = pInputVertexList->end(); iter != iterEnd; ++iter)
    {
        try
        {
            Vertex *const pVertex(*iter);
            const float figureOfMerit(this->GetFigureOfMerit(pVertex, m_maxHitVertexDisplacement,m_hitDeweightingPower));
            vertexScoreList.push_back(VertexScore(pVertex, figureOfMerit));
        }
        catch (StatusCodeException &statusCodeException)
        {
            if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                throw statusCodeException;
        }
    }

    std::sort(vertexScoreList.begin(), vertexScoreList.end());

    VertexScoreList selectedVertexScoreList;
    this->SelectTopScoreVertices(vertexScoreList, selectedVertexScoreList);

    Vertex *pFinalVertex(NULL);
    this->SelectFinalVertex(selectedVertexScoreList, pFinalVertex);

    if (NULL != pFinalVertex)
    {
        VertexList selectedVertexList;
        selectedVertexList.insert(pFinalVertex);

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_outputVertexListName, selectedVertexList));

        if (m_replaceCurrentVertexList)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_outputVertexListName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VertexSelectionAlgorithm::GetFigureOfMerit(const Vertex *const pVertex, const float maxHitVertexDisplacement, const float hitDeweightingPower) const
{
    Histogram histogramU(m_histogramNPhiBins, m_histogramPhiMin, m_histogramPhiMax);
    Histogram histogramV(m_histogramNPhiBins, m_histogramPhiMin, m_histogramPhiMax);
    Histogram histogramW(m_histogramNPhiBins, m_histogramPhiMin, m_histogramPhiMax);

    const bool isVertexOnHitU(this->FillHistogram(pVertex, TPC_VIEW_U, m_inputClusterListNameU, maxHitVertexDisplacement, hitDeweightingPower, histogramU));
    const bool isVertexOnHitV(this->FillHistogram(pVertex, TPC_VIEW_V, m_inputClusterListNameV, maxHitVertexDisplacement, hitDeweightingPower, histogramV));
    const bool isVertexOnHitW(this->FillHistogram(pVertex, TPC_VIEW_W, m_inputClusterListNameW, maxHitVertexDisplacement, hitDeweightingPower, histogramW));

    if (!isVertexOnHitU || !isVertexOnHitV || !isVertexOnHitW)
        throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);

    return this->GetFigureOfMerit(histogramU, histogramV, histogramW);
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

bool VertexSelectionAlgorithm::FillHistogram(const Vertex *const pVertex, const HitType hitType, const std::string &clusterListName,
    const float maxHitVertexDisplacement, const float hitDeweightingPower, Histogram &histogram) const
{
    bool isVertexOnHit(false);
    const ClusterList *pClusterList = NULL;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));

    for (ClusterList::const_iterator cIter = pClusterList->begin(), cIterEnd = pClusterList->end(); cIter != cIterEnd; ++cIter)
    {
        const Cluster *const pCluster(*cIter);
        const HitType clusterHitType(LArClusterHelper::GetClusterHitType(pCluster));

        if (hitType != clusterHitType)
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), hitType));
        isVertexOnHit |= this->FillHistogram(vertexPosition2D, pCluster, maxHitVertexDisplacement, hitDeweightingPower, histogram);
    }

    return isVertexOnHit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexSelectionAlgorithm::FillHistogram(const CartesianVector &vertexPosition2D, const Cluster *const pCluster,
    const float maxHitVertexDisplacement, const float hitDeweightingPower, Histogram &histogram) const
{
    bool isVertexOnHit(false);
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hIter = iter->second->begin(), hIterEnd = iter->second->end(); hIter != hIterEnd; ++hIter)
        {
            const CaloHit *const pCaloHit(*hIter);

            const CartesianVector displacement(pCaloHit->GetPositionVector() - vertexPosition2D);
            const float magnitude(displacement.GetMagnitude());

            if (magnitude < m_maxOnHitDisplacement)
                isVertexOnHit = true;

            if (magnitude > maxHitVertexDisplacement)
                continue;

            const float phi(std::atan2(displacement.GetZ(), displacement.GetX()));
            histogram.Fill(phi, std::pow(magnitude, hitDeweightingPower));
        }
    }

    return isVertexOnHit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionAlgorithm::SelectTopScoreVertices(const VertexScoreList &vertexScoreList, VertexScoreList &selectedVertexScoreList) const
{
    // ATTN Assumes sorted vertex score list
    unsigned int nVerticesConsidered(0);

    for (VertexScoreList::const_iterator iter = vertexScoreList.begin(), iterEnd = vertexScoreList.end(); iter != iterEnd; ++iter)
    {
        if (++nVerticesConsidered > m_maxTopScoreCandidates)
            break;

        if (selectedVertexScoreList.size() >= m_maxTopScoreSelections)
            break;

        if (!selectedVertexScoreList.empty() && !this->AcceptVertexLocation(iter->GetVertex(), selectedVertexScoreList))
            continue;

        if (!selectedVertexScoreList.empty() && !this->AcceptVertexScore(iter->GetScore(), selectedVertexScoreList))
            continue;

        selectedVertexScoreList.push_back(*iter);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexSelectionAlgorithm::AcceptVertexLocation(const Vertex *const pVertex, const VertexScoreList &vertexScoreList) const
{
    const CartesianVector position(pVertex->GetPosition());

    for (VertexScoreList::const_iterator iter = vertexScoreList.begin(), iterEnd = vertexScoreList.end(); iter != iterEnd; ++iter)
    {
        const float displacement3D((position - iter->GetVertex()->GetPosition()).GetMagnitude());

        if (displacement3D < m_minCandidateDisplacement)
            return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexSelectionAlgorithm::AcceptVertexScore(const float score, const VertexScoreList &vertexScoreList) const
{
    for (VertexScoreList::const_iterator iter = vertexScoreList.begin(), iterEnd = vertexScoreList.end(); iter != iterEnd; ++iter)
    {
        if (score < m_minCandidateScoreFraction * iter->GetScore())
            return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionAlgorithm::SelectFinalVertex(const VertexScoreList &vertexScoreList, Vertex *&pFinalVertex) const
{
    pFinalVertex = NULL;

    for (VertexScoreList::const_iterator iter = vertexScoreList.begin(), iterEnd = vertexScoreList.end(); iter != iterEnd; ++iter)
    {
        Vertex *const pVertex(iter->GetVertex());

        if (NULL == pFinalVertex)
        {
            pFinalVertex = pVertex;
            continue;
        }

        unsigned int nWinsCurrent(0), nWinsNew(0);

        for (FloatVector::const_iterator dIter = m_hitVertexDisplacementList.begin(), dIterEnd = m_hitVertexDisplacementList.end(); dIter != dIterEnd; ++dIter)
        {
            for (FloatVector::const_iterator pIter = m_hitDeweightingPowerList.begin(), pIterEnd = m_hitDeweightingPowerList.end(); pIter != pIterEnd; ++pIter)
            {
                if (this->GetFigureOfMerit(pFinalVertex, *dIter, *pIter) > this->GetFigureOfMerit(pVertex, *dIter, *pIter))
                {
                    ++nWinsCurrent;
                }
                else
                {
                    ++nWinsNew;
                }
            }
        }

        if (nWinsNew > nWinsCurrent)
            pFinalVertex = pVertex;
    }
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
        "MaxOnHitDisplacement", m_maxOnHitDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "HitDeweightingPower", m_hitDeweightingPower));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTopScoreCandidates", m_maxTopScoreCandidates));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTopScoreSelections", m_maxTopScoreSelections));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCandidateDisplacement", m_minCandidateDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCandidateScoreFraction", m_minCandidateScoreFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "HitVertexDisplacementList", m_hitVertexDisplacementList));

    if (m_hitVertexDisplacementList.empty())
    {
        const float values[4] = {5.f, 10.f, 20.f, std::numeric_limits<float>::max()};
        m_hitVertexDisplacementList.insert(m_hitVertexDisplacementList.begin(), values, values + 4);
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "HitDeweightingPowerList", m_hitDeweightingPowerList));

    if (m_hitDeweightingPowerList.empty())
    {
        const float values[4] = {0.f, -0.25f, -0.5f, -0.75f};
        m_hitDeweightingPowerList.insert(m_hitDeweightingPowerList.begin(), values, values + 4);
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
