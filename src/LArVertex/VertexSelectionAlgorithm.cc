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

#include "LArUtility/KDTreeLinkerAlgoT.h"

#include "LArVertex/VertexSelectionAlgorithm.h"

using namespace pandora;

namespace lar_content
{

VertexSelectionAlgorithm::VertexSelectionAlgorithm() :
    m_replaceCurrentVertexList(true),
    m_beamMode(false),
    m_selectSingleVertex(true),
    m_histogramNPhiBins(200),
    m_histogramPhiMin(-1.1f * M_PI),
    m_histogramPhiMax(+1.1f * M_PI),
    m_maxOnHitDisplacement(1.f),
    m_maxHitVertexDisplacement1D(25.f),
    m_maxTopScoreCandidates(50),
    m_maxTopScoreSelections(3),
    m_maxBeamTopScoreCandidates(50),
    m_maxBeamTopScoreSelections(3),
    m_minCandidateDisplacement(2.f),
    m_minCandidateScoreFraction(0.5f),
    m_minBeamCandidateScoreFraction(0.5f),
    m_nDecayLengthsInZSpan(2.f),
    m_bestScoreMultiplier(0.75f),
    m_bestBeamScoreMultiplier(1.1f),
    m_mustUseBeamScoreMultiplier(1.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexSelectionAlgorithm::Run()
{
    const VertexList *pInputVertexList(NULL);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pInputVertexList));

    HitKDTree2D kdTreeU, kdTreeV, kdTreeW;
    this->InitializeKDTrees(kdTreeU, kdTreeV, kdTreeW);

    VertexScoreList vertexScoreList;
    float minZCoordinate(std::numeric_limits<float>::max()), maxZCoordinate(-std::numeric_limits<float>::max());

    for (const Vertex *const pVertex : *pInputVertexList)
    {
        try
        {
            float figureOfMerit(this->GetFigureOfMerit(pVertex, kdTreeU, kdTreeV, kdTreeW));
            vertexScoreList.push_back(VertexScore(pVertex, figureOfMerit));

            if (pVertex->GetPosition().GetZ() < minZCoordinate)
                minZCoordinate = pVertex->GetPosition().GetZ();

            if (pVertex->GetPosition().GetZ() > maxZCoordinate)
                maxZCoordinate = pVertex->GetPosition().GetZ();
        }
        catch (StatusCodeException &statusCodeException)
        {
            if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                throw statusCodeException;
        }
    }

    std::sort(vertexScoreList.begin(), vertexScoreList.end());
    const float zSpan(maxZCoordinate - minZCoordinate);
    const float decayConstant((zSpan < std::numeric_limits<float>::epsilon()) ? 0.f : (m_nDecayLengthsInZSpan / zSpan));

    VertexScoreList selectedVertexScoreList;
    this->SelectTopScoreVertices(vertexScoreList, selectedVertexScoreList);
    this->SelectTopScoreBeamVertices(vertexScoreList, minZCoordinate, decayConstant, selectedVertexScoreList);

    VertexList finalVertexList;
    this->SelectFinalVertices(selectedVertexScoreList, minZCoordinate, decayConstant, finalVertexList);

    if (!finalVertexList.empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_outputVertexListName, finalVertexList));

        if (m_replaceCurrentVertexList)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_outputVertexListName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionAlgorithm::InitializeKDTrees(HitKDTree2D &kdTreeU, HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW) const
{
    this->InitializeKDTree(m_inputClusterListNameU, kdTreeU);
    this->InitializeKDTree(m_inputClusterListNameV, kdTreeV);
    this->InitializeKDTree(m_inputClusterListNameW, kdTreeW);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionAlgorithm::InitializeKDTree(const std::string &clusterListName, HitKDTree2D &kdTree) const
{
    // TODO Decide whether to find vertices using either i) all hits, or ii) all hits in existing clusters, as below.
    const ClusterList *pClusterList = NULL;
    const StatusCode statusCode(PandoraContentApi::GetList(*this, clusterListName, pClusterList));

    if (STATUS_CODE_NOT_INITIALIZED == statusCode)
    {
        std::cout << "VertexSelectionAlgorithm: cluster list not found " << clusterListName << std::endl;
        return;
    }

    if (STATUS_CODE_SUCCESS != statusCode)
        throw StatusCodeException(statusCode);

    CaloHitList caloHitList;

    for (const Cluster *const pCluster : *pClusterList)
    {
        pCluster->GetOrderedCaloHitList().GetCaloHitList(caloHitList);
    }

    HitKDNode2DList hitKDNode2DList;
    KDTreeBox hitsBoundingRegion2D = fill_and_bound_2d_kd_tree(this, caloHitList, hitKDNode2DList, true);
    kdTree.build(hitKDNode2DList, hitsBoundingRegion2D);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VertexSelectionAlgorithm::GetFigureOfMerit(const Vertex *const pVertex, HitKDTree2D &kdTreeU, HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW) const
{
    Histogram histogramU(m_histogramNPhiBins, m_histogramPhiMin, m_histogramPhiMax);
    Histogram histogramV(m_histogramNPhiBins, m_histogramPhiMin, m_histogramPhiMax);
    Histogram histogramW(m_histogramNPhiBins, m_histogramPhiMin, m_histogramPhiMax);

    const bool isVertexOnHitU(this->FillHistogram(pVertex, TPC_VIEW_U, kdTreeU, histogramU));
    const bool isVertexOnHitV(this->FillHistogram(pVertex, TPC_VIEW_V, kdTreeV, histogramV));
    const bool isVertexOnHitW(this->FillHistogram(pVertex, TPC_VIEW_W, kdTreeW, histogramW));

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

bool VertexSelectionAlgorithm::FillHistogram(const Vertex *const pVertex, const HitType hitType, HitKDTree2D &kdTree, Histogram &histogram) const
{
    bool isVertexOnHit(false);
    const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), hitType));

    // Get nearby hits from kd tree
    CaloHitList nearbyHits;
    KDTreeBox searchRegionHits = build_2d_kd_search_region(vertexPosition2D, m_maxHitVertexDisplacement1D, m_maxHitVertexDisplacement1D);

    HitKDNode2DList found;
    kdTree.search(searchRegionHits, found);

    for (const auto &hit : found)
    {
        nearbyHits.insert(hit.data);
    }

    for (const CaloHit *const pCaloHit : nearbyHits)
    {
        const CartesianVector displacement(pCaloHit->GetPositionVector() - vertexPosition2D);
        const float magnitude(displacement.GetMagnitude());

        if (magnitude < m_maxOnHitDisplacement)
            isVertexOnHit = true;

        if (magnitude < std::numeric_limits<float>::epsilon())
            continue;

        const float phi(std::atan2(displacement.GetZ(), displacement.GetX()));
        const float weight(1.f / std::sqrt(magnitude));
        histogram.Fill(phi, weight);
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

        if (!selectedVertexScoreList.empty() && !this->AcceptVertexScore(iter->GetScore(), m_minCandidateScoreFraction, selectedVertexScoreList))
            continue;

        selectedVertexScoreList.push_back(*iter);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionAlgorithm::SelectTopScoreBeamVertices(const VertexScoreList &vertexScoreList, const float minZCoordinate,
    const float decayConstant, VertexScoreList &selectedVertexScoreList) const
{
    if (!m_beamMode)
        return;

    VertexScoreList beamVertexScoreList;

    for (VertexScoreList::const_iterator iter = vertexScoreList.begin(), iterEnd = vertexScoreList.end(); iter != iterEnd; ++iter)
    {
        const float beamScore(iter->GetScore() * std::exp(-(iter->GetVertex()->GetPosition().GetZ() - minZCoordinate) * decayConstant));
        beamVertexScoreList.push_back(VertexScore(iter->GetVertex(), beamScore));
    }

    unsigned int nVerticesConsidered(0), nVerticesAdded(0);
    std::sort(beamVertexScoreList.begin(), beamVertexScoreList.end());

    for (VertexScoreList::const_iterator iter = beamVertexScoreList.begin(), iterEnd = beamVertexScoreList.end(); iter != iterEnd; ++iter)
    {
        if (++nVerticesConsidered > m_maxBeamTopScoreCandidates)
            break;

        if (nVerticesAdded >= m_maxBeamTopScoreSelections)
            break;

        if (!selectedVertexScoreList.empty() && !this->AcceptVertexLocation(iter->GetVertex(), selectedVertexScoreList))
            continue;

        const float nonBeamScore(iter->GetScore() / std::exp(-(iter->GetVertex()->GetPosition().GetZ() - minZCoordinate) * decayConstant));

        if (!selectedVertexScoreList.empty() && !this->AcceptVertexScore(nonBeamScore, m_minBeamCandidateScoreFraction, selectedVertexScoreList))
            continue;

        ++nVerticesAdded;
        selectedVertexScoreList.push_back(VertexScore(iter->GetVertex(), nonBeamScore));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexSelectionAlgorithm::AcceptVertexLocation(const Vertex *const pVertex, const VertexScoreList &vertexScoreList) const
{
    const CartesianVector position(pVertex->GetPosition());

    for (VertexScoreList::const_iterator iter = vertexScoreList.begin(), iterEnd = vertexScoreList.end(); iter != iterEnd; ++iter)
    {
        if (iter->GetVertex() == pVertex)
            return false;

        const float displacement3D((position - iter->GetVertex()->GetPosition()).GetMagnitude());

        if (displacement3D < m_minCandidateDisplacement)
            return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexSelectionAlgorithm::AcceptVertexScore(const float score, const float minScoreFraction, const VertexScoreList &vertexScoreList) const
{
    for (VertexScoreList::const_iterator iter = vertexScoreList.begin(), iterEnd = vertexScoreList.end(); iter != iterEnd; ++iter)
    {
        if (score < minScoreFraction * iter->GetScore())
            return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionAlgorithm::SelectFinalVertices(const VertexScoreList &vertexScoreList, const float minZCoordinate, const float decayConstant,
    VertexList &finalVertexList) const
{
    const Vertex *pFinalVertex(NULL);
    float bestScore(0.f), bestBeamScore(0.f);

    for (VertexScoreList::const_iterator iter = vertexScoreList.begin(), iterEnd = vertexScoreList.end(); iter != iterEnd; ++iter)
    {
        if (!m_selectSingleVertex)
        {
            finalVertexList.insert(iter->GetVertex());
            continue;
        }

        if (!m_beamMode && (iter->GetScore() < bestScore))
            continue;

        const float beamScore(iter->GetScore() * std::exp(-(iter->GetVertex()->GetPosition().GetZ() - minZCoordinate) * decayConstant));

        if (m_beamMode && (beamScore < m_bestBeamScoreMultiplier * bestBeamScore))
            continue;

        if (m_beamMode && (iter->GetScore() < m_bestScoreMultiplier * bestScore) && (beamScore < m_mustUseBeamScoreMultiplier * bestBeamScore))
            continue;

        pFinalVertex = iter->GetVertex();
        bestScore = iter->GetScore();
        bestBeamScore = beamScore;
    }

    if (m_selectSingleVertex && (NULL != pFinalVertex))
        finalVertexList.insert(pFinalVertex);
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
        "BeamMode", m_beamMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SelectSingleVertex", m_selectSingleVertex));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "HistogramNPhiBins", m_histogramNPhiBins));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "HistogramPhiMin", m_histogramPhiMin));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "HistogramPhiMax", m_histogramPhiMax));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxOnHitDisplacement", m_maxOnHitDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxHitVertexDisplacement1D", m_maxHitVertexDisplacement1D));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTopScoreCandidates", m_maxTopScoreCandidates));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTopScoreSelections", m_maxTopScoreSelections));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxBeamTopScoreCandidates", m_maxBeamTopScoreCandidates));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxBeamTopScoreSelections", m_maxBeamTopScoreSelections));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCandidateDisplacement", m_minCandidateDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCandidateScoreFraction", m_minCandidateScoreFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinBeamCandidateScoreFraction", m_minBeamCandidateScoreFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NDecayLengthsInZSpan", m_nDecayLengthsInZSpan));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "BestScoreMultiplier", m_bestScoreMultiplier));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "BestBeamScoreMultiplier", m_bestBeamScoreMultiplier));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MustUseBeamScoreMultiplier", m_mustUseBeamScoreMultiplier));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
