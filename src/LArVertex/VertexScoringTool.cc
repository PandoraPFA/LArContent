/**
 *  @file   larpandoracontent/LArVertex/VertexScoringTool.cc
 * 
 *  @brief  Implementation of the vertex selection algorithm class.
 * 
 *  $Log: $
 */
#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"

#include "LArUtility/KDTreeLinkerAlgoT.h"

#include "LArVertex/VertexScoringTool.h"
#include "LArHelpers/LArMCParticleHelper.h"

using namespace pandora;

namespace lar_content
{

VertexScoringTool::VertexScoringTool() :
    m_replaceCurrentVertexList(true),
    m_fastScoreCheck(true),
    m_fastScoreOnly(false),
    m_fullScore(false),
    m_beamMode(false),
    m_nDecayLengthsInZSpan(2.f),
    m_kappa(0.42f),
    m_selectSingleVertex(true),
    m_maxTopScoreSelections(3),
    m_kernelEstimateSigma(0.048f),
    m_minFastScoreFraction(0.8f),
    m_fastHistogramNPhiBins(200),
    m_fastHistogramPhiMin(-1.1f * M_PI),
    m_fastHistogramPhiMax(+1.1f * M_PI),
    m_maxOnHitDisplacement(1.f),
    m_maxHitVertexDisplacement1D(100.f),
    m_minCandidateDisplacement(2.f),
    m_minCandidateScoreFraction(0.5f),
    m_enableFolding(true),
    m_useDetectorGaps(true),
    m_gapTolerance(0.f),
    m_isEmptyViewAcceptable(true),
    m_minVertexAcceptableViews(3),
    m_nSelectedVerticesPerCluster(1)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexScoringTool::ScoreVertices(const Algorithm *const pAlgorithm, const VertexList* pInputVertexList, std::vector<const VertexList*> vertexListVector, VertexScoreList &vertexScoreList)
{
    HitKDTree2D kdTreeU, kdTreeV, kdTreeW;
    this->InitializeKDTrees(pAlgorithm, kdTreeU, kdTreeV, kdTreeW);

    std::cout << "Number of vertex candidates in event: " << pInputVertexList->size() << std::endl;

    VertexList filteredVertexList;
    this->FilterVertexList(pInputVertexList, kdTreeU, kdTreeV, kdTreeW, filteredVertexList);

    if (filteredVertexList.empty())
        return;

    std::cout << "Number of filtered vertex candidates in event: " << filteredVertexList.size() << std::endl;

    BeamConstants beamConstants;
    this->GetBeamConstants(filteredVertexList, beamConstants);

    float bestFastScore(0.f);
    this->GetGlobalBestFastScore(filteredVertexList, bestFastScore, beamConstants, kdTreeU, kdTreeV, kdTreeW);

    std::vector<VertexScoreList> scoredClusterCollection;

    for (const VertexList *const pVertexList : vertexListVector)
    {
        VertexList filteredCluster;
        for (const Vertex *const pVertex: (*pVertexList))
        {
            if (filteredVertexList.count(pVertex) != 0)
                filteredCluster.insert(pVertex);
        }

        VertexScoreList clusterVertexScoreList;
        this->GetVertexScoreList(filteredCluster, bestFastScore, beamConstants, kdTreeU, kdTreeV, kdTreeW, clusterVertexScoreList);

        std::sort(clusterVertexScoreList.begin(), clusterVertexScoreList.end());
        scoredClusterCollection.push_back(clusterVertexScoreList);
        
//        if (clusterVertexScoreList.size() >= m_nSelectedVerticesPerCluster)
//            vertexScoreList.insert(vertexScoreList.begin(), clusterVertexScoreList.begin(), std::next(clusterVertexScoreList.begin(), m_nSelectedVerticesPerCluster) );
//        else
//        {
//            for (VertexScore &vertexScore : clusterVertexScoreList)
//                vertexScoreList.push_back(vertexScore);
//        }

    }

    std::sort(scoredClusterCollection.begin(), scoredClusterCollection.end(), SortClustersByScore);

    unsigned int clusterCounter(0);

    for (VertexScoreList &thisVertexScoreList : scoredClusterCollection)
    {
        if (clusterCounter == 5 || clusterCounter == scoredClusterCollection.size() || thisVertexScoreList.size() == 0)
            break;

        std::sort(thisVertexScoreList.begin(), thisVertexScoreList.end());
    
        float totalClusterScore(0.f);
        int vertexCounter(0);

        for (VertexScore &thisVertexScore : thisVertexScoreList)
        {
            if (vertexCounter == 4)
                break;

            totalClusterScore += thisVertexScore.GetScore();
            vertexCounter++;
        }

        if (thisVertexScoreList.size() != 0)
        {        
            std::cout << "Total cluster score is: " << totalClusterScore << std::endl;
            std::cout << "Number of vertices in cluster: " << thisVertexScoreList.size() << std::endl;
            std::cout << "************************************" << std::endl;
        }

        vertexScoreList.push_back(*(thisVertexScoreList.begin()));

        clusterCounter++;
    }

   
    std::sort(vertexScoreList.begin(), vertexScoreList.end());
   
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexScoringTool::InitializeKDTrees(const Algorithm *const pAlgorithm, HitKDTree2D &kdTreeU, HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW) const
{
    this->InitializeKDTree(pAlgorithm, m_inputCaloHitListNameU, kdTreeU);
    this->InitializeKDTree(pAlgorithm, m_inputCaloHitListNameV, kdTreeV);
    this->InitializeKDTree(pAlgorithm, m_inputCaloHitListNameW, kdTreeW);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexScoringTool::InitializeKDTree(const Algorithm *const pAlgorithm, const std::string &caloHitListName, HitKDTree2D &kdTree) const
{
    const CaloHitList *pCaloHitList = NULL;
    PandoraContentApi::GetList(*pAlgorithm, caloHitListName, pCaloHitList);

    if (!pCaloHitList || pCaloHitList->empty())
    {
        if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
            std::cout << "VertexScoringTool: unable to find calo hit list " << caloHitListName << std::endl;

        return;
    }

    HitKDNode2DList hitKDNode2DList;
    KDTreeBox hitsBoundingRegion2D = fill_and_bound_2d_kd_tree(pAlgorithm, *pCaloHitList, hitKDNode2DList, true);
    kdTree.build(hitKDNode2DList, hitsBoundingRegion2D);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexScoringTool::FilterVertexList(const VertexList *const pInputVertexList, HitKDTree2D &kdTreeU, HitKDTree2D &kdTreeV,
    HitKDTree2D &kdTreeW, VertexList &filteredVertexList) const
{
    for (const Vertex *const pVertex : *pInputVertexList)
    {
        unsigned int nAcceptableViews(0);

        if ((m_isEmptyViewAcceptable && kdTreeU.empty()) || this->IsVertexOnHit(pVertex, TPC_VIEW_U, kdTreeU) || this->IsVertexInGap(pVertex, TPC_VIEW_U))
            ++nAcceptableViews;

        if ((m_isEmptyViewAcceptable && kdTreeV.empty()) || this->IsVertexOnHit(pVertex, TPC_VIEW_V, kdTreeV) || this->IsVertexInGap(pVertex, TPC_VIEW_V))
            ++nAcceptableViews;

        if ((m_isEmptyViewAcceptable && kdTreeW.empty()) || this->IsVertexOnHit(pVertex, TPC_VIEW_W, kdTreeW) || this->IsVertexInGap(pVertex, TPC_VIEW_W))
            ++nAcceptableViews;

        if (nAcceptableViews >= m_minVertexAcceptableViews)
            (void) filteredVertexList.insert(pVertex);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexScoringTool::GetBeamConstants(const VertexList &vertexList, BeamConstants &beamConstants) const
{
//    if (!m_beamMode)
//        return;

    if (vertexList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    float minZCoordinate(std::numeric_limits<float>::max()), maxZCoordinate(-std::numeric_limits<float>::max());

    for (const Vertex *const pVertex : vertexList)
    {
        if (pVertex->GetPosition().GetZ() < minZCoordinate)
            minZCoordinate = pVertex->GetPosition().GetZ();

        if (pVertex->GetPosition().GetZ() > maxZCoordinate)
            maxZCoordinate = pVertex->GetPosition().GetZ();
    }

    const float zSpan(maxZCoordinate - minZCoordinate);
    const float decayConstant((zSpan < std::numeric_limits<float>::epsilon()) ? 0.f : (m_nDecayLengthsInZSpan / zSpan));

    beamConstants.SetConstants(minZCoordinate, decayConstant);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexScoringTool::GetVertexScoreList(const VertexList &vertexList, const float &bestFastScore, const BeamConstants &beamConstants, HitKDTree2D &kdTreeU,
    HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW, VertexScoreList &vertexScoreList) const
{
    VertexVector vertexVector(vertexList.begin(), vertexList.end());
    std::sort(vertexVector.begin(), vertexVector.end(), SortByVertexZPosition);

    for (const Vertex *const pVertex : vertexVector)
    {
        KernelEstimate kernelEstimateU(m_kernelEstimateSigma);
        KernelEstimate kernelEstimateV(m_kernelEstimateSigma);
        KernelEstimate kernelEstimateW(m_kernelEstimateSigma);

        this->FillKernelEstimate(pVertex, TPC_VIEW_U, kdTreeU, kernelEstimateU);
        this->FillKernelEstimate(pVertex, TPC_VIEW_V, kdTreeV, kernelEstimateV);
        this->FillKernelEstimate(pVertex, TPC_VIEW_W, kdTreeW, kernelEstimateW);

        const float vertexMinZ(std::max(pVertex->GetPosition().GetZ(), beamConstants.GetMinZCoordinate()));
        const float multiplier(!m_beamMode ? 1.f : std::exp(-(vertexMinZ - beamConstants.GetMinZCoordinate()) * beamConstants.GetDecayConstant()));
        
        if (m_fastScoreCheck || m_fastScoreOnly)
        {
            const float fastScore(multiplier * this->GetFastScore(kernelEstimateU, kernelEstimateV, kernelEstimateW));

            if (m_fastScoreOnly)
            {
                vertexScoreList.push_back(VertexScore(pVertex, fastScore));
                continue;
            }

            if (fastScore < m_minFastScoreFraction * bestFastScore)
                continue;
        }

        const float finalScore(multiplier * (m_fullScore ? this->GetFullScore(kernelEstimateU, kernelEstimateV, kernelEstimateW) :
            this->GetMidwayScore(kernelEstimateU, kernelEstimateV, kernelEstimateW)));

        vertexScoreList.push_back(VertexScore(pVertex, finalScore));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexScoringTool::GetGlobalBestFastScore(const VertexList &vertexList, float &bestFastScore, const BeamConstants &beamConstants, HitKDTree2D &kdTreeU,
    HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW)
{
    VertexVector vertexVector(vertexList.begin(), vertexList.end());
    std::sort(vertexVector.begin(), vertexVector.end(), SortByVertexZPosition);

    for (const Vertex *const pVertex : vertexVector)
    {
        KernelEstimate kernelEstimateU(m_kernelEstimateSigma);
        KernelEstimate kernelEstimateV(m_kernelEstimateSigma);
        KernelEstimate kernelEstimateW(m_kernelEstimateSigma);

        this->FillKernelEstimate(pVertex, TPC_VIEW_U, kdTreeU, kernelEstimateU);
        this->FillKernelEstimate(pVertex, TPC_VIEW_V, kdTreeV, kernelEstimateV);
        this->FillKernelEstimate(pVertex, TPC_VIEW_W, kdTreeW, kernelEstimateW);

        const float vertexMinZ(std::max(pVertex->GetPosition().GetZ(), beamConstants.GetMinZCoordinate()));
        const float multiplier(!m_beamMode ? 1.f : std::exp(-(vertexMinZ - beamConstants.GetMinZCoordinate()) * beamConstants.GetDecayConstant()));
        
        const float fastScore(multiplier * this->GetFastScore(kernelEstimateU, kernelEstimateV, kernelEstimateW));

        if (fastScore > bestFastScore)
            bestFastScore = fastScore;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VertexScoringTool::GetFastScore(const KernelEstimate &kernelEstimateU, const KernelEstimate &kernelEstimateV,
    const KernelEstimate &kernelEstimateW) const
{
    Histogram histogramU(m_fastHistogramNPhiBins, m_fastHistogramPhiMin, m_fastHistogramPhiMax);
    Histogram histogramV(m_fastHistogramNPhiBins, m_fastHistogramPhiMin, m_fastHistogramPhiMax);
    Histogram histogramW(m_fastHistogramNPhiBins, m_fastHistogramPhiMin, m_fastHistogramPhiMax);

    for (const KernelEstimate::ContributionList::value_type &contribution : kernelEstimateU.GetContributionList())
        histogramU.Fill(contribution.first, contribution.second);

    for (const KernelEstimate::ContributionList::value_type &contribution : kernelEstimateV.GetContributionList())
        histogramV.Fill(contribution.first, contribution.second);

    for (const KernelEstimate::ContributionList::value_type &contribution : kernelEstimateW.GetContributionList())
        histogramW.Fill(contribution.first, contribution.second);

    // ATTN Need to renormalise histograms if ever want to directly compare fast and full scores
    float figureOfMerit(0.f);

    for (int xBin = 0; xBin < histogramU.GetNBinsX(); ++xBin)
    {
        const float binContentU(histogramU.GetBinContent(xBin));
        const float binContentV(histogramV.GetBinContent(xBin));
        const float binContentW(histogramW.GetBinContent(xBin));
        figureOfMerit += binContentU * binContentU + binContentV * binContentV + binContentW * binContentW;
    }

    return figureOfMerit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VertexScoringTool::GetMidwayScore(const KernelEstimate &kernelEstimateU, const KernelEstimate &kernelEstimateV,
    const KernelEstimate &kernelEstimateW) const
{
    Histogram histogramU(m_fastHistogramNPhiBins, m_fastHistogramPhiMin, m_fastHistogramPhiMax);
    Histogram histogramV(m_fastHistogramNPhiBins, m_fastHistogramPhiMin, m_fastHistogramPhiMax);
    Histogram histogramW(m_fastHistogramNPhiBins, m_fastHistogramPhiMin, m_fastHistogramPhiMax);

    for (const KernelEstimate::ContributionList::value_type &contribution : kernelEstimateU.GetContributionList())
        histogramU.Fill(contribution.first, contribution.second);

    for (const KernelEstimate::ContributionList::value_type &contribution : kernelEstimateV.GetContributionList())
        histogramV.Fill(contribution.first, contribution.second);

    for (const KernelEstimate::ContributionList::value_type &contribution : kernelEstimateW.GetContributionList())
        histogramW.Fill(contribution.first, contribution.second);

    float figureOfMerit(0.f);

    for (int xBin = 0; xBin < histogramU.GetNBinsX(); ++xBin)
    {
        const float binCenter(histogramU.GetXLow() + (static_cast<float>(xBin) + 0.5f) * histogramU.GetXBinWidth());
        figureOfMerit += histogramU.GetBinContent(xBin) * kernelEstimateU.Sample(binCenter);
        figureOfMerit += histogramV.GetBinContent(xBin) * kernelEstimateV.Sample(binCenter);
        figureOfMerit += histogramW.GetBinContent(xBin) * kernelEstimateW.Sample(binCenter);
    }

    return figureOfMerit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VertexScoringTool::GetFullScore(const KernelEstimate &kernelEstimateU, const KernelEstimate &kernelEstimateV,
    const KernelEstimate &kernelEstimateW) const
{
    float figureOfMerit(0.f);

    for (const KernelEstimate::ContributionList::value_type &contribution : kernelEstimateU.GetContributionList())
        figureOfMerit += contribution.second * kernelEstimateU.Sample(contribution.first);

    for (const KernelEstimate::ContributionList::value_type &contribution : kernelEstimateV.GetContributionList())
        figureOfMerit += contribution.second * kernelEstimateV.Sample(contribution.first);

    for (const KernelEstimate::ContributionList::value_type &contribution : kernelEstimateW.GetContributionList())
        figureOfMerit += contribution.second * kernelEstimateW.Sample(contribution.first);

    return figureOfMerit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexScoringTool::IsVertexOnHit(const Vertex *const pVertex, const HitType hitType, HitKDTree2D &kdTree) const
{
    const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), hitType));
    KDTreeBox searchRegionHits = build_2d_kd_search_region(vertexPosition2D, m_maxOnHitDisplacement, m_maxOnHitDisplacement);

    HitKDNode2DList found;
    kdTree.search(searchRegionHits, found);

    return (!found.empty());
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexScoringTool::IsVertexInGap(const Vertex *const pVertex, const HitType hitType) const
{
    if (!m_useDetectorGaps)
        return false;

    return LArGeometryHelper::IsInGap3D(this->GetPandora(), pVertex->GetPosition(), hitType, m_gapTolerance);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexScoringTool::FillKernelEstimate(const Vertex *const pVertex, const HitType hitType, HitKDTree2D &kdTree, KernelEstimate &kernelEstimate) const
{
    const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), hitType));
    KDTreeBox searchRegionHits = build_2d_kd_search_region(vertexPosition2D, m_maxHitVertexDisplacement1D, m_maxHitVertexDisplacement1D);

    HitKDNode2DList found;
    kdTree.search(searchRegionHits, found);

    for (const auto &hit : found)
    {
        const CartesianVector displacement(hit.data->GetPositionVector() - vertexPosition2D);
        const float magnitude(displacement.GetMagnitude());

        if (magnitude < std::numeric_limits<float>::epsilon())
            continue;

        float phi(this->atan2Fast(displacement.GetZ(), displacement.GetX()));
        float weight(1.f / (std::sqrt(magnitude + std::fabs(m_kappa))));

        if (m_enableFolding && (phi < 0.f))
        {
            phi += M_PI;
            weight *= -1.f;
        }

        kernelEstimate.AddContribution(phi, weight);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexScoringTool::AcceptVertexLocation(const Vertex *const pVertex, const VertexList &selectedVertexList) const
{
    const CartesianVector &position(pVertex->GetPosition());
    const float minCandidateDisplacementSquared(m_minCandidateDisplacement * m_minCandidateDisplacement);

    for (const Vertex *const pSelectedVertex : selectedVertexList)
    {
        if (pVertex == pSelectedVertex)
            return false;

        if ((position - pSelectedVertex->GetPosition()).GetMagnitudeSquared() < minCandidateDisplacementSquared)
            return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VertexScoringTool::atan2Fast(const float y, const float x) const
{
    const float ONE_QTR_PI(0.25f * M_PI);
    const float THR_QTR_PI(0.75f * M_PI);

    const float abs_y(std::max(std::fabs(y), std::numeric_limits<float>::epsilon()));
    const float abs_x(std::fabs(x));

    const float r((x < 0.f) ? (x + abs_y) / (abs_y + abs_x) : (abs_x - abs_y) / (abs_x + abs_y));
    const float angle(((x < 0.f) ? THR_QTR_PI : ONE_QTR_PI) + (0.1963f * r * r - 0.9817f) * r);

    return ((y < 0.f) ? -angle : angle); // negate if in quad III or IV
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexScoringTool::SortByVertexZPosition(const pandora::Vertex *const pLhs, const pandora::Vertex *const pRhs)
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

bool VertexScoringTool::SortClustersByScore(VertexScoreList &firstScoreList, VertexScoreList &secondScoreList)
{
    float firstTotalClusterScore(0.f);
    for (VertexScore &firstTempVertexScore : firstScoreList)
        firstTotalClusterScore += firstTempVertexScore.GetScore();

    float secondTotalClusterScore(0.f);
    for (VertexScore &secondTempVertexScore : secondScoreList)
        secondTotalClusterScore += secondTempVertexScore.GetScore();

    return firstTotalClusterScore > secondTotalClusterScore;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

float VertexScoringTool::KernelEstimate::Sample(const float x) const
{
    const ContributionList &contributionList(this->GetContributionList());
    ContributionList::const_iterator lowerIter(contributionList.lower_bound(x - 3.f * m_sigma));
    ContributionList::const_iterator upperIter(contributionList.upper_bound(x + 3.f * m_sigma));

    float sample(0.f);
    const float gaussConstant(1.f / std::sqrt(2.f * M_PI * m_sigma * m_sigma));

    for (ContributionList::const_iterator iter = lowerIter; iter != upperIter; ++iter)
    {
        const float deltaSigma((x - iter->first) / m_sigma);
        const float gaussian(gaussConstant * std::exp(-0.5f * deltaSigma * deltaSigma));
        sample += iter->second * gaussian;
    }

    return sample;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexScoringTool::KernelEstimate::AddContribution(const float x, const float weight)
{
    m_contributionList.insert(ContributionList::value_type(x, weight));
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexScoringTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListNameU", m_inputCaloHitListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListNameV", m_inputCaloHitListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListNameW", m_inputCaloHitListNameW));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ReplaceCurrentVertexList", m_replaceCurrentVertexList));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FastScoreCheck", m_fastScoreCheck));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FastScoreOnly", m_fastScoreOnly));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FullScore", m_fullScore));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "BeamMode", m_beamMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NDecayLengthsInZSpan", m_nDecayLengthsInZSpan));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "Kappa", m_kappa));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SelectSingleVertex", m_selectSingleVertex));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTopScoreSelections", m_maxTopScoreSelections));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "KernelEstimateSigma", m_kernelEstimateSigma));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinFastScoreFraction", m_minFastScoreFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FastHistogramNPhiBins", m_fastHistogramNPhiBins));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FastHistogramPhiMin", m_fastHistogramPhiMin));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FastHistogramPhiMax", m_fastHistogramPhiMax));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxOnHitDisplacement", m_maxOnHitDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxHitVertexDisplacement1D", m_maxHitVertexDisplacement1D));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCandidateDisplacement", m_minCandidateDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCandidateScoreFraction", m_minCandidateScoreFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "EnableFolding", m_enableFolding));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseDetectorGaps", m_useDetectorGaps));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "GapTolerance", m_gapTolerance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "IsEmptyViewAcceptable", m_isEmptyViewAcceptable));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinVertexAcceptableViews", m_minVertexAcceptableViews));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SelectedVerticesPerCluster", m_nSelectedVerticesPerCluster));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
