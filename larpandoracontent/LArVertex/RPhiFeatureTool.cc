/**
 *  @file   larpandoracontent/LArVertex/RPhiFeatureTool.cc
 *
 *  @brief  Implementation of the r/phi feature tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArVertex/RPhiFeatureTool.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

using namespace pandora;

namespace lar_content
{

RPhiFeatureTool::RPhiFeatureTool() :
    m_fastScoreCheck(true),
    m_fastScoreOnly(false),
    m_fullScore(false),
    m_kernelEstimateSigma(0.048f),
    m_kappa(0.42f),
    m_maxHitVertexDisplacement1D(100.f),
    m_minFastScoreFraction(0.8f),
    m_fastHistogramNPhiBins(200),
    m_fastHistogramPhiMin(-1.1f * M_PI),
    m_fastHistogramPhiMax(+1.1f * M_PI),
    m_enableFolding(true)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void RPhiFeatureTool::Run(SupportVectorMachine::DoubleVector &featureVector, const VertexSelectionBaseAlgorithm *const pAlgorithm, const Vertex * const pVertex,
    const VertexSelectionBaseAlgorithm::SlidingFitDataListMap &, const VertexSelectionBaseAlgorithm::ClusterListMap &,
    const VertexSelectionBaseAlgorithm::KDTreeMap &kdTreeMap, const VertexSelectionBaseAlgorithm::ShowerClusterListMap &,
    const float beamDeweightingScore, float &bestFastScore)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    KernelEstimate kernelEstimateU(m_kernelEstimateSigma);
    KernelEstimate kernelEstimateV(m_kernelEstimateSigma);
    KernelEstimate kernelEstimateW(m_kernelEstimateSigma);

    this->FillKernelEstimate(pVertex, TPC_VIEW_U, kdTreeMap.at(TPC_VIEW_U), kernelEstimateU);
    this->FillKernelEstimate(pVertex, TPC_VIEW_V, kdTreeMap.at(TPC_VIEW_V), kernelEstimateV);
    this->FillKernelEstimate(pVertex, TPC_VIEW_W, kdTreeMap.at(TPC_VIEW_W), kernelEstimateW);

    const float expBeamDeweightingScore = std::exp(beamDeweightingScore);

    if (m_fastScoreCheck || m_fastScoreOnly)
    {
        const float fastScore(this->GetFastScore(kernelEstimateU, kernelEstimateV, kernelEstimateW));

        if (m_fastScoreOnly)
        {
            featureVector.push_back(fastScore);
            return;
        }

        if (expBeamDeweightingScore * fastScore < m_minFastScoreFraction * bestFastScore)
        {
            featureVector.push_back(0.);
            return;
        }

        if (expBeamDeweightingScore * fastScore > bestFastScore)
            bestFastScore = expBeamDeweightingScore * fastScore;
    }

    featureVector.push_back(m_fullScore ? this->GetFullScore(kernelEstimateU, kernelEstimateV, kernelEstimateW) :
        this->GetMidwayScore(kernelEstimateU, kernelEstimateV, kernelEstimateW));
}

//------------------------------------------------------------------------------------------------------------------------------------------

float RPhiFeatureTool::GetFastScore(const KernelEstimate &kernelEstimateU, const KernelEstimate &kernelEstimateV,
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

    // Is the below correct?
    histogramU.Scale(1.f/histogramU.GetCumulativeSum());
    histogramV.Scale(1.f/histogramV.GetCumulativeSum());
    histogramW.Scale(1.f/histogramW.GetCumulativeSum());

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

float RPhiFeatureTool::GetMidwayScore(const KernelEstimate &kernelEstimateU, const KernelEstimate &kernelEstimateV,
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

float RPhiFeatureTool::GetFullScore(const KernelEstimate &kernelEstimateU, const KernelEstimate &kernelEstimateV,
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

void RPhiFeatureTool::FillKernelEstimate(const Vertex *const pVertex, const HitType hitType, VertexSelectionBaseAlgorithm::HitKDTree2D &kdTree,
    KernelEstimate &kernelEstimate) const
{
    const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), hitType));
    KDTreeBox searchRegionHits = build_2d_kd_search_region(vertexPosition2D, m_maxHitVertexDisplacement1D, m_maxHitVertexDisplacement1D);

    VertexSelectionBaseAlgorithm::HitKDNode2DList found;
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

float RPhiFeatureTool::atan2Fast(const float y, const float x) const
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
//------------------------------------------------------------------------------------------------------------------------------------------

float RPhiFeatureTool::KernelEstimate::Sample(const float x) const
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

void RPhiFeatureTool::KernelEstimate::AddContribution(const float x, const float weight)
{
    m_contributionList.insert(ContributionList::value_type(x, weight));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode RPhiFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FastScoreCheck", m_fastScoreCheck));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FastScoreOnly", m_fastScoreOnly));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FullScore", m_fullScore));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "KernelEstimateSigma", m_kernelEstimateSigma));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "Kappa", m_kappa));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxHitVertexDisplacement1D", m_maxHitVertexDisplacement1D));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinFastScoreFraction", m_minFastScoreFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FastHistogramNPhiBins", m_fastHistogramNPhiBins));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FastHistogramPhiMin", m_fastHistogramPhiMin));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FastHistogramPhiMax", m_fastHistogramPhiMax));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "EnableFolding", m_enableFolding));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
