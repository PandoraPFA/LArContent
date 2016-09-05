/**
 *  @file   larpandoracontent/LArVertex/HitAngleVertexSelectionAlgorithm.cc
 * 
 *  @brief  Implementation of the hit angle vertex selection algorithm class.
 * 
 *  $Log: $
 */
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

#include "larpandoracontent/LArVertex/HitAngleVertexSelectionAlgorithm.h"

using namespace pandora;

namespace lar_content
{

HitAngleVertexSelectionAlgorithm::HitAngleVertexSelectionAlgorithm() :
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

void HitAngleVertexSelectionAlgorithm::GetVertexScoreList(const VertexVector &vertexVector, const BeamConstants &beamConstants, HitKDTree2D &kdTreeU,
    HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW, VertexScoreList &vertexScoreList) const
{
    float bestFastScore(0.f);

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

            if (fastScore > bestFastScore)
                bestFastScore = fastScore;
        }

        const float finalScore(multiplier * (m_fullScore ? this->GetFullScore(kernelEstimateU, kernelEstimateV, kernelEstimateW) :
            this->GetMidwayScore(kernelEstimateU, kernelEstimateV, kernelEstimateW)));

        vertexScoreList.push_back(VertexScore(pVertex, finalScore));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float HitAngleVertexSelectionAlgorithm::GetFastScore(const KernelEstimate &kernelEstimateU, const KernelEstimate &kernelEstimateV,
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

float HitAngleVertexSelectionAlgorithm::GetMidwayScore(const KernelEstimate &kernelEstimateU, const KernelEstimate &kernelEstimateV,
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

float HitAngleVertexSelectionAlgorithm::GetFullScore(const KernelEstimate &kernelEstimateU, const KernelEstimate &kernelEstimateV,
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

void HitAngleVertexSelectionAlgorithm::FillKernelEstimate(const Vertex *const pVertex, const HitType hitType, HitKDTree2D &kdTree, KernelEstimate &kernelEstimate) const
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

float HitAngleVertexSelectionAlgorithm::atan2Fast(const float y, const float x) const
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

float HitAngleVertexSelectionAlgorithm::KernelEstimate::Sample(const float x) const
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

void HitAngleVertexSelectionAlgorithm::KernelEstimate::AddContribution(const float x, const float weight)
{
    m_contributionList.insert(ContributionList::value_type(x, weight));
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HitAngleVertexSelectionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
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

    return VertexSelectionBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
