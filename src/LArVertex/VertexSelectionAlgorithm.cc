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
    m_nDecayLengthsInZSpan(2.f),
    m_selectSingleVertex(true),
    m_maxTopScoreSelections(3),
    m_kernelEstimateSigma(0.026f),
    m_maxOnHitDisplacement(1.f),
    m_maxHitVertexDisplacement1D(100.f),
    m_minCandidateDisplacement(2.f),
    m_minCandidateScoreFraction(0.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexSelectionAlgorithm::Run()
{
    const VertexList *pInputVertexList(NULL);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pInputVertexList));

    HitKDTree2D kdTreeU, kdTreeV, kdTreeW;
    this->InitializeKDTrees(kdTreeU, kdTreeV, kdTreeW);

    VertexList filteredVertexList;
    this->FilterVertexList(pInputVertexList, kdTreeU, kdTreeV, kdTreeW, filteredVertexList);

    BeamConstants beamConstants;
    this->GetBeamConstants(filteredVertexList, beamConstants);

    VertexScoreList vertexScoreList;
    this->GetVertexScoreList(filteredVertexList, beamConstants, kdTreeU, kdTreeV, kdTreeW, vertexScoreList);

    VertexList selectedVertexList;
    this->SelectTopScoreVertices(vertexScoreList, selectedVertexList);

    if (!selectedVertexList.empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_outputVertexListName, selectedVertexList));

        if (m_replaceCurrentVertexList)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_outputVertexListName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionAlgorithm::InitializeKDTrees(HitKDTree2D &kdTreeU, HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW) const
{
    this->InitializeKDTree(m_inputCaloHitListNameU, kdTreeU);
    this->InitializeKDTree(m_inputCaloHitListNameV, kdTreeV);
    this->InitializeKDTree(m_inputCaloHitListNameW, kdTreeW);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionAlgorithm::InitializeKDTree(const std::string &caloHitListName, HitKDTree2D &kdTree) const
{
    const CaloHitList *pCaloHitList = NULL;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, caloHitListName, pCaloHitList));

    if (!pCaloHitList || pCaloHitList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "VertexSelectionAlgorithm: unable to find calo hit list " << caloHitListName << std::endl;

        return;
    }

    HitKDNode2DList hitKDNode2DList;
    KDTreeBox hitsBoundingRegion2D = fill_and_bound_2d_kd_tree(this, *pCaloHitList, hitKDNode2DList, true);
    kdTree.build(hitKDNode2DList, hitsBoundingRegion2D);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionAlgorithm::FilterVertexList(const VertexList *const pInputVertexList, HitKDTree2D &kdTreeU, HitKDTree2D &kdTreeV,
    HitKDTree2D &kdTreeW, VertexList &filteredVertexList) const
{
    for (const Vertex *const pVertex : *pInputVertexList)
    {
        if ((!kdTreeU.empty() && !this->IsVertexOnHit(pVertex, TPC_VIEW_U, kdTreeU)) ||
            (!kdTreeV.empty() && !this->IsVertexOnHit(pVertex, TPC_VIEW_V, kdTreeV)) ||
            (!kdTreeW.empty() && !this->IsVertexOnHit(pVertex, TPC_VIEW_W, kdTreeW)))
        {
            continue;
        }

        (void) filteredVertexList.insert(pVertex);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionAlgorithm::GetBeamConstants(const VertexList &vertexList, BeamConstants &beamConstants) const
{
    if (!m_beamMode)
        return;

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

void VertexSelectionAlgorithm::GetVertexScoreList(const VertexList &vertexList, const BeamConstants &beamConstants, HitKDTree2D &kdTreeU,
    HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW, VertexScoreList &vertexScoreList) const
{
    for (const Vertex *const pVertex : vertexList)
    {
        KernelEstimate kernelEstimateU(m_kernelEstimateSigma);
        KernelEstimate kernelEstimateV(m_kernelEstimateSigma);
        KernelEstimate kernelEstimateW(m_kernelEstimateSigma);

        this->FillKernelEstimate(pVertex, TPC_VIEW_U, kdTreeU, kernelEstimateU);
        this->FillKernelEstimate(pVertex, TPC_VIEW_V, kdTreeV, kernelEstimateV);
        this->FillKernelEstimate(pVertex, TPC_VIEW_W, kdTreeW, kernelEstimateW);

        // TODO fast score, best score, etc.

        float figureOfMerit(this->GetFigureOfMerit(kernelEstimateU, kernelEstimateV, kernelEstimateW));

        if (m_beamMode)
            figureOfMerit *= std::exp(-(pVertex->GetPosition().GetZ() - beamConstants.GetMinZCoordinate()) * beamConstants.GetDecayConstant());

        vertexScoreList.push_back(VertexScore(pVertex, figureOfMerit));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VertexSelectionAlgorithm::GetFigureOfMerit(const KernelEstimate &kernelEstimateU, const KernelEstimate &kernelEstimateV,
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

bool VertexSelectionAlgorithm::IsVertexOnHit(const Vertex *const pVertex, const HitType hitType, HitKDTree2D &kdTree) const
{
    const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), hitType));
    KDTreeBox searchRegionHits = build_2d_kd_search_region(vertexPosition2D, m_maxOnHitDisplacement, m_maxOnHitDisplacement);

    HitKDNode2DList found;
    kdTree.search(searchRegionHits, found);

    return (!found.empty());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionAlgorithm::FillKernelEstimate(const Vertex *const pVertex, const HitType hitType, HitKDTree2D &kdTree, KernelEstimate &kernelEstimate) const
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

        const float phi(this->atan2Fast(displacement.GetZ(), displacement.GetX()));
        const float weight(1.f / std::sqrt(magnitude));
        kernelEstimate.AddContribution(phi, weight);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionAlgorithm::SelectTopScoreVertices(VertexScoreList &vertexScoreList, VertexList &selectedVertexList) const
{
    float bestScore(0.f);
    std::sort(vertexScoreList.begin(), vertexScoreList.end());

    for (const VertexScore &vertexScore : vertexScoreList)
    {
        if (selectedVertexList.size() >= m_maxTopScoreSelections)
            break;

        if (!selectedVertexList.empty() && !this->AcceptVertexLocation(vertexScore.GetVertex(), selectedVertexList))
            continue;

        if (!selectedVertexList.empty() && (vertexScore.GetScore() < m_minCandidateScoreFraction * bestScore))
            continue;

        selectedVertexList.insert(vertexScore.GetVertex());

        if (m_selectSingleVertex)
            return;

        if (vertexScore.GetScore() > bestScore)
            bestScore = vertexScore.GetScore();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexSelectionAlgorithm::AcceptVertexLocation(const Vertex *const pVertex, const VertexList &selectedVertexList) const
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

float VertexSelectionAlgorithm::atan2Fast(const float y, const float x) const
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

float VertexSelectionAlgorithm::KernelEstimate::Sample(const float x) const
{
    const float bandwidth(this->GetBandwidth());
    const float weightSum(this->GetWeightSum());

    if ((bandwidth < std::numeric_limits<float>::epsilon()) || (weightSum < std::numeric_limits<float>::epsilon()))
        return 0.f;

    const ContributionList &contributionList(this->GetContributionList());
    ContributionList::const_iterator lowerIter(contributionList.lower_bound(x - 3.f * bandwidth));
    ContributionList::const_iterator upperIter(contributionList.upper_bound(x + 3.f * bandwidth));

    float sample(0.f);
    const float gaussConstant(1.f / std::sqrt(2.f * M_PI * bandwidth * bandwidth));

    for (ContributionList::const_iterator iter = lowerIter; iter != upperIter; ++iter)
    {
        const float deltaSigma((x - iter->first) / bandwidth);
        const float gaussian(gaussConstant * std::exp(-0.5f * deltaSigma * deltaSigma));
        sample += iter->second * gaussian;
    }

    return (sample / weightSum);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VertexSelectionAlgorithm::KernelEstimate::GetBandwidth() const
{
    if (!m_bandWidth.IsInitialized())
        m_bandWidth = 1.06f * m_sigma * std::pow(m_weightSum, -0.2f);

    return m_bandWidth.Get();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionAlgorithm::KernelEstimate::AddContribution(const float x, const float weight)
{
    m_contributionList.insert(ContributionList::value_type(x, weight));
    m_weightSum += weight;
    m_bandWidth.Reset();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexSelectionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListNameU", m_inputCaloHitListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListNameV", m_inputCaloHitListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListNameW", m_inputCaloHitListNameW));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputVertexListName", m_outputVertexListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ReplaceCurrentVertexList", m_replaceCurrentVertexList));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "BeamMode", m_beamMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NDecayLengthsInZSpan", m_nDecayLengthsInZSpan));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SelectSingleVertex", m_selectSingleVertex));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTopScoreSelections", m_maxTopScoreSelections));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "KernelEstimateSigma", m_kernelEstimateSigma));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxOnHitDisplacement", m_maxOnHitDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxHitVertexDisplacement1D", m_maxHitVertexDisplacement1D));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCandidateDisplacement", m_minCandidateDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCandidateScoreFraction", m_minCandidateScoreFraction));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
