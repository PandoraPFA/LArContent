/**
 *  @file   larpandoracontent/LArVertex/VertexSelectionBaseAlgorithm.cc
 *
 *  @brief  Implementation of the vertex selection base algorithm class.
 *
 *  $Log: $
 */
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

#include "larpandoracontent/LArVertex/VertexSelectionBaseAlgorithm.h"

using namespace pandora;

namespace lar_content
{

VertexSelectionBaseAlgorithm::VertexSelectionBaseAlgorithm() :
    m_inputVertexListName(""),
    m_replaceCurrentVertexList(true),
    m_beamMode(true),
    m_nDecayLengthsInZSpan(2.f),
    m_selectSingleVertex(true),
    m_maxTopScoreSelections(3),
    m_maxOnHitDisplacement(1.f),
    m_minCandidateDisplacement(2.f),
    m_minCandidateScoreFraction(0.5f),
    m_useDetectorGaps(true),
    m_gapTolerance(0.f),
    m_isEmptyViewAcceptable(true),
    m_minVertexAcceptableViews(3)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionBaseAlgorithm::FilterVertexList(const VertexList *const pInputVertexList, HitKDTree2D &kdTreeU, HitKDTree2D &kdTreeV,
    HitKDTree2D &kdTreeW, VertexVector &filteredVertices) const
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
            filteredVertices.push_back(pVertex);
    }

    std::sort(filteredVertices.begin(), filteredVertices.end(), SortByVertexZPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionBaseAlgorithm::GetBeamConstants(const VertexVector &vertexVector, BeamConstants &beamConstants) const
{
    if (!m_beamMode)
        return;

    if (vertexVector.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    float minZCoordinate(std::numeric_limits<float>::max()), maxZCoordinate(-std::numeric_limits<float>::max());

    for (const Vertex *const pVertex : vertexVector)
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

void VertexSelectionBaseAlgorithm::GetClusterLists(
    const StringVector &inputClusterListNames, ClusterList &clusterListU, ClusterList &clusterListV, ClusterList &clusterListW) const
{
    for (const std::string &clusterListName : inputClusterListNames)
    {
        const ClusterList *pClusterList(NULL);
        PANDORA_THROW_RESULT_IF_AND_IF(
            STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));

        if (!pClusterList || pClusterList->empty())
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "EnergyKickVertexSelectionAlgorithm: unable to find cluster list " << clusterListName << std::endl;

            continue;
        }

        const HitType hitType(LArClusterHelper::GetClusterHitType(*(pClusterList->begin())));

        if ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        ClusterList &clusterList((TPC_VIEW_U == hitType) ? clusterListU : (TPC_VIEW_V == hitType) ? clusterListV : clusterListW);
        clusterList.insert(clusterList.end(), pClusterList->begin(), pClusterList->end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionBaseAlgorithm::CalculateClusterSlidingFits(const ClusterList &inputClusterList, const unsigned int minClusterCaloHits,
    const unsigned int slidingFitWindow, SlidingFitDataList &slidingFitDataList) const
{
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));

    ClusterVector sortedClusters(inputClusterList.begin(), inputClusterList.end());
    std::sort(sortedClusters.begin(), sortedClusters.end(), LArClusterHelper::SortByNHits);

    for (const Cluster *const pCluster : sortedClusters)
    {
        if (pCluster->GetNCaloHits() < minClusterCaloHits)
            continue;

        // Make sure the window size is such that there are not more layers than hits (following TwoDSlidingLinearFit calculation).
        const unsigned int newSlidingFitWindow(
            std::min(static_cast<int>(pCluster->GetNCaloHits()), static_cast<int>(slidingFitPitch * slidingFitWindow)));
        slidingFitDataList.emplace_back(pCluster, newSlidingFitWindow, slidingFitPitch);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexSelectionBaseAlgorithm::Run()
{
    const VertexList *pInputVertexList(NULL);

    if (m_inputVertexListName == "")
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pInputVertexList));
    }
    else
    {
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_inputVertexListName, pInputVertexList));
    }

    if (!pInputVertexList || pInputVertexList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "VertexSelectionBaseAlgorithm: unable to find current vertex list " << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    HitKDTree2D kdTreeU, kdTreeV, kdTreeW;
    this->InitializeKDTrees(kdTreeU, kdTreeV, kdTreeW);

    VertexVector filteredVertices;
    this->FilterVertexList(pInputVertexList, kdTreeU, kdTreeV, kdTreeW, filteredVertices);

    if (filteredVertices.empty())
        return STATUS_CODE_SUCCESS;

    BeamConstants beamConstants;
    this->GetBeamConstants(filteredVertices, beamConstants);

    VertexScoreList vertexScoreList;
    this->GetVertexScoreList(filteredVertices, beamConstants, kdTreeU, kdTreeV, kdTreeW, vertexScoreList);

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

void VertexSelectionBaseAlgorithm::InitializeKDTrees(HitKDTree2D &kdTreeU, HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW) const
{
    for (const std::string &caloHitListName : m_inputCaloHitListNames)
    {
        const CaloHitList *pCaloHitList = NULL;
        PANDORA_THROW_RESULT_IF_AND_IF(
            STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, caloHitListName, pCaloHitList));

        if (!pCaloHitList || pCaloHitList->empty())
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "VertexSelectionBaseAlgorithm: unable to find calo hit list " << caloHitListName << std::endl;

            continue;
        }

        const HitType hitType((*(pCaloHitList->begin()))->GetHitType());

        if ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        HitKDTree2D &kdTree((TPC_VIEW_U == hitType) ? kdTreeU : (TPC_VIEW_V == hitType) ? kdTreeV : kdTreeW);

        if (!kdTree.empty())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        HitKDNode2DList hitKDNode2DList;
        KDTreeBox hitsBoundingRegion2D(fill_and_bound_2d_kd_tree(*pCaloHitList, hitKDNode2DList));
        kdTree.build(hitKDNode2DList, hitsBoundingRegion2D);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexSelectionBaseAlgorithm::IsVertexOnHit(const Vertex *const pVertex, const HitType hitType, HitKDTree2D &kdTree) const
{
    const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), hitType));
    KDTreeBox searchRegionHits = build_2d_kd_search_region(vertexPosition2D, m_maxOnHitDisplacement, m_maxOnHitDisplacement);

    HitKDNode2DList found;
    kdTree.search(searchRegionHits, found);

    return (!found.empty());
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexSelectionBaseAlgorithm::IsVertexInGap(const Vertex *const pVertex, const HitType hitType) const
{
    if (!m_useDetectorGaps)
        return false;

    return LArGeometryHelper::IsInGap3D(this->GetPandora(), pVertex->GetPosition(), hitType, m_gapTolerance);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VertexSelectionBaseAlgorithm::GetVertexEnergy(const Vertex *const pVertex, const KDTreeMap &kdTreeMap) const
{
    float totalEnergy(0.f);

    if (!this->IsVertexInGap(pVertex, TPC_VIEW_U))
        totalEnergy += this->VertexHitEnergy(pVertex, TPC_VIEW_U, kdTreeMap.at(TPC_VIEW_U));

    if (!this->IsVertexInGap(pVertex, TPC_VIEW_V))
        totalEnergy += this->VertexHitEnergy(pVertex, TPC_VIEW_V, kdTreeMap.at(TPC_VIEW_V));

    if (!this->IsVertexInGap(pVertex, TPC_VIEW_W))
        totalEnergy += this->VertexHitEnergy(pVertex, TPC_VIEW_W, kdTreeMap.at(TPC_VIEW_W));

    return totalEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VertexSelectionBaseAlgorithm::VertexHitEnergy(const Vertex *const pVertex, const HitType hitType, HitKDTree2D &kdTree) const
{
    const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), hitType));
    KDTreeBox searchRegionHits = build_2d_kd_search_region(vertexPosition2D, m_maxOnHitDisplacement, m_maxOnHitDisplacement);

    HitKDNode2DList foundHits;
    kdTree.search(searchRegionHits, foundHits);

    float dr(std::numeric_limits<float>::max());
    float energy(0);

    for (auto hit : foundHits)
    {
        const float diff = (vertexPosition2D - hit.data->GetPositionVector()).GetMagnitude();
        if (diff < dr)
        {
            dr = diff;
            energy = hit.data->GetElectromagneticEnergy();
        }
    }
    return energy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexSelectionBaseAlgorithm::SelectTopScoreVertices(VertexScoreList &vertexScoreList, VertexList &selectedVertexList) const
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

        selectedVertexList.push_back(vertexScore.GetVertex());

        if (m_selectSingleVertex)
            return;

        if (vertexScore.GetScore() > bestScore)
            bestScore = vertexScore.GetScore();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool VertexSelectionBaseAlgorithm::AcceptVertexLocation(const Vertex *const pVertex, const VertexList &selectedVertexList) const
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

bool VertexSelectionBaseAlgorithm::SortByVertexZPosition(const pandora::Vertex *const pLhs, const pandora::Vertex *const pRhs)
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

VertexSelectionBaseAlgorithm::SlidingFitData::SlidingFitData(const pandora::Cluster *const pCluster, const int slidingFitWindow, const float slidingFitPitch) :
    m_minLayerDirection(0.f, 0.f, 0.f),
    m_maxLayerDirection(0.f, 0.f, 0.f),
    m_minLayerPosition(0.f, 0.f, 0.f),
    m_maxLayerPosition(0.f, 0.f, 0.f),
    m_pCluster(pCluster)
{
    const TwoDSlidingFitResult slidingFitResult(pCluster, slidingFitWindow, slidingFitPitch);
    m_minLayerDirection = slidingFitResult.GetGlobalMinLayerDirection();
    m_maxLayerDirection = slidingFitResult.GetGlobalMaxLayerDirection();
    m_minLayerPosition = slidingFitResult.GetGlobalMinLayerPosition();
    m_maxLayerPosition = slidingFitResult.GetGlobalMaxLayerPosition();
}

//------------------------------------------------------------------------------------------------------------------------------------------

VertexSelectionBaseAlgorithm::ShowerCluster::ShowerCluster(const pandora::ClusterList &clusterList, const int slidingFitWindow, const float slidingFitPitch) :
    m_clusterList(clusterList),
    m_coordinateVector(this->GetClusterListCoordinateVector(clusterList)),
    m_twoDSlidingFitResult(&m_coordinateVector, slidingFitWindow, slidingFitPitch)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::CartesianPointVector VertexSelectionBaseAlgorithm::ShowerCluster::GetClusterListCoordinateVector(const pandora::ClusterList &clusterList) const
{
    CartesianPointVector coordinateVector;

    for (const Cluster *const pCluster : clusterList)
    {
        CartesianPointVector clusterCoordinateVector;
        LArClusterHelper::GetCoordinateVector(pCluster, clusterCoordinateVector);

        coordinateVector.insert(coordinateVector.end(), std::make_move_iterator(clusterCoordinateVector.begin()),
            std::make_move_iterator(clusterCoordinateVector.end()));
    }

    return coordinateVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexSelectionBaseAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputCaloHitListNames", m_inputCaloHitListNames));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputVertexListName", m_inputVertexListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputVertexListName", m_outputVertexListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ReplaceCurrentVertexList", m_replaceCurrentVertexList));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "BeamMode", m_beamMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NDecayLengthsInZSpan", m_nDecayLengthsInZSpan));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SelectSingleVertex", m_selectSingleVertex));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxTopScoreSelections", m_maxTopScoreSelections));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxOnHitDisplacement", m_maxOnHitDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinCandidateDisplacement", m_minCandidateDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinCandidateScoreFraction", m_minCandidateScoreFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UseDetectorGaps", m_useDetectorGaps));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "GapTolerance", m_gapTolerance));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "IsEmptyViewAcceptable", m_isEmptyViewAcceptable));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinVertexAcceptableViews", m_minVertexAcceptableViews));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
