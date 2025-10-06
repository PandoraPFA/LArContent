/**
 *  @file   larpandoracontent/LArVertex/CandidateVertexCreationAlgorithm.cc
 *
 *  @brief  Implementation of the candidate vertex creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArVertex/CandidateVertexCreationAlgorithm.h"

#include <utility>

using namespace pandora;

namespace lar_content
{

CandidateVertexCreationAlgorithm::CandidateVertexCreationAlgorithm() :
    m_replaceCurrentVertexList(true),
    m_slidingFitWindow(20),
    m_minClusterCaloHits(5),
    m_minClusterLengthSquared(3.f * 3.f),
    m_chiSquaredCut(2.f),
    m_enableEndpointCandidates(true),
    m_maxEndpointXDiscrepancy(4.f),
    m_enableCrossingCandidates(false),
    m_nMaxCrossingCandidates(500),
    m_maxCrossingXDiscrepancy(0.5f),
    m_extrapolationNSteps(200),
    m_extrapolationStepSize(0.1f),
    m_maxCrossingSeparationSquared(2.f * 2.f),
    m_minNearbyCrossingDistanceSquared(0.5f * 0.5f),
    m_reducedCandidates(false),
    m_selectionCutFactorMax(2.f),
    m_nClustersPassingMaxCutsPar(26.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CandidateVertexCreationAlgorithm::Run()
{
    try
    {
        // INFO: See if there is already a vertex, and quit early if there is.
        //       The vertex has likely already been defined by another algorithm.
        const VertexList *pVertexList(nullptr);
        const StatusCode status{PandoraContentApi::GetCurrentList(*this, pVertexList)};
        if (status == STATUS_CODE_SUCCESS && pVertexList != nullptr && !pVertexList->empty())
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "CandidateVertexCreationAlgorithm: Vertex already defined, skipping" << std::endl;

            return STATUS_CODE_SUCCESS;
        }

        ClusterVector clusterVectorU, clusterVectorV, clusterVectorW;
        this->SelectClusters(clusterVectorU, clusterVectorV, clusterVectorW);

        std::string temporaryListName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, temporaryListName));

        if (m_enableEndpointCandidates)
        {
            this->CreateEndpointCandidates(clusterVectorU, clusterVectorV);
            this->CreateEndpointCandidates(clusterVectorU, clusterVectorW);
            this->CreateEndpointCandidates(clusterVectorV, clusterVectorW);
        }

        if (m_enableCrossingCandidates)
            this->CreateCrossingCandidates(clusterVectorU, clusterVectorV, clusterVectorW);

        if (!m_inputVertexListName.empty())
            this->AddInputVertices();

        if (!pVertexList->empty())
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_outputVertexListName));

            if (m_replaceCurrentVertexList)
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_outputVertexListName));
        }
    }
    catch (StatusCodeException &statusCodeException)
    {
        this->TidyUp();
        throw statusCodeException;
    }

    this->TidyUp();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::SelectClusters(ClusterVector &clusterVectorU, ClusterVector &clusterVectorV, ClusterVector &clusterVectorW)
{
    for (const std::string &clusterListName : m_inputClusterListNames)
    {
        const ClusterList *pClusterList(NULL);
        PANDORA_THROW_RESULT_IF_AND_IF(
            STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));

        if (!pClusterList || pClusterList->empty())
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "CandidateVertexCreationAlgorithm: unable to find cluster list " << clusterListName << std::endl;

            continue;
        }

        const HitType hitType(LArClusterHelper::GetClusterHitType(*(pClusterList->begin())));

        if ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        ClusterVector &selectedClusterVector((TPC_VIEW_U == hitType) ? clusterVectorU
                : (TPC_VIEW_V == hitType)                            ? clusterVectorV
                                                                     : clusterVectorW);

        if (!selectedClusterVector.empty())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        ClusterVector sortedClusters(pClusterList->begin(), pClusterList->end());
        std::sort(sortedClusters.begin(), sortedClusters.end(), LArClusterHelper::SortByNHits);

        unsigned int nClustersPassingMaxCuts(0);
        if (m_reducedCandidates)
        {
            for (const Cluster *const pCluster : sortedClusters)
            {
                float selectionCutFactor(1.f);

                if (pCluster->GetParticleId() == E_MINUS)
                    selectionCutFactor = m_selectionCutFactorMax;

                if (pCluster->GetNCaloHits() < m_minClusterCaloHits * selectionCutFactor)
                    continue;

                if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLengthSquared * selectionCutFactor * selectionCutFactor)
                    continue;

                nClustersPassingMaxCuts++;
            }
        }

        for (const Cluster *const pCluster : sortedClusters)
        {
            float selectionCutFactor(1.f);

            if (pCluster->GetParticleId() == E_MINUS && m_reducedCandidates)
            {
                selectionCutFactor = (m_selectionCutFactorMax + 1.f) * 0.5f +
                    (m_selectionCutFactorMax - 1.f) * 0.5f * std::tanh(static_cast<float>(nClustersPassingMaxCuts) - m_nClustersPassingMaxCutsPar);
            }

            if (pCluster->GetNCaloHits() < m_minClusterCaloHits * selectionCutFactor)
                continue;

            if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLengthSquared * selectionCutFactor * selectionCutFactor)
                continue;

            try
            {
                this->AddToSlidingFitCache(pCluster);
                selectedClusterVector.push_back(pCluster);
            }
            catch (StatusCodeException &statusCodeException)
            {
                if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                    throw statusCodeException;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::CreateEndpointCandidates(const ClusterVector &clusterVector1, const ClusterVector &clusterVector2) const
{
    for (const Cluster *const pCluster1 : clusterVector1)
    {
        const HitType hitType1(LArClusterHelper::GetClusterHitType(pCluster1));

        const TwoDSlidingFitResult &fitResult1(this->GetCachedSlidingFitResult(pCluster1));
        const CartesianVector minLayerPosition1(fitResult1.GetGlobalMinLayerPosition());
        const CartesianVector maxLayerPosition1(fitResult1.GetGlobalMaxLayerPosition());

        for (const Cluster *const pCluster2 : clusterVector2)
        {
            const HitType hitType2(LArClusterHelper::GetClusterHitType(pCluster2));

            const TwoDSlidingFitResult &fitResult2(this->GetCachedSlidingFitResult(pCluster2));
            const CartesianVector minLayerPosition2(fitResult2.GetGlobalMinLayerPosition());
            const CartesianVector maxLayerPosition2(fitResult2.GetGlobalMaxLayerPosition());

            this->CreateEndpointVertex(maxLayerPosition1, hitType1, fitResult2);
            this->CreateEndpointVertex(minLayerPosition1, hitType1, fitResult2);
            this->CreateEndpointVertex(maxLayerPosition2, hitType2, fitResult1);
            this->CreateEndpointVertex(minLayerPosition2, hitType2, fitResult1);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::CreateEndpointVertex(
    const CartesianVector &position1, const HitType hitType1, const TwoDSlidingFitResult &fitResult2) const
{
    const CartesianVector minLayerPosition2(fitResult2.GetGlobalMinLayerPosition());
    const CartesianVector maxLayerPosition2(fitResult2.GetGlobalMaxLayerPosition());

    if ((((position1.GetX() < minLayerPosition2.GetX()) && (position1.GetX() < maxLayerPosition2.GetX())) ||
            ((position1.GetX() > minLayerPosition2.GetX()) && (position1.GetX() > maxLayerPosition2.GetX()))) &&
        (std::fabs(position1.GetX() - minLayerPosition2.GetX()) > m_maxEndpointXDiscrepancy) &&
        (std::fabs(position1.GetX() - maxLayerPosition2.GetX()) > m_maxEndpointXDiscrepancy))
    {
        return;
    }

    CartesianVector position2(0.f, 0.f, 0.f);
    if (STATUS_CODE_SUCCESS != fitResult2.GetExtrapolatedPositionAtX(position1.GetX(), position2))
        return;

    const HitType hitType2(LArClusterHelper::GetClusterHitType(fitResult2.GetCluster()));

    float chiSquared(0.f);
    CartesianVector position3D(0.f, 0.f, 0.f);
    LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), hitType1, hitType2, position1, position2, position3D, chiSquared);

    if (chiSquared > m_chiSquaredCut)
        return;

    PandoraContentApi::Vertex::Parameters parameters;
    parameters.m_position = position3D;
    parameters.m_vertexLabel = VERTEX_INTERACTION;
    parameters.m_vertexType = VERTEX_3D;

    const Vertex *pVertex(NULL);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertex));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::CreateCrossingCandidates(
    const ClusterVector &clusterVectorU, const ClusterVector &clusterVectorV, const ClusterVector &clusterVectorW) const
{
    CartesianPointVector crossingsU, crossingsV, crossingsW;
    this->FindCrossingPoints(clusterVectorU, crossingsU);
    this->FindCrossingPoints(clusterVectorV, crossingsV);
    this->FindCrossingPoints(clusterVectorW, crossingsW);

    unsigned int nCrossingCandidates(0);
    this->CreateCrossingVertices(crossingsU, crossingsV, TPC_VIEW_U, TPC_VIEW_V, nCrossingCandidates);
    this->CreateCrossingVertices(crossingsU, crossingsW, TPC_VIEW_U, TPC_VIEW_W, nCrossingCandidates);
    this->CreateCrossingVertices(crossingsV, crossingsW, TPC_VIEW_V, TPC_VIEW_W, nCrossingCandidates);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::FindCrossingPoints(const ClusterVector &clusterVector, CartesianPointVector &crossingPoints) const
{
    ClusterToSpacepointsMap clusterToSpacepointsMap;

    for (const Cluster *const pCluster : clusterVector)
    {
        ClusterToSpacepointsMap::iterator mapIter(clusterToSpacepointsMap.emplace(pCluster, CartesianPointVector()).first);
        this->GetSpacepoints(pCluster, mapIter->second);
    }

    for (const Cluster *const pCluster1 : clusterVector)
    {
        for (const Cluster *const pCluster2 : clusterVector)
        {
            if (pCluster1 == pCluster2)
                continue;

            this->FindCrossingPoints(clusterToSpacepointsMap.at(pCluster1), clusterToSpacepointsMap.at(pCluster2), crossingPoints);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::GetSpacepoints(const Cluster *const pCluster, CartesianPointVector &spacepoints) const
{
    LArClusterHelper::GetCoordinateVector(pCluster, spacepoints);

    const TwoDSlidingFitResult &fitResult(this->GetCachedSlidingFitResult(pCluster));
    const float minLayerRL(fitResult.GetL(fitResult.GetMinLayer()));
    const float maxLayerRL(fitResult.GetL(fitResult.GetMaxLayer()));

    for (unsigned int iStep = 0; iStep < m_extrapolationNSteps; ++iStep)
    {
        const float deltaRL(static_cast<float>(iStep) * m_extrapolationStepSize);

        CartesianVector positionPositive(0.f, 0.f, 0.f), positionNegative(0.f, 0.f, 0.f);
        fitResult.GetExtrapolatedPosition(maxLayerRL + deltaRL, positionPositive);
        fitResult.GetExtrapolatedPosition(minLayerRL - deltaRL, positionNegative);

        spacepoints.push_back(positionPositive);
        spacepoints.push_back(positionNegative);
    }

    std::sort(spacepoints.begin(), spacepoints.end(), LArClusterHelper::SortCoordinatesByPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::FindCrossingPoints(
    const CartesianPointVector &spacepoints1, const CartesianPointVector &spacepoints2, CartesianPointVector &crossingPoints) const
{
    bool bestCrossingFound(false);
    float bestSeparationSquared(m_maxCrossingSeparationSquared);
    CartesianVector bestPosition1(0.f, 0.f, 0.f), bestPosition2(0.f, 0.f, 0.f);

    for (const CartesianVector &position1 : spacepoints1)
    {
        for (const CartesianVector &position2 : spacepoints2)
        {
            const float separationSquared((position1 - position2).GetMagnitudeSquared());

            if (separationSquared < bestSeparationSquared)
            {
                bestCrossingFound = true;
                bestSeparationSquared = separationSquared;
                bestPosition1 = position1;
                bestPosition2 = position2;
            }
        }
    }

    if (bestCrossingFound)
    {
        bool alreadyPopulated(false);

        for (const CartesianVector &existingPosition : crossingPoints)
        {
            if (((existingPosition - bestPosition1).GetMagnitudeSquared() < m_minNearbyCrossingDistanceSquared) ||
                ((existingPosition - bestPosition2).GetMagnitudeSquared() < m_minNearbyCrossingDistanceSquared))
            {
                alreadyPopulated = true;
                break;
            }
        }

        if (!alreadyPopulated)
        {
            crossingPoints.push_back(bestPosition1);
            crossingPoints.push_back(bestPosition2);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::CreateCrossingVertices(const CartesianPointVector &crossingPoints1,
    const CartesianPointVector &crossingPoints2, const HitType hitType1, const HitType hitType2, unsigned int &nCrossingCandidates) const
{

    for (const CartesianVector &position1 : crossingPoints1)
    {
        for (const CartesianVector &position2 : crossingPoints2)
        {
            if (nCrossingCandidates > m_nMaxCrossingCandidates)
                return;

            if (std::fabs(position1.GetX() - position2.GetX()) > m_maxCrossingXDiscrepancy)
                continue;

            float chiSquared(0.f);
            CartesianVector position3D(0.f, 0.f, 0.f);
            LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), hitType1, hitType2, position1, position2, position3D, chiSquared);

            if (chiSquared > m_chiSquaredCut)
                continue;

            PandoraContentApi::Vertex::Parameters parameters;
            parameters.m_position = position3D;
            parameters.m_vertexLabel = VERTEX_INTERACTION;
            parameters.m_vertexType = VERTEX_3D;

            const Vertex *pVertex(NULL);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertex));
            ++nCrossingCandidates;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::AddInputVertices() const
{
    const VertexList *pInputVertexList{nullptr};
    try
    { // ATTN - No guarantee the list has been initialised, but silent failure here is ok
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, m_inputVertexListName, pInputVertexList) || !pInputVertexList)
            return;

        for (const Vertex *pInputVertex : *pInputVertexList)
        {
            PandoraContentApi::Vertex::Parameters parameters;
            parameters.m_position = pInputVertex->GetPosition();
            parameters.m_vertexLabel = VERTEX_INTERACTION;
            parameters.m_vertexType = VERTEX_3D;

            const Vertex *pVertex(nullptr);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertex));
        }
    }
    catch (const StatusCodeException &)
    {
        return;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::AddToSlidingFitCache(const Cluster *const pCluster)
{
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const TwoDSlidingFitResult slidingFitResult(pCluster, m_slidingFitWindow, slidingFitPitch);

    if (!m_slidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pCluster, slidingFitResult)).second)
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const TwoDSlidingFitResult &CandidateVertexCreationAlgorithm::GetCachedSlidingFitResult(const Cluster *const pCluster) const
{
    TwoDSlidingFitResultMap::const_iterator iter = m_slidingFitResultMap.find(pCluster);

    if (m_slidingFitResultMap.end() == iter)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CandidateVertexCreationAlgorithm::TidyUp()
{
    m_slidingFitResultMap.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CandidateVertexCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputClusterListNames", m_inputClusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InputVertexListName", m_inputVertexListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputVertexListName", m_outputVertexListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ReplaceCurrentVertexList", m_replaceCurrentVertexList));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterCaloHits", m_minClusterCaloHits));

    float minClusterLength = std::sqrt(m_minClusterLengthSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterLength", minClusterLength));
    m_minClusterLengthSquared = minClusterLength * minClusterLength;

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ChiSquaredCut", m_chiSquaredCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "EnableEndpointCandidates", m_enableEndpointCandidates));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxEndpointXDiscrepancy", m_maxEndpointXDiscrepancy));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "EnableCrossingCandidates", m_enableCrossingCandidates));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NMaxCrossingCandidates", m_nMaxCrossingCandidates));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxCrossingXDiscrepancy", m_maxCrossingXDiscrepancy));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ExtrapolationNSteps", m_extrapolationNSteps));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ExtrapolationStepSize", m_extrapolationStepSize));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ReducedCandidates", m_reducedCandidates));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SelectionCutFactorMax", m_selectionCutFactorMax));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "NClustersPassingMaxCutsPar", m_nClustersPassingMaxCutsPar));

    float maxCrossingSeparation = std::sqrt(m_maxCrossingSeparationSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxCrossingSeparation", maxCrossingSeparation));
    m_maxCrossingSeparationSquared = maxCrossingSeparation * maxCrossingSeparation;

    float minNearbyCrossingDistance = std::sqrt(m_minNearbyCrossingDistanceSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinNearbyCrossingDistance", minNearbyCrossingDistance));
    m_minNearbyCrossingDistanceSquared = minNearbyCrossingDistance * minNearbyCrossingDistance;

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
