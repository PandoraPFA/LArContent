/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/VertexSplittingAlgorithm.cc
 *
 *  @brief  Implementation of the vertex splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/VertexSplittingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

VertexSplittingAlgorithm::VertexSplittingAlgorithm() :
    m_splitDisplacementSquared(4.f * 4.f),
    m_vertexDisplacementSquared(1.f * 1.f)
{
    // ATTN Some default values differ from base class
    m_minClusterLength = 1.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexSplittingAlgorithm::FindBestSplitPosition(const TwoDSlidingFitResult &slidingFitResult, CartesianVector &splitPosition) const
{
    // Identify event vertex
    const VertexList *pVertexList(NULL);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pVertexList));

    if (pVertexList->empty())
        return STATUS_CODE_NOT_INITIALIZED;

    if (pVertexList->size() != 1)
        return STATUS_CODE_OUT_OF_RANGE;

    const Cluster *const pCluster(slidingFitResult.GetCluster());
    const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));

    const Vertex *const pSelectedVertex(*(pVertexList->begin()));

    if (VERTEX_3D != pSelectedVertex->GetVertexType())
        return STATUS_CODE_INVALID_PARAMETER;

    const CartesianVector theVertex2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pSelectedVertex->GetPosition(), hitType));

    const CartesianVector innerVertex2D(slidingFitResult.GetGlobalMinLayerPosition());
    const CartesianVector outerVertex2D(slidingFitResult.GetGlobalMaxLayerPosition());

    if ((outerVertex2D - innerVertex2D).GetMagnitudeSquared() < 4.f * m_vertexDisplacementSquared)
        return STATUS_CODE_NOT_FOUND;

    bool foundSplit(false);
    const StatusCode statusCode(slidingFitResult.GetGlobalFitProjection(theVertex2D, splitPosition));

    if (STATUS_CODE_SUCCESS != statusCode)
        return statusCode;

    const float splitDisplacementSquared((splitPosition - theVertex2D).GetMagnitudeSquared());
    const float vertexDisplacementSquared(
        std::min((splitPosition - innerVertex2D).GetMagnitudeSquared(), (splitPosition - outerVertex2D).GetMagnitudeSquared()));

    if ((splitDisplacementSquared < m_splitDisplacementSquared) && (vertexDisplacementSquared > m_vertexDisplacementSquared) &&
        (splitDisplacementSquared < vertexDisplacementSquared))
    {
        foundSplit = true;
    }

    if (!foundSplit)
        return STATUS_CODE_NOT_FOUND;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    float splitDisplacement = std::sqrt(m_splitDisplacementSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SplitDisplacement", splitDisplacement));
    m_splitDisplacementSquared = splitDisplacement * splitDisplacement;

    float vertexDisplacement = std::sqrt(m_vertexDisplacementSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexDisplacement", vertexDisplacement));
    m_vertexDisplacementSquared = vertexDisplacement * vertexDisplacement;

    return TwoDSlidingFitSplittingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
