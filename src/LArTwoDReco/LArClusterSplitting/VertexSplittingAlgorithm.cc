/**
 *  @file   LArContent/src/LArTwoDReco/LArClusterSplitting/VertexSplittingAlgorithm.cc
 * 
 *  @brief  Implementation of the vertex splitting algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"

#include "LArTwoDReco/LArClusterSplitting/VertexSplittingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

VertexSplittingAlgorithm::VertexSplittingAlgorithm() :
    m_splitDisplacementSquared(4.f * 4.f),
    m_vertexDisplacementSquared(2.f * 2.f)
{
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

    const Cluster *pCluster(slidingFitResult.GetCluster());
    const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));

    const Vertex *pSelectedVertex(*(pVertexList->begin()));

    if (VERTEX_3D != pSelectedVertex->GetVertexType())
        return STATUS_CODE_INVALID_PARAMETER;

    const CartesianVector theVertex2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pSelectedVertex->GetPosition(), hitType));

    const CartesianVector innerVertex2D(slidingFitResult.GetGlobalMinLayerPosition());
    const CartesianVector outerVertex2D(slidingFitResult.GetGlobalMaxLayerPosition());

    if ((outerVertex2D - innerVertex2D).GetMagnitudeSquared() < 4.f * m_vertexDisplacementSquared)
        return STATUS_CODE_NOT_FOUND;

    bool foundSplit(false);

    try
    {
        slidingFitResult.GetGlobalFitProjection(theVertex2D, splitPosition); 

        const float splitDisplacementSquared((splitPosition - theVertex2D).GetMagnitudeSquared());
        const float vertexDisplacementSquared(std::min((splitPosition - innerVertex2D).GetMagnitudeSquared(), (splitPosition - outerVertex2D).GetMagnitudeSquared()));

        if ((splitDisplacementSquared < m_splitDisplacementSquared) && (vertexDisplacementSquared > m_vertexDisplacementSquared) &&
            (splitDisplacementSquared < vertexDisplacementSquared))
        {
            foundSplit = true;
        }
    }
    catch (StatusCodeException &)
    {
    }

    if (!foundSplit)
        return STATUS_CODE_NOT_FOUND;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    float splitDisplacement = std::sqrt(m_splitDisplacementSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SplitDisplacement", splitDisplacement));
    m_splitDisplacementSquared = splitDisplacement * splitDisplacement;

    float vertexDisplacement = std::sqrt(m_vertexDisplacementSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VertexDisplacement", vertexDisplacement));
    m_vertexDisplacementSquared = vertexDisplacement * vertexDisplacement;

    return TwoDSlidingFitSplittingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
