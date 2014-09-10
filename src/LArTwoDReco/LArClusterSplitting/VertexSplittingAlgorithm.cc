/**
 *  @file   LArContent/src/LArTwoDReco/LArClusterSplitting/VertexSplittingAlgorithm.cc
 * 
 *  @brief  Implementation of the vertex splitting algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArTwoDReco/LArClusterSplitting/VertexSplittingAlgorithm.h"

using namespace pandora;

namespace lar_content
{
 
StatusCode VertexSplittingAlgorithm::FindBestSplitPosition(const TwoDSlidingFitResult &slidingFitResult, CartesianVector &splitPosition) const
{
    // Identify event vertex
    const VertexList *pVertexList(NULL);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pVertexList));

    Vertex *pSelectedVertex(NULL);
    for (VertexList::const_iterator iter = pVertexList->begin(), iterEnd = pVertexList->end(); iter != iterEnd; ++iter)
    {
        // TODO Vertex selection
        return STATUS_CODE_NOT_INITIALIZED;
    }

    if (!pSelectedVertex)
        return STATUS_CODE_SUCCESS;

    const CartesianVector &theVertex(pSelectedVertex->GetPosition());
    const CartesianVector innerVertex(slidingFitResult.GetGlobalMinLayerPosition());
    const CartesianVector outerVertex(slidingFitResult.GetGlobalMaxLayerPosition());

    if ((outerVertex - innerVertex).GetMagnitudeSquared() < 4.f * m_vertexDisplacementSquared)
        return STATUS_CODE_NOT_FOUND;

    bool foundSplit(false);

    try
    {
        slidingFitResult.GetGlobalFitProjection(theVertex, splitPosition); 

        const float splitDisplacementSquared((splitPosition - theVertex).GetMagnitudeSquared());
        const float vertexDisplacementSquared(std::min((splitPosition - innerVertex).GetMagnitudeSquared(), (splitPosition - outerVertex).GetMagnitudeSquared()));

        if (splitDisplacementSquared < m_splitDisplacementSquared &&
            vertexDisplacementSquared > m_vertexDisplacementSquared &&
            splitDisplacementSquared < vertexDisplacementSquared )
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
    m_splitDisplacement = 4.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SplitDisplacement", m_splitDisplacement));
    m_splitDisplacementSquared = m_splitDisplacement * m_splitDisplacement;

    m_vertexDisplacement = 2.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VertexDisplacement", m_vertexDisplacement));
    m_vertexDisplacementSquared = m_vertexDisplacement * m_vertexDisplacement;

    return TwoDSlidingFitSplittingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
