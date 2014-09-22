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
    {
        // TODO Vertex selection
        std::cout << "VertexSplittingAlgorithm: vertex selection not yet implemented " << std::endl;
        return STATUS_CODE_OUT_OF_RANGE;
    }

    const Vertex *pSelectedVertex(*(pVertexList->begin()));
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
