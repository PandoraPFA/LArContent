/**
 *  @file   LArContent/src/LArVertex/VertexSplittingAlgorithm.cc
 * 
 *  @brief  Implementation of the vertex splitting algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArClusterSplitting/VertexSplittingAlgorithm.h"

#include "LArHelpers/LArVertexHelper.h"

using namespace pandora;

namespace lar
{
 
StatusCode VertexSplittingAlgorithm::FindBestSplitPosition(const LArClusterHelper::TwoDSlidingFitResult &slidingFitResult,
    CartesianVector &splitPosition) const
{
    
    const CartesianVector &theVertex(LArVertexHelper::GetCurrentVertex());

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

} // namespace lar
