/**
 *  @file   LArContent/src/LArVertex/VertexSplittingAlgorithm.cc
 * 
 *  @brief  Implementation of the vertex splitting algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"
#include "LArHelpers/LArVertexHelper.h"

#include "LArClusterSplitting/VertexSplittingAlgorithm.h"

using namespace pandora;

namespace lar
{

bool VertexSplittingAlgorithm::IsPossibleSplit(const Cluster *const pCluster) const
{
    if (LArClusterHelper::GetLengthSquared(pCluster) < 4.f * m_vertexDisplacementSquared)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexSplittingAlgorithm::FindBestSplitPosition(const Cluster *const pCluster, CartesianVector &splitPosition) const
{
    bool foundSplit(false);

    try
    {
        const CartesianVector &theVertex(LArVertexHelper::GetCurrentVertex());

        // Project vertex onto sliding window fit
        LArClusterHelper::TwoDSlidingFitResult twoDSlidingFitResult;
        LArClusterHelper::LArTwoDSlidingFit(pCluster, m_slidingFitLayerHalfWindow, twoDSlidingFitResult);  

        const CartesianVector innerVertex(twoDSlidingFitResult.GetGlobalMinLayerPosition());
        const CartesianVector outerVertex(twoDSlidingFitResult.GetGlobalMaxLayerPosition());

        twoDSlidingFitResult.GetGlobalFitProjection(theVertex, splitPosition); 

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
    m_slidingFitLayerHalfWindow = 10;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitLayerHalfWindow", m_slidingFitLayerHalfWindow));

    m_splitDisplacement = 4.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SplitDisplacement", m_splitDisplacement));
    m_splitDisplacementSquared = m_splitDisplacement * m_splitDisplacement;

    m_vertexDisplacement = 2.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VertexDisplacement", m_vertexDisplacement));
    m_vertexDisplacementSquared = m_vertexDisplacement * m_vertexDisplacement;

    return ClusterSplittingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
