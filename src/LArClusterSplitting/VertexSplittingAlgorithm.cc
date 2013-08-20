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

#include "LArPseudoLayerCalculator.h"

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
/*
// KEEP OLD METHOD FOR NOW
StatusCode VertexSplittingAlgorithm::FindBestSplitLayer(const Cluster *const pCluster, unsigned int &splitLayer) const
{

    if (!LArVertexHelper::DoesCurrentVertexExist())
        return STATUS_CODE_NOT_FOUND;

    // Find nearest hit to vertex
    const CartesianVector &theVertex(LArVertexHelper::GetCurrentVertex());
    const CartesianVector innerCentroid(pCluster->GetCentroid(pCluster->GetInnerPseudoLayer()));
    const CartesianVector outerCentroid(pCluster->GetCentroid(pCluster->GetOuterPseudoLayer()));

    bool foundSplit(false);
    float minDisplacementSquared(m_splitDisplacementSquared);
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iterI = orderedCaloHitList.begin(), iterEndI = orderedCaloHitList.end(); iterI != iterEndI; ++iterI)
    {
        const unsigned int thisLayer(iterI->first);
        const CaloHitList *pCaloHitList(iterI->second);

        for (CaloHitList::const_iterator iterJ = pCaloHitList->begin(), iterEndJ = pCaloHitList->end(); iterJ != iterEndJ; ++iterJ)
        {
            const CaloHit *pCaloHit = *iterJ;
            const float thisDisplacementSquared((pCaloHit->GetPositionVector() - theVertex).GetMagnitudeSquared());

            if (thisDisplacementSquared > minDisplacementSquared)
                continue;

            foundSplit = false;
            minDisplacementSquared = thisDisplacementSquared;

            if (((pCaloHit->GetPositionVector() - innerCentroid).GetMagnitudeSquared() < m_vertexDisplacementSquared) ||
                ((pCaloHit->GetPositionVector() - outerCentroid).GetMagnitudeSquared() < m_vertexDisplacementSquared) )
            {
                continue;
            }

            foundSplit = true;
            splitLayer = thisLayer;
        }
    }

    if (!foundSplit)
        return STATUS_CODE_NOT_FOUND;

    return STATUS_CODE_SUCCESS;
}
*/

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexSplittingAlgorithm::FindBestSplitLayer(const Cluster *const pCluster, unsigned int &splitLayer) const
{
    if (!LArVertexHelper::DoesCurrentVertexExist())
        return STATUS_CODE_NOT_FOUND;

    const CartesianVector &theVertex(LArVertexHelper::GetCurrentVertex());

    // Project vertex onto sliding window fit
    LArClusterHelper::TwoDSlidingFitResult twoDSlidingFitResult;
    LArClusterHelper::LArTwoDSlidingFit(pCluster, 10, twoDSlidingFitResult);  

    bool foundSplit(false);
    CartesianVector splitPosition(0.f,0.f,0.f);

    try{        
        const CartesianVector innerVertex(twoDSlidingFitResult.GetGlobalMinLayerPosition());
        const CartesianVector outerVertex(twoDSlidingFitResult.GetGlobalMaxLayerPosition());

        twoDSlidingFitResult.GetGlobalFitProjection(theVertex, splitPosition); 

        const float splitDisplacementSquared((splitPosition - theVertex).GetMagnitudeSquared());
        const float vertexDisplacementSquared(std::min((splitPosition - innerVertex).GetMagnitudeSquared(),
                                                       (splitPosition - outerVertex).GetMagnitudeSquared()));

        if ( splitDisplacementSquared < m_splitDisplacementSquared &&
             vertexDisplacementSquared > m_vertexDisplacementSquared &&
             splitDisplacementSquared < vertexDisplacementSquared )
	{
	    splitLayer = GeometryHelper::GetPseudoLayer(splitPosition);
	    foundSplit = true;
	}
    }
    catch (StatusCodeException &)
    {
    }

    if (!foundSplit)
        return STATUS_CODE_NOT_FOUND;

ClusterList tempList;
Cluster* tempCluster = (Cluster*)pCluster;
tempList.insert(tempCluster);
PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
PandoraMonitoringApi::VisualizeClusters(&tempList, "Cluster", GREEN);
PandoraMonitoringApi::AddMarkerToVisualization(&theVertex, "Vertex", RED, 1.75); 
PandoraMonitoringApi::AddMarkerToVisualization(&splitPosition, "Split", BLUE, 1.75);
PandoraMonitoringApi::ViewEvent();

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

    return ClusterSplittingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
