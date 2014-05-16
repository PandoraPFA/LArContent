/**
 *  @file   LArContent/src/LArTwoDReco/LArCosmicRay/CosmicRaySplittingAlgorithm.cc
 * 
 *  @brief  Implementation of the cosmic ray splitting algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArTwoDReco/LArCosmicRay/CosmicRaySplittingAlgorithm.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArPointingClusterHelper.h"

using namespace pandora;

namespace lar
{

StatusCode CosmicRaySplittingAlgorithm::Run()
{
   std::cout << " --- CosmicRaySplittingAlgorithm::Run() --- " << std::endl;


    const ClusterList *pClusterList = NULL;    
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));


    ClusterVector clusterVector;
    this->GetListOfCleanClusters(pClusterList, clusterVector);
    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);


    TwoDSlidingFitResultList slidingFitResultList;
    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
    {
         TwoDSlidingFitResult slidingFitResult;
         LArClusterHelper::LArTwoDSlidingFit(*iter, m_halfWindowLayers, slidingFitResult);
         slidingFitResultList.push_back(slidingFitResult);
    }


    for (TwoDSlidingFitResultList::const_iterator iter1 = slidingFitResultList.begin(), iterEnd1 = slidingFitResultList.end(); 
        iter1 != iterEnd1; ++iter1)
    {
        const TwoDSlidingFitResult &slidingFitResult1 = *iter1;

        CartesianVector splitPosition(0.f,0.f,0.f);
        CartesianVector splitDirection1(0.f,0.f,0.f);
        CartesianVector splitDirection2(0.f,0.f,0.f);

        if (STATUS_CODE_SUCCESS != this->FindBestSplitPosition(slidingFitResult1, splitPosition, splitDirection1, splitDirection2))
            continue;

        bool foundConfirmation(false);

        for (TwoDSlidingFitResultList::const_iterator iter2 = slidingFitResultList.begin(), iterEnd2 = slidingFitResultList.end(); iter2 != iterEnd2; ++iter2)
        {
            const TwoDSlidingFitResult &slidingFitResult2 = *iter2;

            if (slidingFitResult1.GetCluster() == slidingFitResult2.GetCluster())
                continue;

            if (STATUS_CODE_SUCCESS != this->ConfirmSplitPosition(slidingFitResult2, splitPosition, splitDirection1, splitDirection2))
                continue;

// --- BEGIN EVENT DISPLAY ---
ClusterList tempClusterList1, tempClusterList2, tempClusterList3;
Cluster* pCluster1 = (Cluster*)slidingFitResult1.GetCluster();
Cluster* pCluster2 = (Cluster*)slidingFitResult2.GetCluster();
tempClusterList1.insert(pCluster1);
tempClusterList2.insert(pCluster2);
for (ClusterList::const_iterator iter3 = pClusterList->begin(), iterEnd3 = pClusterList->end(); iter3 != iterEnd3; ++iter3)
{
Cluster *pCluster3 = *iter3;
if( pCluster1 != pCluster3 && pCluster2 !=pCluster3)
tempClusterList3.insert(pCluster3);
}
PandoraMonitoringApi::SetEveDisplayParameters(false, DETECTOR_VIEW_XZ);
PandoraMonitoringApi::VisualizeClusters(&tempClusterList1, "ClusterToSplit", RED);
PandoraMonitoringApi::VisualizeClusters(&tempClusterList2, "ClusterToMerge", BLUE);
PandoraMonitoringApi::VisualizeClusters(&tempClusterList3, "OtherClusters", GREEN);
PandoraMonitoringApi::AddMarkerToVisualization(&splitPosition,"SplitPosition",BLACK,2.5);
PandoraMonitoringApi::ViewEvent();
// --- END EVENT DISPLAY --- 

            foundConfirmation = true;
            break;
        }

        if (!foundConfirmation)
            continue;

        // SPLIT...

    }


    return STATUS_CODE_SUCCESS;  
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRaySplittingAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (LArClusterHelper::GetLengthSquared(pCluster) < m_clusterMinLength * m_clusterMinLength)
            continue;

        clusterVector.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRaySplittingAlgorithm::FindBestSplitPosition(const TwoDSlidingFitResult &slidingFitResult, CartesianVector &splitPosition,
    CartesianVector &splitDirection1, CartesianVector &splitDirection2) const
{

  
    const float m_maxScatterCosTheta(0.99); // <-------- CUT!

        
    float splitCosTheta(m_maxScatterCosTheta);
    bool foundSplit(false);

    const CartesianVector &minPosition(slidingFitResult.GetGlobalMinLayerPosition());
    const CartesianVector &maxPosition(slidingFitResult.GetGlobalMaxLayerPosition());
    const float halfWindowLength(slidingFitResult.GetLayerFitHalfWindowLength());

    float minL(0.f), maxL(0.f), minT(0.f), maxT(0.f);
    slidingFitResult.GetLocalPosition(minPosition, minL, minT);
    slidingFitResult.GetLocalPosition(maxPosition, maxL, maxT);

    const unsigned int nSamplingPoints = static_cast<unsigned int>((maxL - minL)/ m_samplingPitch);

    for (unsigned int n = 0; n < nSamplingPoints; ++n)
    {
        const float alpha((0.5f + static_cast<float>(n)) / static_cast<float>(nSamplingPoints));
        const float rL(minL + (maxL - minL) * alpha);

        try
        { 
            CartesianVector centralPosition(0.f,0.f,0.f);
            CartesianVector forwardDirection(0.f,0.f,0.f);
            CartesianVector backwardDirection(0.f,0.f,0.f);

            slidingFitResult.GetGlobalFitPosition(rL, centralPosition);
            slidingFitResult.GetGlobalFitDirection(rL + halfWindowLength, forwardDirection);
            slidingFitResult.GetGlobalFitDirection(rL - halfWindowLength, backwardDirection);

            const float cosTheta(forwardDirection.GetDotProduct(backwardDirection));
            
            if (cosTheta < splitCosTheta)
            {
                CartesianVector forwardPosition(0.f,0.f,0.f);
                CartesianVector backwardPosition(0.f,0.f,0.f);

                slidingFitResult.GetGlobalFitPosition(rL + halfWindowLength, forwardPosition);
                slidingFitResult.GetGlobalFitPosition(rL - halfWindowLength, backwardPosition);
                
                CartesianVector forwardVectorInwards(centralPosition - forwardPosition);
                CartesianVector backwardVectorInwards(centralPosition - backwardPosition); 

                splitPosition = centralPosition;
                splitDirection1 = (forwardDirection.GetDotProduct(forwardVectorInwards) > 0.f) ? forwardDirection : forwardDirection * -1.f;
                splitDirection2 = (backwardDirection.GetDotProduct(backwardVectorInwards) > 0.f) ? backwardDirection : backwardDirection * -1.f;
                splitCosTheta = cosTheta;
                foundSplit = true;
            }  
        }
        catch (StatusCodeException &)
        {

        }
    }

    if (!foundSplit)
        return STATUS_CODE_NOT_FOUND;

    std::cout << " ---> Found Split: cosTheta=" << splitCosTheta << std::endl;

    return STATUS_CODE_SUCCESS;
}      
 
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRaySplittingAlgorithm::ConfirmSplitPosition(const TwoDSlidingFitResult &slidingFitResult, const CartesianVector &splitPosition,
    const CartesianVector &splitDirection1, const CartesianVector &splitDirection2) const
{
    const CartesianVector minPosition(slidingFitResult.GetGlobalMinLayerPosition());
    const CartesianVector maxPosition(slidingFitResult.GetGlobalMaxLayerPosition());
    const CartesianVector minDirection(slidingFitResult.GetGlobalMinLayerDirection());
    const CartesianVector maxDirection(slidingFitResult.GetGlobalMaxLayerDirection());

    for (unsigned int iFwd = 0; iFwd < 2; ++iFwd)
    {
        const CartesianVector pAxis((0 == iFwd) ? (maxPosition - minPosition) : (minPosition - maxPosition));
        const CartesianVector vtxPosition((0 == iFwd) ? minPosition : maxPosition);
        const CartesianVector vtxDirection((0 == iFwd) ? (pAxis.GetDotProduct(minDirection) > 0 ? minDirection : minDirection * -1.f) :
            (pAxis.GetDotProduct(maxDirection) > 0 ? maxDirection : maxDirection * -1.f));

        if((vtxPosition - splitPosition).GetMagnitudeSquared() > 0.25 * (maxPosition - minPosition).GetMagnitudeSquared())
            continue;

        if ((pAxis.GetDotProduct(splitDirection1) < 0.f || vtxDirection.GetDotProduct(splitDirection1) < 0.966) && 
            (pAxis.GetDotProduct(splitDirection2) < 0.f || vtxDirection.GetDotProduct(splitDirection2) < 0.966)) // <-------- CUT!
            continue;
                     
        float impactL(0.f), impactT(0.f);
        LArPointingClusterHelper::GetImpactParameters(vtxPosition, vtxDirection, splitPosition, impactL, impactT);
                
        if (impactL < -1.f || impactT > 5.f) // <-------- CUT!
            continue;

        std::cout << " ---> Found confirmation: rL=" << impactL << " rT=" << impactT << std::endl;

        return STATUS_CODE_SUCCESS;
    }
  
    return STATUS_CODE_NOT_FOUND;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRaySplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_clusterMinLength = 10.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterMinLength", m_clusterMinLength));

    m_halfWindowLayers = 30;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitHalfWindow", m_halfWindowLayers));

    m_samplingPitch = 1.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SamplingPitch", m_samplingPitch));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
