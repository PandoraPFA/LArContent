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



   
    

    for (TwoDSlidingFitResultList::const_iterator iter = slidingFitResultList.begin(), iterEnd = slidingFitResultList.end(); iter != iterEnd; ++iter)
    {
        const TwoDSlidingFitResult &slidingFitResult = *iter;
        const TwoDSlidingFitResult::LayerFitResultMap &layerFitResultMap(slidingFitResult.GetLayerFitResultMap());

        const CartesianVector &minPosition(slidingFitResult.GetGlobalMinLayerPosition());
        const CartesianVector &maxPosition(slidingFitResult.GetGlobalMaxLayerPosition());
        const CartesianVector axisDirection((maxPosition - minPosition).GetUnitVector());

        float maxDisplacement(0.f);
        CartesianVector maxDisplacementPosition(0.f,0.f,0.f);

        for (TwoDSlidingFitResult::LayerFitResultMap::const_iterator iter = layerFitResultMap.begin(), iterEnd = layerFitResultMap.end();
            iter != iterEnd; ++iter)
        {
            const float rL(iter->second.GetL());
            const float rT(iter->second.GetFitT());

            CartesianVector thisPosition(0.f, 0.f, 0.f);
            slidingFitResult.GetGlobalPosition(rL, rT, thisPosition);

            const float thisDisplacement(axisDirection.GetCrossProduct(thisPosition - minPosition).GetMagnitude());

    
            if (thisDisplacement > maxDisplacement)
            {
                maxDisplacement = thisDisplacement;
                maxDisplacementPosition = thisPosition;
            }

        }

// --- BEGIN EVENT DISPLAY ---
ClusterList tempClusterList;
Cluster* pCluster = (Cluster*)slidingFitResult.GetCluster();
tempClusterList.insert(pCluster);
PandoraMonitoringApi::VisualizeClusters(&tempClusterList, "Clusters", RED);
PandoraMonitoringApi::AddMarkerToVisualization(&maxDisplacementPosition,"MaxDisplacement",RED,2.5);
PandoraMonitoringApi::ViewEvent();
// --- END EVENT DISPLAY ---


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

StatusCode CosmicRaySplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_clusterMinLength = 5.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterMinLength", m_clusterMinLength));

    m_halfWindowLayers = 25;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitHalfWindow", m_halfWindowLayers));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
