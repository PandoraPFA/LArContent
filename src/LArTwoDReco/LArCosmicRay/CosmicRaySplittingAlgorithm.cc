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



    CartesianVector splitPosition(0.f,0.f,0.f);

    for (TwoDSlidingFitResultList::const_iterator iter1 = slidingFitResultList.begin(), iterEnd1 = slidingFitResultList.end(); iter1 != iterEnd1; ++iter1)
    {
        const TwoDSlidingFitResult &slidingFitResult1 = *iter1;

        for (TwoDSlidingFitResultList::const_iterator iter2 = iter1, iterEnd2 = iterEnd1; iter2 != iterEnd2; ++iter2)
        {
            const TwoDSlidingFitResult &slidingFitResult2 = *iter2;

            if (slidingFitResult1.GetCluster() == slidingFitResult2.GetCluster())
                continue;

            if (STATUS_CODE_SUCCESS == this->SplitCrossedClusters(slidingFitResult1, slidingFitResult2, splitPosition))
            {
// --- BEGIN EVENT DISPLAY ---
ClusterList tempList1, tempList2;
Cluster* pCluster1 = (Cluster*)slidingFitResult1.GetCluster();
Cluster* pCluster2 = (Cluster*)slidingFitResult2.GetCluster();
tempList1.insert(pCluster1);
tempList2.insert(pCluster2);
PandoraMonitoringApi::VisualizeClusters(&tempList1, "Cluster1", RED);
PandoraMonitoringApi::VisualizeClusters(&tempList2, "Cluster2", BLUE);
PandoraMonitoringApi::AddMarkerToVisualization(&splitPosition,"SplitPosition",RED,2.5);
PandoraMonitoringApi::ViewEvent();
// --- END EVENT DISPLAY ---
            }

        }
    }
    



    /*
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
    */






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

StatusCode CosmicRaySplittingAlgorithm::SplitCrossedClusters(const TwoDSlidingFitResult &slidingFitResult1, 
    const TwoDSlidingFitResult &slidingFitResult2, CartesianVector &splitPosition) const
{
    //
    // TODO: SPEED UP!
    //

    const float m_maxSeparationSquared = 2.5f * 2.5f;

    bool foundSplit(false);
    float closestSeparationSquared(m_maxSeparationSquared);

    const float halfWindowLength1(slidingFitResult1.GetLayerFitHalfWindowLength());
    const float halfWindowLength2(slidingFitResult2.GetLayerFitHalfWindowLength());

    const TwoDSlidingFitResult::LayerFitResultMap &layerFitResultMap1(slidingFitResult1.GetLayerFitResultMap());
    const TwoDSlidingFitResult::LayerFitResultMap &layerFitResultMap2(slidingFitResult2.GetLayerFitResultMap());

    for (TwoDSlidingFitResult::LayerFitResultMap::const_iterator iter1 = layerFitResultMap1.begin(), iterEnd1 = layerFitResultMap1.end();
        iter1 != iterEnd1; ++iter1)
    {
        const float rL1(iter1->second.GetL());
        const float rT1(iter1->second.GetFitT());

        CartesianVector R1(0.f, 0.f, 0.f);
        CartesianVector F1(0.f, 0.f, 0.f);
        CartesianVector B1(0.f, 0.f, 0.f);

        try
        {
            slidingFitResult1.GetGlobalPosition(rL1, rT1, R1);
            slidingFitResult1.GetGlobalFitPosition(rL1 + halfWindowLength1, F1);
            slidingFitResult1.GetGlobalFitPosition(rL1 - halfWindowLength1, B1);
        }
        catch (StatusCodeException &)
        {
            continue;
        }
        
        for (TwoDSlidingFitResult::LayerFitResultMap::const_iterator iter2 = layerFitResultMap2.begin(), iterEnd2 = layerFitResultMap2.end();
            iter2 != iterEnd2; ++iter2)
        {
            const float rL2(iter2->second.GetL());
            const float rT2(iter2->second.GetFitT());

            CartesianVector R2(0.f, 0.f, 0.f);
            CartesianVector F2(0.f, 0.f, 0.f);
            CartesianVector B2(0.f, 0.f, 0.f);

            try
            {
                slidingFitResult2.GetGlobalPosition(rL2, rT2, R2);
                slidingFitResult2.GetGlobalFitPosition(rL2 + halfWindowLength2, F2);
                slidingFitResult2.GetGlobalFitPosition(rL2 - halfWindowLength2, B2);
            }
            catch (StatusCodeException &)
            {
                continue;
            }

            if ((R1 - R2).GetMagnitudeSquared() > m_maxSeparationSquared)
                continue;

            const CartesianVector C0((R1 + R2) * 0.5);

            const CartesianVector a1(B1);
            const CartesianVector a2(F1); 

            for (unsigned int iForward = 0; iForward<2; ++iForward)
            {
                CartesianVector b1((0 == iForward) ? F2 : B2);
                CartesianVector b2((0 == iForward) ? B2 : F2);

                if ((b1 - C0).GetDotProduct(a1 - C0) > 0 || (b2 - C0).GetDotProduct(a2 - C0) > 0)
                    continue;

                //
                // ADD MORE SELECTION HERE
                //

                try{
                    float mu1(0.f), mu2(0.f);
                    CartesianVector C1(0.f,0.f,0.f);

                    const CartesianVector p1((b1 - a1).GetUnitVector());
                    const CartesianVector p2((b2 - a2).GetUnitVector());
                    LArPointingClusterHelper::GetIntersection(a1, p1, a2, p2, C1, mu1, mu2);
            
                    const float thisSeparationSquared((C0 - C1).GetMagnitudeSquared());

                    if (thisSeparationSquared < closestSeparationSquared)
                    {
                        closestSeparationSquared = thisSeparationSquared;
                        splitPosition = (C0 + C1) * 0.5;
                        foundSplit = true;
                    }
                }
                catch (StatusCodeException &)
                {
                    
                }

            }
        }
    }

    if( foundSplit )
    {
        ClusterList tempList1, tempList2;
        Cluster* pCluster1 = (Cluster*)slidingFitResult1.GetCluster();
        Cluster* pCluster2 = (Cluster*)slidingFitResult2.GetCluster();
        tempList1.insert(pCluster1);
        tempList2.insert(pCluster2);
        PandoraMonitoringApi::VisualizeClusters(&tempList1, "Cluster1", RED);
        PandoraMonitoringApi::VisualizeClusters(&tempList2, "Cluster2", BLUE);
        PandoraMonitoringApi::AddMarkerToVisualization(&splitPosition,"SplitPosition",BLACK,2.5);
        PandoraMonitoringApi::ViewEvent();
    }

    return STATUS_CODE_NOT_FOUND;  
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
