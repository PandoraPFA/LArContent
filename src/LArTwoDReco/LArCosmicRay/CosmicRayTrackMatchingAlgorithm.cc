/**
 *  @file   LArContent/src/LArTwoDReco/LArCosmicRay/CosmicRayTrackMatchingAlgorithm.cc
 *
 *  @brief  Implementation of the cosmic ray splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArTwoDReco/LArCosmicRay/CosmicRayTrackMatchingAlgorithm.h"

#include "LArHelpers/LArGeometryHelper.h"
#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArThreeDHelper.h"

using namespace pandora;

namespace lar
{

StatusCode CosmicRayTrackMatchingAlgorithm::Run()
{



    ClusterVector availableClustersU, availableClustersV, availableClustersW;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetAvailableClusters(m_inputClusterListNamesU, availableClustersU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetAvailableClusters(m_inputClusterListNamesV, availableClustersV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetAvailableClusters(m_inputClusterListNamesW, availableClustersW));

    ClusterVector cleanClustersU, cleanClustersV, cleanClustersW;
    this->SelectCleanClusters(availableClustersU, cleanClustersU);
    this->SelectCleanClusters(availableClustersV, cleanClustersV);
    this->SelectCleanClusters(availableClustersW, cleanClustersW);

    TwoDSlidingFitResultMap slidingFitResultMap;
    this->AddToSlidingFitResultMap(cleanClustersU, slidingFitResultMap);
    this->AddToSlidingFitResultMap(cleanClustersV, slidingFitResultMap);
    this->AddToSlidingFitResultMap(cleanClustersW, slidingFitResultMap);



    this->SelectMatchedTracks(slidingFitResultMap, cleanClustersU, cleanClustersV, availableClustersW);
    this->SelectMatchedTracks(slidingFitResultMap, cleanClustersV, cleanClustersW, availableClustersU);
    this->SelectMatchedTracks(slidingFitResultMap, cleanClustersW, cleanClustersU, availableClustersV);




    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayTrackMatchingAlgorithm::GetAvailableClusters(const StringVector inputClusterListNames, ClusterVector &clusterVector) const
{
    for (StringVector::const_iterator sIter = inputClusterListNames.begin(), sIterEnd = inputClusterListNames.end(); sIter != sIterEnd; ++sIter)
    {
        const ClusterList *pClusterList = NULL;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, *sIter, pClusterList))

        for (ClusterList::const_iterator cIter = pClusterList->begin(), cIterEnd = pClusterList->end(); cIter != cIterEnd; ++cIter)
        {
            Cluster *pCluster = *cIter;
            if (!pCluster->IsAvailable())
                continue;

            clusterVector.push_back(pCluster);
        }
    }

    if (clusterVector.empty())
        return STATUS_CODE_NOT_FOUND;

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTrackMatchingAlgorithm::SelectCleanClusters(const ClusterVector &inputVector, ClusterVector &outputVector) const
{
    for (ClusterVector::const_iterator iter = inputVector.begin(), iterEnd = inputVector.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (LArClusterHelper::GetLengthSquared(pCluster) < m_clusterMinLength * m_clusterMinLength)
            continue;

        outputVector.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTrackMatchingAlgorithm::AddToSlidingFitResultMap(const ClusterVector &clusterVector,
    TwoDSlidingFitResultMap &slidingFitResultMap) const
{
    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
    {
        if (slidingFitResultMap.end() == slidingFitResultMap.find(*iter))
        {
            TwoDSlidingFitResult slidingFitResult;
            LArClusterHelper::LArTwoDSlidingFit(*iter, m_halfWindowLayers, slidingFitResult);

            if (!slidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(*iter, slidingFitResult)).second)
                throw StatusCodeException(STATUS_CODE_FAILURE);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTrackMatchingAlgorithm::SelectMatchedTracks(const TwoDSlidingFitResultMap &slidingFitResultMap,
    const ClusterVector &clusterVector1, const ClusterVector &clusterVector2, const ClusterVector &clusterVector3)
{
    if (clusterVector1.empty() || clusterVector2.empty() || clusterVector3.empty())
        return;

    const HitType hitType1(LArThreeDHelper::GetClusterHitType(*clusterVector1.begin()));
    const HitType hitType2(LArThreeDHelper::GetClusterHitType(*clusterVector2.begin()));
    const HitType hitType3(LArThreeDHelper::GetClusterHitType(*clusterVector3.begin()));

    if (hitType1 == hitType2 || hitType2 == hitType3 || hitType3 == hitType1)
        return;

    for (ClusterVector::const_iterator iter1 = clusterVector1.begin(), iterEnd1 = clusterVector1.end(); iter1 != iterEnd1; ++iter1)
    {
        Cluster* pCluster1 = *iter1;

        TwoDSlidingFitResultMap::const_iterator sIter1 = slidingFitResultMap.find(pCluster1);
        if (slidingFitResultMap.end() == sIter1)
            continue;

        const TwoDSlidingFitResult &slidingFitResult1(sIter1->second);

        for (ClusterVector::const_iterator iter2 = clusterVector2.begin(), iterEnd2 = clusterVector2.end(); iter2 != iterEnd2; ++iter2)
        {
            Cluster* pCluster2 = *iter2;

            TwoDSlidingFitResultMap::const_iterator sIter2 = slidingFitResultMap.find(pCluster2);
            if (slidingFitResultMap.end() == sIter2)
                continue;

            const TwoDSlidingFitResult &slidingFitResult2(sIter2->second);



            this->SelectMatchedTracks(slidingFitResult1, slidingFitResult2, clusterVector3);




        }
    }

}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTrackMatchingAlgorithm::SelectMatchedTracks(const TwoDSlidingFitResult &slidingFitResult1,
    const TwoDSlidingFitResult &slidingFitResult2, const ClusterVector &availableClusters3)
{


    const Cluster* pCluster1(slidingFitResult1.GetCluster());
    const Cluster* pCluster2(slidingFitResult2.GetCluster());

    const HitType hitType1(LArThreeDHelper::GetClusterHitType(pCluster1));
    const HitType hitType2(LArThreeDHelper::GetClusterHitType(pCluster2));

    float xMin1(0.f), xMax1(0.f), xMin2(0.f), xMax2(0.f);
    LArClusterHelper::GetClusterSpanX(pCluster1, xMin1, xMax1);
    LArClusterHelper::GetClusterSpanX(pCluster2, xMin2, xMax2);

    const float xOverlap(std::min(xMax1,xMax2) - std::max(xMin1,xMin2));
    const float xSpan(std::max(xMax1,xMax2) - std::min(xMin1,xMin2));

    if (xOverlap < m_minXOverlap || xOverlap/xSpan < m_minXOverlapFraction)
        return;

std::cout << " *** CosmicRayTrackMatchingAlgorithm *** " << std::endl;
std::cout << "   xOverlap=" << xOverlap << " xOverlapFraction=" << xOverlap/xSpan << std::endl;

// --- BEGIN EVENT DISPLAY ---
 ClusterList tempList1, tempList2, tempList3;
tempList1.insert((Cluster*)pCluster1);
tempList2.insert((Cluster*)pCluster2);
tempList3.insert(availableClusters3.begin(), availableClusters3.end());
PandoraMonitoringApi::SetEveDisplayParameters(false, DETECTOR_VIEW_XZ);
PandoraMonitoringApi::VisualizeClusters(&tempList1, "Cluster1", RED);
PandoraMonitoringApi::VisualizeClusters(&tempList2, "Cluster2", BLUE);
PandoraMonitoringApi::VisualizeClusters(&tempList3, "Cluster3", GREEN);


    // Sampling in x
    const float xMin(std::max(xMin1,xMin2));
    const float xMax(std::min(xMax1,xMax2));

    unsigned int nSamplingPoints(100);

    for (unsigned int n = 0; n < nSamplingPoints; ++n)
    {
        const float alpha((0.5f + static_cast<float>(n)) / static_cast<float>(nSamplingPoints));
        const float x(xMin + alpha * (xMax - xMin));

        try
        {
            float chi2(0.f);
            CartesianVector position1(0.f, 0.f, 0.f), position2(0.f, 0.f, 0.f), position3(0.f, 0.f, 0.f);
            slidingFitResult1.GetGlobalFitPosition(x, true, position1);
            slidingFitResult2.GetGlobalFitPosition(x, true, position2);
            LArGeometryHelper::MergeTwoPositions(hitType1, hitType2, position1, position2, position3, chi2);

            PandoraMonitoringApi::AddMarkerToVisualization(&position3, "Projection", GREEN, 2.f);
        }
        catch (StatusCodeException &)
        {
        }
    }

PandoraMonitoringApi::ViewEvent();
// --- END EVENT DISPLAY ---

}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayTrackMatchingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "InputClusterListNamesU", m_inputClusterListNamesU));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "InputClusterListNamesV", m_inputClusterListNamesV));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "InputClusterListNamesW", m_inputClusterListNamesW));

    m_clusterMinLength = 10.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterMinLength", m_clusterMinLength));

    m_halfWindowLayers = 15;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitHalfWindow", m_halfWindowLayers));

    m_minXOverlap = 3.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinXOverlap", m_minXOverlap));

    m_minXOverlapFraction = 0.8f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinXOverlapFraction", m_minXOverlapFraction));


    return STATUS_CODE_SUCCESS;
}

} // namespace lar
