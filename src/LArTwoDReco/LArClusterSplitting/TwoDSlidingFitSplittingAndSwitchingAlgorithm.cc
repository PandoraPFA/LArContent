/**
 *  @file   LArContent/src/LArTwoDReco/ClusterSplitting/TwoDSlidingFitSplittingAndSwitchingAlgorithm.cc
 *
 *  @brief  Implementation of the two dimensional sliding fit splitting and splicing algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"

#include "LArPlugins/LArTransformationPlugin.h"

#include "LArTwoDReco/LArClusterSplitting/TwoDSlidingFitSplittingAndSwitchingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode TwoDSlidingFitSplittingAndSwitchingAlgorithm::Run()
{
    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

    // Get ordered list of clean clusters
    ClusterVector clusterVector;
    this->GetListOfCleanClusters(pClusterList, clusterVector);

    // Calculate sliding fit results for clean clusters
    TwoDSlidingFitResultMap slidingFitResultMap;
    this->BuildSlidingFitResultMap(clusterVector, slidingFitResultMap);

    // Loop over clusters, identify split positions, perform splits
    for (ClusterVector::iterator iter1 = clusterVector.begin(), iterEnd1 = clusterVector.end(); iter1 != iterEnd1; ++iter1)
    {
        if (NULL == *iter1)
            continue;

        TwoDSlidingFitResultMap::iterator sIter1 = slidingFitResultMap.find(*iter1);

        if (slidingFitResultMap.end() == sIter1)
            continue;

        const TwoDSlidingFitResult &slidingFitResult1(sIter1->second);

        for (ClusterVector::iterator iter2 = iter1, iterEnd2 = iterEnd1; iter2 != iterEnd2; ++iter2)
        {
            if (NULL == *iter2)
                continue;

            TwoDSlidingFitResultMap::iterator sIter2 = slidingFitResultMap.find(*iter2);

            if (slidingFitResultMap.end() == sIter2)
                continue;

            const TwoDSlidingFitResult &slidingFitResult2(sIter2->second);

            if (slidingFitResult1.GetCluster() == slidingFitResult2.GetCluster())
                continue;

            CartesianVector splitPosition(0.f,0.f,0.f);
            CartesianVector firstDirection(0.f,0.f,0.f);
            CartesianVector secondDirection(0.f,0.f,0.f);

            if (STATUS_CODE_SUCCESS != this->FindBestSplitPosition(slidingFitResult1, slidingFitResult2, splitPosition, firstDirection, secondDirection))
                continue;

            Cluster *pCluster1 = const_cast<Cluster*>(slidingFitResult1.GetCluster());
            Cluster *pCluster2 = const_cast<Cluster*>(slidingFitResult2.GetCluster());

            if (STATUS_CODE_SUCCESS != this->ReplaceClusters(pCluster1, pCluster2, splitPosition, firstDirection, secondDirection))
                continue;

            slidingFitResultMap.erase(sIter1);
            slidingFitResultMap.erase(sIter2);

            *iter1 = NULL;
            *iter2 = NULL;

            break;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitSplittingAndSwitchingAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLength * m_minClusterLength)
            continue;

        clusterVector.push_back(pCluster);
    }

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitSplittingAndSwitchingAlgorithm::BuildSlidingFitResultMap(const ClusterVector &clusterVector,
    TwoDSlidingFitResultMap &slidingFitResultMap) const
{
    const float slidingFitPitch(LArGeometryHelper::GetLArTransformationPlugin(this->GetPandora())->GetWireZPitch());

    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
    {
        if (slidingFitResultMap.end() == slidingFitResultMap.find(*iter))
        {
            try
            {
                const TwoDSlidingFitResult slidingFitResult(*iter, m_halfWindowLayers, slidingFitPitch);

                if (!slidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(*iter, slidingFitResult)).second)
                    throw StatusCodeException(STATUS_CODE_FAILURE);
            }
            catch (StatusCodeException &statusCodeException)
            {
                if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                    throw statusCodeException;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitSplittingAndSwitchingAlgorithm::SplitCluster(Cluster *const pCluster, const CartesianVector &splitPosition,
    const CartesianVector &splitDirection, CaloHitList &firstCaloHitList, CaloHitList &secondCaloHitList) const
{
    CaloHitList caloHitsToDistribute;
    pCluster->GetOrderedCaloHitList().GetCaloHitList(caloHitsToDistribute);

    for (CaloHitList::const_iterator iter = caloHitsToDistribute.begin(), iterEnd = caloHitsToDistribute.end(); iter != iterEnd; ++iter)
    {
        CaloHit *pCaloHit = *iter;

        if (splitDirection.GetDotProduct((pCaloHit->GetPositionVector() - splitPosition)) > 0.f)
        {
            firstCaloHitList.insert(pCaloHit);
        }
        else
        {
            secondCaloHitList.insert(pCaloHit);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitSplittingAndSwitchingAlgorithm::ReplaceClusters(Cluster *const pCluster1, Cluster *const pCluster2,
    const CartesianVector &splitPosition, const CartesianVector &firstDirection, const CartesianVector &secondDirection) const
{
    // Split cluster into two hit lists (note the convention for 'firstDirection' and 'secondDirection')
    PandoraContentApi::Cluster::Parameters firstParameters, secondParameters;

    this->SplitCluster(pCluster1, splitPosition, firstDirection, firstParameters.m_caloHitList, secondParameters.m_caloHitList);
    this->SplitCluster(pCluster2, splitPosition, secondDirection, secondParameters.m_caloHitList, firstParameters.m_caloHitList);

    if (firstParameters.m_caloHitList.empty() || secondParameters.m_caloHitList.empty())
        return STATUS_CODE_NOT_ALLOWED;

    // Begin cluster fragmentation operations
    ClusterList clusterList;
    clusterList.insert(pCluster1);
    clusterList.insert(pCluster2);

    std::string clusterListToSaveName, clusterListToDeleteName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*this, clusterList,
        clusterListToDeleteName, clusterListToSaveName));

    // Create new clusters
    Cluster *pFirstCluster(NULL), *pSecondCluster(NULL);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, firstParameters, pFirstCluster));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, secondParameters, pSecondCluster));

    // End cluster fragmentation operations
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, clusterListToSaveName, clusterListToDeleteName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitSplittingAndSwitchingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_halfWindowLayers = 25;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "HalfWindowLayers", m_halfWindowLayers));

    m_minClusterLength = 10.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLength", m_minClusterLength));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
