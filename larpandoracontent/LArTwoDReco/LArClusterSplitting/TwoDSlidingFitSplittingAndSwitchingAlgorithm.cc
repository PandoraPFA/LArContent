/**
 *  @file   larpandoracontent/LArTwoDReco/ClusterSplitting/TwoDSlidingFitSplittingAndSwitchingAlgorithm.cc
 *
 *  @brief  Implementation of the two dimensional sliding fit splitting and splicing algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/TwoDSlidingFitSplittingAndSwitchingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

TwoDSlidingFitSplittingAndSwitchingAlgorithm::TwoDSlidingFitSplittingAndSwitchingAlgorithm() :
    m_halfWindowLayers(25),
    m_minClusterLength(10.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

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

    // May choose to cache information here, for subsequent expensive calculations
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->PreparationStep(clusterVector));

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

            CartesianVector splitPosition(0.f, 0.f, 0.f);
            CartesianVector firstDirection(0.f, 0.f, 0.f);
            CartesianVector secondDirection(0.f, 0.f, 0.f);

            if (STATUS_CODE_SUCCESS != this->FindBestSplitPosition(slidingFitResult1, slidingFitResult2, splitPosition, firstDirection, secondDirection))
                continue;

            const Cluster *const pCluster1 = slidingFitResult1.GetCluster();
            const Cluster *const pCluster2 = slidingFitResult2.GetCluster();

            if (STATUS_CODE_SUCCESS != this->ReplaceClusters(pCluster1, pCluster2, splitPosition, firstDirection, secondDirection))
                continue;

            slidingFitResultMap.erase(sIter1);
            slidingFitResultMap.erase(sIter2);

            *iter1 = NULL;
            *iter2 = NULL;

            break;
        }
    }

    return this->TidyUpStep();
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitSplittingAndSwitchingAlgorithm::PreparationStep(const ClusterVector & /*clusterVector*/)
{
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitSplittingAndSwitchingAlgorithm::TidyUpStep()
{
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitSplittingAndSwitchingAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pCluster = *iter;

        if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLength * m_minClusterLength)
            continue;

        clusterVector.push_back(pCluster);
    }

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitSplittingAndSwitchingAlgorithm::BuildSlidingFitResultMap(
    const ClusterVector &clusterVector, TwoDSlidingFitResultMap &slidingFitResultMap) const
{
    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
    {
        if (slidingFitResultMap.end() == slidingFitResultMap.find(*iter))
        {
            try
            {
                const float slidingFitPitch(LArGeometryHelper::GetWirePitch(this->GetPandora(), LArClusterHelper::GetClusterHitType(*iter)));
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

void TwoDSlidingFitSplittingAndSwitchingAlgorithm::SplitCluster(const Cluster *const pCluster, const CartesianVector &splitPosition,
    const CartesianVector &splitDirection, CaloHitList &firstCaloHitList, CaloHitList &secondCaloHitList) const
{
    CaloHitList caloHitsToDistribute;
    pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitsToDistribute);

    for (CaloHitList::const_iterator iter = caloHitsToDistribute.begin(), iterEnd = caloHitsToDistribute.end(); iter != iterEnd; ++iter)
    {
        const CaloHit *const pCaloHit = *iter;

        if (splitDirection.GetDotProduct((pCaloHit->GetPositionVector() - splitPosition)) > 0.f)
        {
            firstCaloHitList.push_back(pCaloHit);
        }
        else
        {
            secondCaloHitList.push_back(pCaloHit);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitSplittingAndSwitchingAlgorithm::ReplaceClusters(const Cluster *const pCluster1, const Cluster *const pCluster2,
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
    clusterList.push_back(pCluster1);
    clusterList.push_back(pCluster2);

    std::string clusterListToSaveName, clusterListToDeleteName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        PandoraContentApi::InitializeFragmentation(*this, clusterList, clusterListToDeleteName, clusterListToSaveName));

    // Create new clusters
    const Cluster *pFirstCluster(NULL), *pSecondCluster(NULL);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, firstParameters, pFirstCluster));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, secondParameters, pSecondCluster));

    // End cluster fragmentation operations
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, clusterListToSaveName, clusterListToDeleteName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitSplittingAndSwitchingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "HalfWindowLayers", m_halfWindowLayers));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterLength", m_minClusterLength));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
