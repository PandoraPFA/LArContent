/**
 *  @file   LArContent/src/LArTwoDReco/LArClusterSplitting/TwoDSlidingFitSplittingAlgorithm.cc
 *
 *  @brief  Implementation of the two dimensional sliding fit splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArTwoDReco/LArClusterSplitting/TwoDSlidingFitSplittingAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode TwoDSlidingFitSplittingAlgorithm::SplitCluster(const Cluster *const pCluster, CaloHitList &firstHitList, CaloHitList &secondHitList) const
{
    if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLength * m_minClusterLength)
        return STATUS_CODE_NOT_FOUND;

    TwoDSlidingFitResult slidingFitResult;
    LArClusterHelper::LArTwoDSlidingFit(pCluster, m_slidingFitHalfWindow, slidingFitResult);

    CartesianVector splitPosition(0.f, 0.f, 0.f);

    if (STATUS_CODE_SUCCESS == this->FindBestSplitPosition(slidingFitResult, splitPosition))
        return this->SplitCluster(slidingFitResult, splitPosition, firstHitList, secondHitList);

    return STATUS_CODE_NOT_FOUND;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitSplittingAlgorithm::SplitCluster(const TwoDSlidingFitResult &slidingFitResult, const CartesianVector &splitPosition, 
    CaloHitList &firstCaloHitList, CaloHitList &secondCaloHitList) const
{
    float rL(0.f), rT(0.f);
    slidingFitResult.GetLocalPosition(splitPosition, rL, rT);

    const Cluster* pCluster = slidingFitResult.GetCluster();
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(); iter != orderedCaloHitList.end(); ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            CaloHit *pCaloHit = *hitIter;

            float thisL(0.f), thisT(0.f);
            slidingFitResult.GetLocalPosition(pCaloHit->GetPositionVector(), thisL, thisT);

            if (thisL < rL)
            {
                firstCaloHitList.insert(pCaloHit);
            }
            else
            {
                secondCaloHitList.insert(pCaloHit);
            }
        }
    }

    if (firstCaloHitList.empty() || secondCaloHitList.empty())
        return STATUS_CODE_NOT_FOUND;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_slidingFitHalfWindow = 20;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitHalfWindow", m_slidingFitHalfWindow));

    m_minClusterLength = 10.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLength", m_minClusterLength));

    return ClusterSplittingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
