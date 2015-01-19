/**
 *  @file   LArContent/src/LArTwoDReco/LArClusterSplitting/TwoDSlidingFitSplittingAlgorithm.cc
 *
 *  @brief  Implementation of the two dimensional sliding fit splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"

#include "LArTwoDReco/LArClusterSplitting/TwoDSlidingFitSplittingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

TwoDSlidingFitSplittingAlgorithm::TwoDSlidingFitSplittingAlgorithm() :
    m_slidingFitHalfWindow(20),
    m_minClusterLength(10.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitSplittingAlgorithm::DivideCaloHits(const Cluster *const pCluster, CaloHitList &firstHitList, CaloHitList &secondHitList) const
{
    if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLength * m_minClusterLength)
        return STATUS_CODE_NOT_FOUND;

    try
    {
        const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));

        const TwoDSlidingFitResult slidingFitResult(pCluster, m_slidingFitHalfWindow, slidingFitPitch);
        CartesianVector splitPosition(0.f, 0.f, 0.f);

        if (STATUS_CODE_SUCCESS == this->FindBestSplitPosition(slidingFitResult, splitPosition))
        {
            return this->DivideCaloHits(slidingFitResult, splitPosition, firstHitList, secondHitList);
        }
    }
    catch (StatusCodeException &statusCodeException)
    {
        if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
            throw statusCodeException;
    }

    return STATUS_CODE_NOT_FOUND;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitSplittingAlgorithm::DivideCaloHits(const TwoDSlidingFitResult &slidingFitResult, const CartesianVector &splitPosition, 
    CaloHitList &firstCaloHitList, CaloHitList &secondCaloHitList) const
{
    float rL(0.f), rT(0.f);
    slidingFitResult.GetLocalPosition(splitPosition, rL, rT);

    const Cluster *pCluster(slidingFitResult.GetCluster());
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
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitHalfWindow", m_slidingFitHalfWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLength", m_minClusterLength));

    return ClusterSplittingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
