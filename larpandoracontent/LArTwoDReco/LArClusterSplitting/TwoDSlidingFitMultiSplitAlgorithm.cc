/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/TwoDSlidingFitMultiSplitAlgorithm.cc
 *
 *  @brief  Implementation of the 2D sliding fit multi-splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/TwoDSlidingFitMultiSplitAlgorithm.h"

using namespace pandora;

namespace lar_content
{

TwoDSlidingFitMultiSplitAlgorithm::TwoDSlidingFitMultiSplitAlgorithm() :
    m_slidingFitHalfWindow(15),
    m_inputClusterList("")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitMultiSplitAlgorithm::Run()
{
    std::string originalListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*this, originalListName));

    if (!m_inputClusterList.empty())
    {
        const StatusCode statusCode(PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_inputClusterList));

        if (STATUS_CODE_NOT_FOUND == statusCode)
        {
            std::cout << "TwoDSlidingFitMultiSplitAlgorithm: cluster list not found " << m_inputClusterList << std::endl;
            return STATUS_CODE_SUCCESS;
        }

        if (STATUS_CODE_SUCCESS != statusCode)
            return statusCode;
    }

    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

    // Get ordered list of candidate clusters
    ClusterVector clusterVector;
    this->GetListOfCleanClusters(pClusterList, clusterVector);

    // Build a set of sliding fit results for clean clusters
    TwoDSlidingFitResultMap slidingFitResultMap;
    this->BuildSlidingFitResultMap(clusterVector, m_slidingFitHalfWindow, slidingFitResultMap);

    // Find best split positions for each cluster
    ClusterPositionMap clusterSplittingMap;
    this->FindBestSplitPositions(slidingFitResultMap, clusterSplittingMap);

    // Perform splits
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->SplitClusters(slidingFitResultMap, clusterSplittingMap));

    if (!m_inputClusterList.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, originalListName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitMultiSplitAlgorithm::BuildSlidingFitResultMap(
    const ClusterVector &clusterVector, const unsigned int halfWindowLayers, TwoDSlidingFitResultMap &slidingFitResultMap) const
{
    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
    {
        if (slidingFitResultMap.end() == slidingFitResultMap.find(*iter))
        {
            try
            {
                const float slidingFitPitch(LArGeometryHelper::GetWirePitch(this->GetPandora(), LArClusterHelper::GetClusterHitType(*iter)));
                const TwoDSlidingFitResult slidingFitResult(*iter, halfWindowLayers, slidingFitPitch);

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

StatusCode TwoDSlidingFitMultiSplitAlgorithm::SplitClusters(
    const TwoDSlidingFitResultMap &slidingFitResultMap, const ClusterPositionMap &clusterSplittingMap) const
{
    ClusterList clusterList;
    for (const auto &mapEntry : clusterSplittingMap)
        clusterList.push_back(mapEntry.first);
    clusterList.sort(LArClusterHelper::SortByNHits);

    for (const Cluster *const pCluster : clusterList)
    {
        const CartesianPointVector &splitPositionVector(clusterSplittingMap.at(pCluster));

        if (splitPositionVector.empty())
            continue;

        TwoDSlidingFitResultMap::const_iterator sIter = slidingFitResultMap.find(pCluster);
        if (slidingFitResultMap.end() == sIter)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        const TwoDSlidingFitResult &slidingFitResult = sIter->second;

        StatusCode statusCode(this->SplitCluster(slidingFitResult, splitPositionVector));

        if (STATUS_CODE_SUCCESS != statusCode)
            return statusCode;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitMultiSplitAlgorithm::SplitCluster(const TwoDSlidingFitResult &slidingFitResult, const CartesianPointVector &splitPositionVector) const
{
    const Cluster *const pCluster = slidingFitResult.GetCluster();

    // Get split positions for this cluster
    FloatVector displacementVector;

    for (CartesianPointVector::const_iterator pIter = splitPositionVector.begin(), pIterEnd = splitPositionVector.end(); pIter != pIterEnd; ++pIter)
    {
        const CartesianVector &splitPosition = *pIter;

        float rL(0.f), rT(0.f);
        slidingFitResult.GetLocalPosition(splitPosition, rL, rT);
        displacementVector.push_back(rL);
    }

    const float bigL(2.f * slidingFitResult.GetL(slidingFitResult.GetMaxLayer()));
    displacementVector.push_back(-bigL);
    displacementVector.push_back(+bigL);

    std::sort(displacementVector.begin(), displacementVector.end());

    // Begin cluster fragmentation operations
    const ClusterList clusterList(1, pCluster);
    std::string clusterListToSave, clusterListToDelete;

    PANDORA_RETURN_RESULT_IF(
        STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*this, clusterList, clusterListToDelete, clusterListToSave));

    CaloHitList oldCaloHitList;
    pCluster->GetOrderedCaloHitList().FillCaloHitList(oldCaloHitList);

    bool foundPreviousL(false);
    float prevL(0.f);

    for (FloatVector::const_iterator fIter = displacementVector.begin(), fIterEnd = displacementVector.end(); fIter != fIterEnd; ++fIter)
    {
        const float nextL(*fIter);

        if (foundPreviousL)
        {
            // Select hits for new cluster
            CaloHitList newCaloHitList;

            for (CaloHitList::const_iterator hIter = oldCaloHitList.begin(), hIterEnd = oldCaloHitList.end(); hIter != hIterEnd; ++hIter)
            {
                const CaloHit *const pCaloHit = *hIter;

                float rL(0.f), rT(0.f);
                slidingFitResult.GetLocalPosition(pCaloHit->GetPositionVector(), rL, rT);

                if (rL >= prevL && rL < nextL)
                    newCaloHitList.push_back(pCaloHit);
            }

            if (newCaloHitList.empty())
                continue;

            // Create new cluster
            PandoraContentApi::Cluster::Parameters newParameters;
            newParameters.m_caloHitList = newCaloHitList;

            const Cluster *pNewCluster(NULL);
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, newParameters, pNewCluster));
        }

        prevL = nextL;
        foundPreviousL = true;
    }

    // End cluster fragmentation operations
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, clusterListToSave, clusterListToDelete));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitMultiSplitAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListName", m_inputClusterList));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingFitHalfWindow", m_slidingFitHalfWindow));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
