/**
 *  @file
 * larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewThreeDKinkBaseTool.cc
 *
 *  @brief  Implementation of the three d kink base tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArObjects/LArPointingCluster.h"

#include "larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewThreeDKinkBaseTool.h"

using namespace pandora;

namespace lar_content
{

TwoViewThreeDKinkBaseTool::TwoViewThreeDKinkBaseTool(const unsigned int nCommonClusters) :
    m_nCommonClusters(nCommonClusters),
    m_majorityRulesMode(false),
    m_minXOverlapFraction(0.1f),
    m_minMatchingScore(0.95f),
    m_minLocallyMatchedFraction(0.3f),
    m_minLongitudinalImpactParameter(-1.f),
    m_nLayersForKinkSearch(10),
    m_additionalXStepForKinkSearch(0.01f)
{
    if (!(1 == m_nCommonClusters))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

TwoViewThreeDKinkBaseTool::~TwoViewThreeDKinkBaseTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewThreeDKinkBaseTool::PassesElementCuts(MatrixType::ElementList::const_iterator eIter, const ClusterSet &usedClusters) const
{
    if (usedClusters.count(eIter->GetCluster1()) || usedClusters.count(eIter->GetCluster2()))
        return false;

    if (eIter->GetOverlapResult().GetTwoViewXOverlap().GetXOverlapFraction0() - m_minXOverlapFraction < -1.f * std::numeric_limits<float>::epsilon())
        return false;
    if (eIter->GetOverlapResult().GetTwoViewXOverlap().GetXOverlapFraction1() - m_minXOverlapFraction < -1.f * std::numeric_limits<float>::epsilon())
        return false;

    if (eIter->GetOverlapResult().GetMatchingScore() - m_minMatchingScore < std::numeric_limits<float>::epsilon())
        return false;

    if (eIter->GetOverlapResult().GetLocallyMatchedFraction() - m_minLocallyMatchedFraction < std::numeric_limits<float>::epsilon())
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float TwoViewThreeDKinkBaseTool::GetXSamplingPoint(const CartesianVector &splitPosition1, const bool isForwardInX,
    const TwoDSlidingFitResult &fitResult1, const TwoDSlidingFitResult &fitResult2) const
{
    // Nearest common x position
    float xMin1(std::numeric_limits<float>::max()), xMax1(-std::numeric_limits<float>::max());
    float xMin2(std::numeric_limits<float>::max()), xMax2(-std::numeric_limits<float>::max());
    fitResult1.GetMinAndMaxX(xMin1, xMax1);
    fitResult2.GetMinAndMaxX(xMin2, xMax2);

    const float commonX(isForwardInX ? std::max(xMin1, xMin2) : std::min(xMax1, xMax2));

    if (isForwardInX && ((commonX > xMax1) || (commonX > xMax2)))
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    if (!isForwardInX && ((commonX < xMin1) || (commonX < xMin2)))
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    // Layer step x position
    float rL1(0.f), rT1(0.f);
    fitResult1.GetLocalPosition(splitPosition1, rL1, rT1);
    const int splitLayer(fitResult1.GetLayer(rL1));

    const int lowLayer(std::max(fitResult1.GetMinLayer(), std::min(fitResult1.GetMaxLayer(), splitLayer - m_nLayersForKinkSearch)));
    const int highLayer(std::max(fitResult1.GetMinLayer(), std::min(fitResult1.GetMaxLayer(), splitLayer + m_nLayersForKinkSearch)));

    CartesianVector minus(0.f, 0.f, 0.f), plus(0.f, 0.f, 0.f);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, fitResult1.GetGlobalFitPosition(fitResult1.GetL(lowLayer), minus));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, fitResult1.GetGlobalFitPosition(fitResult1.GetL(highLayer), plus));

    if (minus.GetX() > plus.GetX())
    {
        CartesianVector temporary(minus);
        minus = plus;
        plus = temporary;
    }

    const float layerStepX(isForwardInX ? plus.GetX() : minus.GetX());

    // Final x position selection
    const float chosenX(isForwardInX ? std::max(layerStepX, commonX) : std::min(layerStepX, commonX));
    const float finalX(isForwardInX ? chosenX + m_additionalXStepForKinkSearch : chosenX - m_additionalXStepForKinkSearch);
    return finalX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewThreeDKinkBaseTool::IsALowestInX(const LArPointingCluster &pointingClusterA, const LArPointingCluster &pointingClusterB)
{
    if ((pointingClusterA.GetInnerVertex().GetPosition().GetX() < pointingClusterB.GetInnerVertex().GetPosition().GetX()) &&
        (pointingClusterA.GetInnerVertex().GetPosition().GetX() < pointingClusterB.GetOuterVertex().GetPosition().GetX()))
    {
        return true;
    }

    if ((pointingClusterA.GetOuterVertex().GetPosition().GetX() < pointingClusterB.GetInnerVertex().GetPosition().GetX()) &&
        (pointingClusterA.GetOuterVertex().GetPosition().GetX() < pointingClusterB.GetOuterVertex().GetPosition().GetX()))
    {
        return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewThreeDKinkBaseTool::Run(TwoViewTransverseTracksAlgorithm *const pAlgorithm, MatrixType &overlapMatrix)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    ModificationList modificationList;
    this->GetModifications(pAlgorithm, overlapMatrix, modificationList);
    const bool changesMade(this->ApplyChanges(pAlgorithm, modificationList));

    return changesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewThreeDKinkBaseTool::GetModifications(
    TwoViewTransverseTracksAlgorithm *const pAlgorithm, const MatrixType &overlapMatrix, ModificationList &modificationList) const
{
    ClusterSet usedClusters;
    ClusterVector sortedKeyClusters;
    overlapMatrix.GetSortedKeyClusters(sortedKeyClusters);

    for (const Cluster *const pKeyCluster : sortedKeyClusters)
    {
        if (!pKeyCluster->IsAvailable())
            continue;

        unsigned int n1(0), n2(0);
        MatrixType::ElementList elementList;
        overlapMatrix.GetConnectedElements(pKeyCluster, true, elementList, n1, n2);

        if (n1 * n2 < 2)
            continue;

        for (MatrixType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
        {
            if (!this->PassesElementCuts(eIter, usedClusters))
                continue;

            IteratorList iteratorList;
            this->SelectMatrixElements(eIter, elementList, usedClusters, iteratorList);

            if (iteratorList.size() < 2)
                continue;

            ModificationList localModificationList;
            this->GetIteratorListModifications(pAlgorithm, iteratorList, localModificationList);

            if (localModificationList.empty())
                continue;

            for (ModificationList::const_iterator mIter = localModificationList.begin(), mIterEnd = localModificationList.end(); mIter != mIterEnd; ++mIter)
            {
                usedClusters.insert(mIter->m_affectedClusters.begin(), mIter->m_affectedClusters.end());
            }

            modificationList.insert(modificationList.end(), localModificationList.begin(), localModificationList.end());
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewThreeDKinkBaseTool::ApplyChanges(TwoViewTransverseTracksAlgorithm *const pAlgorithm, const ModificationList &modificationList) const
{
    ClusterMergeMap consolidatedMergeMap;
    SplitPositionMap consolidatedSplitMap;

    for (const Modification &modification : modificationList)
    {
        ClusterList parentClusters;
        for (const auto &mapEntry : modification.m_clusterMergeMap)
            parentClusters.push_back(mapEntry.first);
        parentClusters.sort(LArClusterHelper::SortByNHits);

        for (const Cluster *const pParentCluster : parentClusters)
        {
            const ClusterList &daughterClusters(modification.m_clusterMergeMap.at(pParentCluster));

            for (const Cluster *const pDaughterCluster : daughterClusters)
            {
                if (consolidatedMergeMap.count(pDaughterCluster))
                    throw StatusCodeException(STATUS_CODE_FAILURE);
            }

            ClusterList &targetClusterList(consolidatedMergeMap[pParentCluster]);
            targetClusterList.insert(targetClusterList.end(), daughterClusters.begin(), daughterClusters.end());
        }

        ClusterList splitClusters;
        for (const auto &mapEntry : modification.m_splitPositionMap)
            splitClusters.push_back(mapEntry.first);
        splitClusters.sort(LArClusterHelper::SortByNHits);

        for (const Cluster *const pSplitCluster : splitClusters)
        {
            const CartesianPointVector &splitPositions(modification.m_splitPositionMap.at(pSplitCluster));

            CartesianPointVector &cartesianPointVector(consolidatedSplitMap[pSplitCluster]);
            cartesianPointVector.insert(cartesianPointVector.end(), splitPositions.begin(), splitPositions.end());
        }
    }

    bool changesMade(false);
    changesMade |= pAlgorithm->MakeClusterMerges(consolidatedMergeMap);
    changesMade |= pAlgorithm->MakeClusterSplits(consolidatedSplitMap);

    return changesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewThreeDKinkBaseTool::SelectMatrixElements(MatrixType::ElementList::const_iterator eIter,
    const MatrixType::ElementList &elementList, const ClusterSet &usedClusters, IteratorList &iteratorList) const
{
    iteratorList.push_back(eIter);

    for (MatrixType::ElementList::const_iterator eIter2 = elementList.begin(); eIter2 != elementList.end(); ++eIter2)
    {
        if (eIter == eIter2)
            continue;

        if (!this->PassesElementCuts(eIter2, usedClusters))
            continue;

        for (IteratorList::const_iterator iIter = iteratorList.begin(); iIter != iteratorList.end(); ++iIter)
        {
            if ((*iIter) == eIter2)
                continue;

            unsigned int nMatchedClusters(0);

            if ((*iIter)->GetCluster1() == eIter2->GetCluster1())
                ++nMatchedClusters;

            if ((*iIter)->GetCluster2() == eIter2->GetCluster2())
                ++nMatchedClusters;

            if (m_nCommonClusters == nMatchedClusters)
            {
                iteratorList.push_back(eIter2);
                return;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoViewThreeDKinkBaseTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MajorityRulesMode", m_majorityRulesMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinMatchingScore",
            m_minMatchingScore)); // from calo matching - clear tracks tool

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinLocallyMatchedFraction",
            m_minLocallyMatchedFraction)); // from calo matching
                                           // - clear tracks tool

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinXOverlapFraction",
            m_minXOverlapFraction)); // from calo matching - clear tracks tool

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinLongitudinalImpactParameter", m_minLongitudinalImpactParameter));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NLayersForKinkSearch", m_nLayersForKinkSearch));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "AdditionalXStepForKinkSearch", m_additionalXStepForKinkSearch));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
