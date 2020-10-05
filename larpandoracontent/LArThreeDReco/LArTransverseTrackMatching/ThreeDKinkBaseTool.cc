/**
 *  @file   larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/ThreeDKinkBaseTool.cc
 *
 *  @brief  Implementation of the three d kink base tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArObjects/LArPointingCluster.h"

#include "larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/ThreeDKinkBaseTool.h"

using namespace pandora;

namespace lar_content
{

ThreeDKinkBaseTool::ThreeDKinkBaseTool(const unsigned int nCommonClusters) :
    m_nCommonClusters(nCommonClusters),
    m_majorityRulesMode(false),
    m_minMatchedFraction(0.75f),
    m_minMatchedSamplingPoints(10),
    m_minLongitudinalImpactParameter(-1.f),
    m_nLayersForKinkSearch(10),
    m_additionalXStepForKinkSearch(0.01f)
{
    if (!((1 == m_nCommonClusters) || (2 == m_nCommonClusters)))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDKinkBaseTool::~ThreeDKinkBaseTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ThreeDKinkBaseTool::PassesElementCuts(TensorType::ElementList::const_iterator eIter, const ClusterSet &usedClusters) const
{
    if (usedClusters.count(eIter->GetClusterU()) || usedClusters.count(eIter->GetClusterV()) || usedClusters.count(eIter->GetClusterW()))
        return false;

    if (eIter->GetOverlapResult().GetMatchedFraction() < m_minMatchedFraction)
        return false;

    if (eIter->GetOverlapResult().GetNMatchedSamplingPoints() < m_minMatchedSamplingPoints)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ThreeDKinkBaseTool::GetXSamplingPoint(const CartesianVector &splitPosition1, const bool isForwardInX, const TwoDSlidingFitResult &fitResult1,
    const TwoDSlidingFitResult &fitResult2, const TwoDSlidingFitResult &fitResult3) const
{
    // Nearest common x position
    float xMin1(std::numeric_limits<float>::max()), xMax1(-std::numeric_limits<float>::max());
    float xMin2(std::numeric_limits<float>::max()), xMax2(-std::numeric_limits<float>::max());
    float xMin3(std::numeric_limits<float>::max()), xMax3(-std::numeric_limits<float>::max());
    fitResult1.GetMinAndMaxX(xMin1, xMax1);
    fitResult2.GetMinAndMaxX(xMin2, xMax2);
    fitResult3.GetMinAndMaxX(xMin3, xMax3);

    const float commonX(isForwardInX ? std::max(xMin1, std::max(xMin2, xMin3)) : std::min(xMax1, std::min(xMax2, xMax3)));

    if (isForwardInX && ((commonX > xMax1) || (commonX > xMax2) || (commonX > xMax3)))
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    if (!isForwardInX && ((commonX < xMin1) || (commonX < xMin2) || (commonX < xMin3)))
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

bool ThreeDKinkBaseTool::IsALowestInX(const LArPointingCluster &pointingClusterA, const LArPointingCluster &pointingClusterB)
{
    if ((pointingClusterA.GetInnerVertex().GetPosition().GetX() < pointingClusterB.GetInnerVertex().GetPosition().GetX()) &&
        (pointingClusterA.GetInnerVertex().GetPosition().GetX() < pointingClusterB.GetOuterVertex().GetPosition().GetX()) )
    {
        return true;
    }

    if ((pointingClusterA.GetOuterVertex().GetPosition().GetX() < pointingClusterB.GetInnerVertex().GetPosition().GetX()) &&
        (pointingClusterA.GetOuterVertex().GetPosition().GetX() < pointingClusterB.GetOuterVertex().GetPosition().GetX()) )
    {
        return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ThreeDKinkBaseTool::Run(ThreeViewTransverseTracksAlgorithm *const pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    ModificationList modificationList;
    this->GetModifications(pAlgorithm, overlapTensor, modificationList);
    const bool changesMade(this->ApplyChanges(pAlgorithm, modificationList));

    return changesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDKinkBaseTool::GetModifications(ThreeViewTransverseTracksAlgorithm *const pAlgorithm, const TensorType &overlapTensor, ModificationList &modificationList) const
{
    ClusterSet usedClusters;
    ClusterVector sortedKeyClusters;
    overlapTensor.GetSortedKeyClusters(sortedKeyClusters);

    for (const Cluster *const pKeyCluster : sortedKeyClusters)
    {
        if (!pKeyCluster->IsAvailable())
            continue;

        unsigned int nU(0), nV(0), nW(0);
        TensorType::ElementList elementList;
        overlapTensor.GetConnectedElements(pKeyCluster, true, elementList, nU, nV, nW);

        if (nU * nV * nW < 2)
            continue;

        for (TensorType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
        {
            if (!this->PassesElementCuts(eIter, usedClusters))
                continue;

            IteratorList iteratorList;
            this->SelectTensorElements(eIter, elementList, usedClusters, iteratorList);

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

bool ThreeDKinkBaseTool::ApplyChanges(ThreeViewTransverseTracksAlgorithm *const pAlgorithm, const ModificationList &modificationList) const
{
    ClusterMergeMap consolidatedMergeMap;
    SplitPositionMap consolidatedSplitMap;

    for (const Modification &modification : modificationList)
    {
        ClusterList parentClusters;
        for (const auto &mapEntry : modification.m_clusterMergeMap) parentClusters.push_back(mapEntry.first);
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
        for (const auto &mapEntry : modification.m_splitPositionMap) splitClusters.push_back(mapEntry.first);
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

void ThreeDKinkBaseTool::SelectTensorElements(TensorType::ElementList::const_iterator eIter, const TensorType::ElementList &elementList,
    const ClusterSet &usedClusters, IteratorList &iteratorList) const
{
    iteratorList.push_back(eIter);

    for (TensorType::ElementList::const_iterator eIter2 = elementList.begin(); eIter2 != elementList.end(); ++eIter2)
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

            if ((*iIter)->GetClusterU() == eIter2->GetClusterU())
                ++nMatchedClusters;

            if ((*iIter)->GetClusterV() == eIter2->GetClusterV())
                ++nMatchedClusters;

            if ((*iIter)->GetClusterW() == eIter2->GetClusterW())
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

StatusCode ThreeDKinkBaseTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MajorityRulesMode", m_majorityRulesMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedFraction", m_minMatchedFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedSamplingPoints", m_minMatchedSamplingPoints));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinLongitudinalImpactParameter", m_minLongitudinalImpactParameter));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NLayersForKinkSearch", m_nLayersForKinkSearch));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AdditionalXStepForKinkSearch", m_additionalXStepForKinkSearch));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
