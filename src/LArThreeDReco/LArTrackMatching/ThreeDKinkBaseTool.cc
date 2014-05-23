/**
 *  @file   LArContent/src/LArThreeDReco/LArTrackMatching/ThreeDKinkBaseTool.cc
 * 
 *  @brief  Implementation of the three d kink base tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArThreeDHelper.h"

#include "LArObjects/LArPointingCluster.h"

#include "LArThreeDReco/LArTrackMatching/ThreeDKinkBaseTool.h"

using namespace pandora;

namespace lar
{

ThreeDKinkBaseTool::ThreeDKinkBaseTool(const unsigned int nCommonClusters) :
    m_nCommonClusters(nCommonClusters)
{
    if (!((1 == m_nCommonClusters) || (2 == m_nCommonClusters)))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDKinkBaseTool::~ThreeDKinkBaseTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ThreeDKinkBaseTool::PassesElementCuts(TensorType::ElementList::const_iterator eIter, const ClusterList &usedClusters) const
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
    fitResult1.GetGlobalFitPosition(fitResult1.GetL(lowLayer), minus);
    fitResult1.GetGlobalFitPosition(fitResult1.GetL(highLayer), plus);

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

bool ThreeDKinkBaseTool::Run(ThreeDTransverseTracksAlgorithm *pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraSettings::ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << m_algorithmToolType << std::endl;

    ModificationList modificationList;
    this->GetModifications(pAlgorithm, overlapTensor, modificationList);
    const bool changesMade(this->ApplyChanges(pAlgorithm, modificationList));

    return changesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDKinkBaseTool::GetModifications(ThreeDTransverseTracksAlgorithm *pAlgorithm, const TensorType &overlapTensor, ModificationList &modificationList) const
{
    ClusterList usedClusters;

    for (TensorType::const_iterator iterU = overlapTensor.begin(), iterUEnd = overlapTensor.end(); iterU != iterUEnd; ++iterU)
    {
        if (!iterU->first->IsAvailable())
            continue;

        unsigned int nU(0), nV(0), nW(0);
        TensorType::ElementList elementList;
        overlapTensor.GetConnectedElements(iterU->first, true, elementList, nU, nV, nW);

        if (nU * nV * nW < 2)
            continue;

        std::sort(elementList.begin(), elementList.end(), ThreeDTransverseTracksAlgorithm::SortByNMatchedSamplingPoints);

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

bool ThreeDKinkBaseTool::ApplyChanges(ThreeDTransverseTracksAlgorithm *pAlgorithm, const ModificationList &modificationList) const
{
    ClusterMergeMap consolidatedMergeMap;
    SplitPositionMap consolidatedSplitMap;

    for (ModificationList::const_iterator iter = modificationList.begin(), iterEnd = modificationList.end(); iter != iterEnd; ++iter)
    {
        for (ClusterMergeMap::const_iterator cIter = iter->m_clusterMergeMap.begin(), cIterEnd = iter->m_clusterMergeMap.end(); cIter != cIterEnd; ++cIter)
        {
            const ClusterList &daughterClusters(cIter->second);

            for (ClusterList::const_iterator dIter = daughterClusters.begin(), dIterEnd = daughterClusters.end(); dIter != dIterEnd; ++dIter)
            {
                if (consolidatedMergeMap.count(*dIter))
                    throw StatusCodeException(STATUS_CODE_FAILURE);
            }

            ClusterList &targetClusterList(consolidatedMergeMap[cIter->first]);
            targetClusterList.insert(daughterClusters.begin(), daughterClusters.end());
        }

        for (SplitPositionMap::const_iterator cIter = iter->m_splitPositionMap.begin(), cIterEnd = iter->m_splitPositionMap.end(); cIter != cIterEnd; ++cIter)
        {
            CartesianPointList &cartesianPointList(consolidatedSplitMap[cIter->first]);
            cartesianPointList.insert(cartesianPointList.end(), cIter->second.begin(), cIter->second.end());
        }
    }

    bool changesMade(false);
    changesMade |= pAlgorithm->MakeClusterMerges(consolidatedMergeMap);
    changesMade |= pAlgorithm->MakeClusterSplits(consolidatedSplitMap);

    return changesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDKinkBaseTool::SelectTensorElements(TensorType::ElementList::const_iterator eIter, const TensorType::ElementList &elementList,
    const ClusterList &usedClusters, IteratorList &iteratorList) const
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
    m_majorityRulesMode = false;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MajorityRulesMode", m_majorityRulesMode));

    m_minMatchedFraction = 0.75f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedFraction", m_minMatchedFraction));

    m_minMatchedSamplingPoints = 10;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedSamplingPoints", m_minMatchedSamplingPoints));

    m_minLongitudinalImpactParameter = -1.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinLongitudinalImpactParameter", m_minLongitudinalImpactParameter));

    m_nLayersForKinkSearch = 10;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NLayersForKinkSearch", m_nLayersForKinkSearch));

    m_additionalXStepForKinkSearch = 0.01f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "AdditionalXStepForKinkSearch", m_additionalXStepForKinkSearch));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
