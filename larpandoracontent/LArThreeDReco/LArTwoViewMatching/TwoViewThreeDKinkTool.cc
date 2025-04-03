/**
 *  @file   larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewThreeDKinkTool.cc
 *
 *  @brief  Implementation of the two view three d kink tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"
#include "larpandoracontent/LArObjects/LArPointingCluster.h"

#include "larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewThreeDKinkTool.h"

using namespace pandora;

namespace lar_content
{

TwoViewThreeDKinkTool::TwoViewThreeDKinkTool() :
    m_minXOverlapFraction(0.1f),
    m_minMatchingScore(0.95f),
    m_minLocallyMatchedFraction(0.3f),
    m_minLongitudinalImpactParameter(-1.f),
    m_nLayersForKinkSearch(10),
    m_additionalXStepForKinkSearch(0.01f),
    m_maxTransverseImpactParameter(5.f),
    m_minImpactParameterCosTheta(0.5f),
    m_cosThetaCutForKinkSearch(0.75f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TwoViewThreeDKinkTool::~TwoViewThreeDKinkTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewThreeDKinkTool::Run(TwoViewTransverseTracksAlgorithm *const pAlgorithm, MatrixType &overlapMatrix)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    ModificationList modificationList;
    this->GetModifications(pAlgorithm, overlapMatrix, modificationList);
    const bool changesMade(this->ApplyChanges(pAlgorithm, modificationList));

    return changesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewThreeDKinkTool::PassesElementCuts(MatrixType::ElementList::const_iterator eIter, const ClusterSet &usedClusters) const
{
    if (usedClusters.count(eIter->GetCluster1()) || usedClusters.count(eIter->GetCluster2()))
        return false;

    if (eIter->GetOverlapResult().GetTwoViewXOverlap().GetXOverlapFraction0() - m_minXOverlapFraction < std::numeric_limits<float>::epsilon())
        return false;

    if (eIter->GetOverlapResult().GetTwoViewXOverlap().GetXOverlapFraction1() - m_minXOverlapFraction < std::numeric_limits<float>::epsilon())
        return false;

    if (eIter->GetOverlapResult().GetMatchingScore() - m_minMatchingScore < std::numeric_limits<float>::epsilon())
        return false;

    if (eIter->GetOverlapResult().GetLocallyMatchedFraction() - m_minLocallyMatchedFraction < std::numeric_limits<float>::epsilon())
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float TwoViewThreeDKinkTool::GetXSamplingPoint(const CartesianVector &splitPosition1, const bool isForwardInX,
    const TwoDSlidingFitResult &fitResult1, const TwoDSlidingFitResult &fitResult2) const
{
    // Nearest common x position
    float xMin1(std::numeric_limits<float>::max()), xMax1(std::numeric_limits<float>::lowest());
    float xMin2(std::numeric_limits<float>::max()), xMax2(std::numeric_limits<float>::lowest());
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
        std::swap(minus, plus);

    const float layerStepX(isForwardInX ? plus.GetX() : minus.GetX());

    // Final x position selection
    const float chosenX(isForwardInX ? std::max(layerStepX, commonX) : std::min(layerStepX, commonX));
    const float finalX(isForwardInX ? chosenX + m_additionalXStepForKinkSearch : chosenX - m_additionalXStepForKinkSearch);
    return finalX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewThreeDKinkTool::IsALowestInX(const LArPointingCluster &pointingClusterA, const LArPointingCluster &pointingClusterB)
{
    const float xMinA{std::min(pointingClusterA.GetInnerVertex().GetPosition().GetX(), pointingClusterA.GetOuterVertex().GetPosition().GetX())};
    const float xMinB{std::min(pointingClusterB.GetInnerVertex().GetPosition().GetX(), pointingClusterB.GetOuterVertex().GetPosition().GetX())};
    return xMinA < xMinB;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewThreeDKinkTool::GetModifications(
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

bool TwoViewThreeDKinkTool::ApplyChanges(TwoViewTransverseTracksAlgorithm *const pAlgorithm, const ModificationList &modificationList) const
{
    ClusterMergeMap consolidatedMergeMap;
    SplitPositionMap consolidatedSplitMap;

    for (const Modification &modification : modificationList)
    {
        ClusterList parentClusters;
        for (const auto &mapEntry : modification.m_clusterMergeMap)
            parentClusters.emplace_back(mapEntry.first);
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
            splitClusters.emplace_back(mapEntry.first);
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

void TwoViewThreeDKinkTool::SelectMatrixElements(MatrixType::ElementList::const_iterator eIter, const MatrixType::ElementList &elementList,
    const ClusterSet &usedClusters, IteratorList &iteratorList) const
{
    iteratorList.emplace_back(eIter);

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

            if (nMatchedClusters > 1)
                continue;

            if (nMatchedClusters)
            {
                iteratorList.emplace_back(eIter2);
                return;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewThreeDKinkTool::GetIteratorListModifications(
    TwoViewTransverseTracksAlgorithm *const pAlgorithm, const IteratorList &iteratorList, ModificationList &modificationList) const
{
    for (IteratorList::const_iterator iIter1 = iteratorList.begin(), iIter1End = iteratorList.end(); iIter1 != iIter1End; ++iIter1)
    {
        for (IteratorList::const_iterator iIter2 = std::next(iIter1); iIter2 != iIter1End; ++iIter2)
        {
            try
            {
                const unsigned int nMatchedReUpsampledSamplingPoints1((*iIter1)->GetOverlapResult().GetNMatchedReUpsampledSamplingPoints());
                const unsigned int nMatchedReUpsampledSamplingPoints2((*iIter2)->GetOverlapResult().GetNMatchedReUpsampledSamplingPoints());
                IteratorList::const_iterator iIterA((nMatchedReUpsampledSamplingPoints1 >= nMatchedReUpsampledSamplingPoints2) ? iIter1 : iIter2);
                IteratorList::const_iterator iIterB((nMatchedReUpsampledSamplingPoints1 >= nMatchedReUpsampledSamplingPoints2) ? iIter2 : iIter1);

                Particle particle(*(*iIterA), *(*iIterB));
                const LArPointingCluster pointingClusterA(pAlgorithm->GetCachedSlidingFitResult(particle.m_pClusterA));
                const LArPointingCluster pointingClusterB(pAlgorithm->GetCachedSlidingFitResult(particle.m_pClusterB));

                LArPointingCluster::Vertex vertexA, vertexB;
                LArPointingClusterHelper::GetClosestVertices(pointingClusterA, pointingClusterB, vertexA, vertexB);

                float transverseAB(std::numeric_limits<float>::max()), transverseBA(std::numeric_limits<float>::max());
                float longitudinalAB(-std::numeric_limits<float>::max()), longitudinalBA(-std::numeric_limits<float>::max());

                LArPointingClusterHelper::GetImpactParameters(vertexA, vertexB, longitudinalAB, transverseAB);
                LArPointingClusterHelper::GetImpactParameters(vertexB, vertexA, longitudinalBA, transverseBA);

                if (std::min(longitudinalAB, longitudinalBA) < m_minLongitudinalImpactParameter)
                    continue;

                if (std::min(transverseAB, transverseBA) > m_maxTransverseImpactParameter)
                    continue;

                const float cosTheta(-vertexA.GetDirection().GetCosOpeningAngle(vertexB.GetDirection()));

                if (cosTheta < m_minImpactParameterCosTheta)
                    continue;

                const bool isALowestInX(this->IsALowestInX(pointingClusterA, pointingClusterB));
                const CartesianVector splitPosition((vertexA.GetPosition() + vertexB.GetPosition()) * 0.5f);
                const bool isThreeDKink(this->IsThreeDKink(pAlgorithm, particle, splitPosition, isALowestInX));

                // Construct the modification object
                Modification modification;

                if (isThreeDKink)
                {
                    const TwoDSlidingFitResult &fitResultCommon(pAlgorithm->GetCachedSlidingFitResult(particle.m_pCommonCluster));

                    CartesianVector splitPositionCommon(0.f, 0.f, 0.f);
                    if (STATUS_CODE_SUCCESS != fitResultCommon.GetGlobalFitPositionAtX(splitPosition.GetX(), splitPositionCommon))
                    {
                        continue;
                    }

                    modification.m_splitPositionMap[particle.m_pCommonCluster].emplace_back(splitPositionCommon);
                }
                else
                {
                    const bool vertexAIsLowX(vertexA.GetPosition().GetX() < vertexB.GetPosition().GetX());
                    const Cluster *const pLowXCluster(vertexAIsLowX ? particle.m_pClusterA : particle.m_pClusterB);
                    const Cluster *const pHighXCluster(vertexAIsLowX ? particle.m_pClusterB : particle.m_pClusterA);
                    modification.m_clusterMergeMap[pLowXCluster].emplace_back(pHighXCluster);
                }

                modification.m_affectedClusters.emplace_back(particle.m_pClusterA);
                modification.m_affectedClusters.emplace_back(particle.m_pClusterB);
                modification.m_affectedClusters.emplace_back(particle.m_pCommonCluster);

                modificationList.emplace_back(modification);
            }
            catch (const StatusCodeException &statusCodeException)
            {
                if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                    throw statusCodeException;

                continue;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewThreeDKinkTool::IsThreeDKink(TwoViewTransverseTracksAlgorithm *const pAlgorithm, const Particle &particle,
    const CartesianVector &splitPosition, const bool isALowestInX) const
{
    try
    {
        const TwoDSlidingFitResult &fitResultCommon(pAlgorithm->GetCachedSlidingFitResult(particle.m_pCommonCluster));
        const TwoDSlidingFitResult &lowXFitResult(isALowestInX ? pAlgorithm->GetCachedSlidingFitResult(particle.m_pClusterA)
                                                               : pAlgorithm->GetCachedSlidingFitResult(particle.m_pClusterB));
        const TwoDSlidingFitResult &highXFitResult(isALowestInX ? pAlgorithm->GetCachedSlidingFitResult(particle.m_pClusterB)
                                                                : pAlgorithm->GetCachedSlidingFitResult(particle.m_pClusterA));

        const float minusX(this->GetXSamplingPoint(splitPosition, false, lowXFitResult, fitResultCommon));
        const float plusX(this->GetXSamplingPoint(splitPosition, true, highXFitResult, fitResultCommon));
        const float splitX(splitPosition.GetX());

        CartesianVector minus1(0.f, 0.f, 0.f), split1(0.f, 0.f, 0.f), plus1(0.f, 0.f, 0.f);
        CartesianVector minus2(0.f, 0.f, 0.f), split2(splitPosition), plus2(0.f, 0.f, 0.f);

        if ((STATUS_CODE_SUCCESS != fitResultCommon.GetGlobalFitPositionAtX(minusX, minus1)) ||
            (STATUS_CODE_SUCCESS != fitResultCommon.GetGlobalFitPositionAtX(splitX, split1)) ||
            (STATUS_CODE_SUCCESS != fitResultCommon.GetGlobalFitPositionAtX(plusX, plus1)) ||
            (STATUS_CODE_SUCCESS != lowXFitResult.GetGlobalFitPositionAtX(minusX, minus2)) ||
            (STATUS_CODE_SUCCESS != highXFitResult.GetGlobalFitPositionAtX(plusX, plus2)))
        {
            return true; // split by default
        }

        // Extract results
        const HitType hitType1(LArClusterHelper::GetClusterHitType(particle.m_pCommonCluster));
        const HitType hitType2(LArClusterHelper::GetClusterHitType(particle.m_pClusterA));

        CartesianVector minus(0.f, 0.f, 0.f), split(0.f, 0.f, 0.f), plus(0.f, 0.f, 0.f);
        float chi2Dummy(std::numeric_limits<float>::max());
        LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), hitType1, hitType2, minus1, minus2, minus, chi2Dummy);
        LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), hitType1, hitType2, split1, split2, split, chi2Dummy);
        LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), hitType1, hitType2, plus1, plus2, plus, chi2Dummy);

        // Apply final cuts
        const CartesianVector minusToSplit((split - minus).GetUnitVector());
        const CartesianVector splitToPlus((plus - split).GetUnitVector());
        const float dotProduct(minusToSplit.GetDotProduct(splitToPlus));

        if (dotProduct < m_cosThetaCutForKinkSearch)
            return true;
    }
    catch (const StatusCodeException &)
    {
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

TwoViewThreeDKinkTool::Particle::Particle(const MatrixType::Element &elementA, const MatrixType::Element &elementB)
{
    m_pClusterA = (elementA.GetCluster1() != elementB.GetCluster1()) ? elementA.GetCluster1() : elementA.GetCluster2();
    m_pClusterB = (elementA.GetCluster1() != elementB.GetCluster1()) ? elementB.GetCluster1() : elementB.GetCluster2();
    m_pCommonCluster = (elementA.GetCluster1() == elementB.GetCluster1()) ? elementA.GetCluster1() : elementA.GetCluster2();

    if (m_pClusterA == m_pClusterB)
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoViewThreeDKinkTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinMatchingScore", m_minMatchingScore));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinLocallyMatchedFraction", m_minLocallyMatchedFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinXOverlapFraction", m_minXOverlapFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinLongitudinalImpactParameter", m_minLongitudinalImpactParameter));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NLayersForKinkSearch", m_nLayersForKinkSearch));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "AdditionalXStepForKinkSearch", m_additionalXStepForKinkSearch));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxTransverseImpactParameter", m_maxTransverseImpactParameter));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinImpactParameterCosTheta", m_minImpactParameterCosTheta));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "CosThetaCutForKinkSearch", m_cosThetaCutForKinkSearch));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
