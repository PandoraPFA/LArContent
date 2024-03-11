/**
 *  @file   larpandoracontent/LArThreeDReco/LArShowerMatching/SplitShowersTool.cc
 *
 *  @brief  Implementation of the split showers tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"

#include "larpandoracontent/LArObjects/LArPointingCluster.h"

#include "larpandoracontent/LArThreeDReco/LArShowerMatching/SplitShowersTool.h"

using namespace pandora;

namespace lar_content
{

SplitShowersTool::SplitShowersTool() :
    m_nCommonClusters(2),
    m_minMatchedFraction(0.25f),
    m_minMatchedSamplingPoints(40),
    m_checkClusterProximities(true),
    m_maxClusterSeparation(25.f),
    m_checkClusterVertexRelations(true),
    m_minVertexLongitudinalDistance(-2.5f),
    m_maxVertexLongitudinalDistance(20.f),
    m_maxVertexTransverseDistance(1.5f),
    m_vertexAngularAllowance(3.f),
    m_maxVertexAssociations(1),
    m_checkClusterSplitPositions(false),
    m_vetoMergeXDifference(2.f),
    m_vetoMergeXOverlap(2.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SplitShowersTool::Run(ThreeViewShowersAlgorithm *const pAlgorithm, TensorType &overlapTensor)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    ClusterMergeMap clusterMergeMap;
    this->FindSplitShowers(pAlgorithm, overlapTensor, clusterMergeMap);

    return this->ApplyChanges(pAlgorithm, clusterMergeMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SplitShowersTool::FindSplitShowers(ThreeViewShowersAlgorithm *const pAlgorithm, const TensorType &overlapTensor, ClusterMergeMap &clusterMergeMap) const
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

            this->FindShowerMerges(pAlgorithm, iteratorList, usedClusters, clusterMergeMap);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SplitShowersTool::PassesElementCuts(TensorType::ElementList::const_iterator eIter, const ClusterSet &usedClusters) const
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

void SplitShowersTool::SelectTensorElements(TensorType::ElementList::const_iterator eIter, const TensorType::ElementList &elementList,
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

void SplitShowersTool::FindShowerMerges(ThreeViewShowersAlgorithm *const pAlgorithm, const IteratorList &iteratorList,
    ClusterSet &usedClusters, ClusterMergeMap &clusterMergeMap) const
{
    for (IteratorList::const_iterator iIter1 = iteratorList.begin(), iIter1End = iteratorList.end(); iIter1 != iIter1End; ++iIter1)
    {
        for (IteratorList::const_iterator iIter2 = iIter1; iIter2 != iIter1End; ++iIter2)
        {
            try
            {
                if (iIter1 == iIter2)
                    continue;

                const TensorType::Element &element1(*(*iIter1));
                const TensorType::Element &element2(*(*iIter2));

                ClusterList clusterListU(1, element1.GetClusterU());
                if (element1.GetClusterU() != element2.GetClusterU())
                    clusterListU.push_back(element2.GetClusterU());

                ClusterList clusterListV(1, element1.GetClusterV());
                if (element1.GetClusterV() != element2.GetClusterV())
                    clusterListV.push_back(element2.GetClusterV());

                ClusterList clusterListW(1, element1.GetClusterW());
                if (element1.GetClusterW() != element2.GetClusterW())
                    clusterListW.push_back(element2.GetClusterW());

                const unsigned int nClustersU(clusterListU.size()), nClustersV(clusterListV.size()), nClustersW(clusterListW.size());
                const unsigned int nClustersProduct(nClustersU * nClustersV * nClustersW);

                if (((1 == m_nCommonClusters) && (4 != nClustersProduct)) || ((2 == m_nCommonClusters) && (2 != nClustersProduct)))
                    throw StatusCodeException(STATUS_CODE_FAILURE);

                if ((1 == m_nCommonClusters) && !((2 == nClustersU) || (2 == nClustersV) || (2 == nClustersW)))
                    throw StatusCodeException(STATUS_CODE_FAILURE);

                if (m_checkClusterProximities &&
                    (!this->CheckClusterProximities(pAlgorithm, clusterListU) || !this->CheckClusterProximities(pAlgorithm, clusterListV) ||
                        !this->CheckClusterProximities(pAlgorithm, clusterListW)))
                {
                    continue;
                }

                if (m_checkClusterVertexRelations &&
                    (!this->CheckClusterVertexRelations(pAlgorithm, clusterListU) || !this->CheckClusterVertexRelations(pAlgorithm, clusterListV) ||
                        !this->CheckClusterVertexRelations(pAlgorithm, clusterListW)))
                {
                    continue;
                }

                if (m_checkClusterSplitPositions && !this->CheckClusterSplitPositions(pAlgorithm, clusterListU, clusterListV, clusterListW))
                {
                    continue;
                }

                this->SpecifyClusterMerges(pAlgorithm, clusterListU, clusterMergeMap);
                this->SpecifyClusterMerges(pAlgorithm, clusterListV, clusterMergeMap);
                this->SpecifyClusterMerges(pAlgorithm, clusterListW, clusterMergeMap);

                usedClusters.insert(clusterListU.begin(), clusterListU.end());
                usedClusters.insert(clusterListV.begin(), clusterListV.end());
                usedClusters.insert(clusterListW.begin(), clusterListW.end());
            }
            catch (StatusCodeException &)
            {
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SplitShowersTool::CheckClusterProximities(ThreeViewShowersAlgorithm * /*const pAlgorithm*/, const ClusterList &clusterList) const
{
    if (1 == clusterList.size())
        return true;

    if (2 != clusterList.size())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const Cluster *const pCluster1(*(clusterList.begin()));
    const Cluster *const pCluster2(*(++(clusterList.begin())));

    const float outer12(LArClusterHelper::GetClosestDistance(pCluster1->GetCentroid(pCluster1->GetOuterPseudoLayer()), pCluster2));
    const float outer21(LArClusterHelper::GetClosestDistance(pCluster2->GetCentroid(pCluster2->GetOuterPseudoLayer()), pCluster1));
    const float inner12(LArClusterHelper::GetClosestDistance(pCluster1->GetCentroid(pCluster1->GetInnerPseudoLayer()), pCluster2));
    const float inner21(LArClusterHelper::GetClosestDistance(pCluster2->GetCentroid(pCluster2->GetInnerPseudoLayer()), pCluster1));

    if ((outer12 > m_maxClusterSeparation) && (outer21 > m_maxClusterSeparation) && (inner12 > m_maxClusterSeparation) && (inner21 > m_maxClusterSeparation))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SplitShowersTool::CheckClusterVertexRelations(ThreeViewShowersAlgorithm *const pAlgorithm, const ClusterList &clusterList) const
{
    const VertexList *pVertexList(NULL);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*pAlgorithm, pVertexList));
    const Vertex *const pVertex(
        ((pVertexList->size() == 1) && (VERTEX_3D == (*(pVertexList->begin()))->GetVertexType())) ? *(pVertexList->begin()) : NULL);

    if (NULL == pVertex)
        return true;

    unsigned int nVertexAssociations(0);

    for (ClusterList::const_iterator iter = clusterList.begin(), iterEnd = clusterList.end(); iter != iterEnd; ++iter)
    {
        try
        {
            const HitType hitType(LArClusterHelper::GetClusterHitType(*iter));
            const CartesianVector vertex2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), hitType));
            const LArPointingCluster pointingCluster(pAlgorithm->GetCachedSlidingFitResult(*iter).GetShowerFitResult());

            if (LArPointingClusterHelper::IsNode(vertex2D, pointingCluster.GetInnerVertex(), m_minVertexLongitudinalDistance, m_maxVertexTransverseDistance) ||
                LArPointingClusterHelper::IsNode(vertex2D, pointingCluster.GetOuterVertex(), m_minVertexLongitudinalDistance, m_maxVertexTransverseDistance) ||
                LArPointingClusterHelper::IsEmission(vertex2D, pointingCluster.GetInnerVertex(), m_minVertexLongitudinalDistance,
                    m_maxVertexLongitudinalDistance, m_maxVertexTransverseDistance, m_vertexAngularAllowance) ||
                LArPointingClusterHelper::IsEmission(vertex2D, pointingCluster.GetOuterVertex(), m_minVertexLongitudinalDistance,
                    m_maxVertexLongitudinalDistance, m_maxVertexTransverseDistance, m_vertexAngularAllowance))
            {
                ++nVertexAssociations;
            }
        }
        catch (StatusCodeException &)
        {
        }
    }

    if (nVertexAssociations > m_maxVertexAssociations)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SplitShowersTool::CheckClusterSplitPositions(ThreeViewShowersAlgorithm *const pAlgorithm, const ClusterList &clusterListU,
    const ClusterList &clusterListV, const ClusterList &clusterListW) const
{
    const unsigned int nClustersU(clusterListU.size()), nClustersV(clusterListV.size()), nClustersW(clusterListW.size());
    const unsigned int nClustersProduct(nClustersU * nClustersV * nClustersW);

    if (2 == nClustersProduct)
        return true;

    const ClusterList &clusterList1((1 == nClustersU) ? clusterListV : clusterListU);
    const ClusterList &clusterList2((1 == nClustersU) ? clusterListW : (1 == nClustersV) ? clusterListW : clusterListV);

    if ((2 != clusterList1.size()) || (2 != clusterList2.size()))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    float splitXPosition1(0.f), overlapX1(0.f);
    this->GetSplitXDetails(pAlgorithm, *(clusterList1.begin()), *(++(clusterList1.begin())), splitXPosition1, overlapX1);

    float splitXPosition2(0.f), overlapX2(0.f);
    this->GetSplitXDetails(pAlgorithm, *(clusterList2.begin()), *(++(clusterList2.begin())), splitXPosition2, overlapX2);

    if ((std::fabs(splitXPosition1 - splitXPosition2) < m_vetoMergeXDifference) && (overlapX1 < m_vetoMergeXOverlap) && (overlapX2 < m_vetoMergeXOverlap))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SplitShowersTool::GetSplitXDetails(ThreeViewShowersAlgorithm *const pAlgorithm, const Cluster *const pClusterA,
    const Cluster *const pClusterB, float &splitXPosition, float &overlapX) const
{
    const TwoDSlidingFitResult &fitResultA(pAlgorithm->GetCachedSlidingFitResult(pClusterA).GetShowerFitResult());
    const TwoDSlidingFitResult &fitResultB(pAlgorithm->GetCachedSlidingFitResult(pClusterB).GetShowerFitResult());

    const float minXA(std::min(fitResultA.GetGlobalMinLayerPosition().GetX(), fitResultA.GetGlobalMaxLayerPosition().GetX()));
    const float maxXA(std::max(fitResultA.GetGlobalMinLayerPosition().GetX(), fitResultA.GetGlobalMaxLayerPosition().GetX()));
    const float minXB(std::min(fitResultB.GetGlobalMinLayerPosition().GetX(), fitResultB.GetGlobalMaxLayerPosition().GetX()));
    const float maxXB(std::max(fitResultB.GetGlobalMinLayerPosition().GetX(), fitResultB.GetGlobalMaxLayerPosition().GetX()));

    FloatVector floatVector;
    floatVector.push_back(minXA);
    floatVector.push_back(maxXA);
    floatVector.push_back(minXB);
    floatVector.push_back(maxXB);
    std::sort(floatVector.begin(), floatVector.end());

    if (4 != floatVector.size())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    splitXPosition = 0.5f * (floatVector.at(1) + floatVector.at(2));
    overlapX = std::max(0.f, std::min(maxXA, maxXB) - std::max(minXA, minXB));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SplitShowersTool::SpecifyClusterMerges(ThreeViewShowersAlgorithm *const pAlgorithm, const ClusterList &clusterList, ClusterMergeMap &clusterMergeMap) const
{
    if (1 == clusterList.size())
        return;

    if (2 != clusterList.size())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const Cluster *const pClusterA(*(clusterList.begin())), *const pClusterB(*(++(clusterList.begin())));
    const TwoDSlidingFitResult &fitResultA(pAlgorithm->GetCachedSlidingFitResult(pClusterA).GetShowerFitResult());
    const TwoDSlidingFitResult &fitResultB(pAlgorithm->GetCachedSlidingFitResult(pClusterB).GetShowerFitResult());

    const float minXA(std::min(fitResultA.GetGlobalMinLayerPosition().GetX(), fitResultA.GetGlobalMaxLayerPosition().GetX()));
    const float minXB(std::min(fitResultB.GetGlobalMinLayerPosition().GetX(), fitResultB.GetGlobalMaxLayerPosition().GetX()));

    const Cluster *const pLowXCluster((minXA < minXB) ? pClusterA : pClusterB);
    const Cluster *const pHighXCluster((minXA < minXB) ? pClusterB : pClusterA);
    clusterMergeMap[pLowXCluster].push_back(pHighXCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SplitShowersTool::ApplyChanges(ThreeViewShowersAlgorithm *const pAlgorithm, const ClusterMergeMap &clusterMergeMap) const
{
    ClusterMergeMap consolidatedMergeMap;

    ClusterList clusterList;
    for (const auto &mapEntry : clusterMergeMap)
        clusterList.push_back(mapEntry.first);
    clusterList.sort(LArClusterHelper::SortByNHits);

    for (const Cluster *const pParentCluster : clusterList)
    {
        const ClusterList &daughterClusters(clusterMergeMap.at(pParentCluster));

        for (const Cluster *const pDaughterCluster : daughterClusters)
        {
            if (consolidatedMergeMap.count(pDaughterCluster))
                throw StatusCodeException(STATUS_CODE_FAILURE);
        }

        ClusterList &targetClusterList(consolidatedMergeMap[pParentCluster]);
        targetClusterList.insert(targetClusterList.end(), daughterClusters.begin(), daughterClusters.end());
    }

    return pAlgorithm->MakeClusterMerges(consolidatedMergeMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SplitShowersTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NCommonClusters", m_nCommonClusters));

    if ((1 != m_nCommonClusters) && (2 != m_nCommonClusters))
    {
        std::cout << "SplitShowersTool: NCommonClusters must be set to either 1 or 2 (provided: " << m_nCommonClusters << ") " << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinMatchedFraction", m_minMatchedFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinMatchedSamplingPoints", m_minMatchedSamplingPoints));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "CheckClusterProximities", m_checkClusterProximities));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxClusterSeparation", m_maxClusterSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "CheckClusterVertexRelations", m_checkClusterVertexRelations));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinVertexLongitudinalDistance", m_minVertexLongitudinalDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxVertexLongitudinalDistance", m_maxVertexLongitudinalDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxVertexTransverseDistance", m_maxVertexTransverseDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexAngularAllowance", m_vertexAngularAllowance));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxVertexAssociations", m_maxVertexAssociations));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "CheckClusterSplitPositions", m_checkClusterSplitPositions));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VetoMergeXDifference", m_vetoMergeXDifference));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VetoMergeXOverlap", m_vetoMergeXOverlap));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
