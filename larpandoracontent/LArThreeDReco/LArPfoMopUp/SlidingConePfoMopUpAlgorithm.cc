/**
 *  @file   larpandoracontent/LArThreeDReco/LArPfoMopUp/SlidingConePfoMopUpAlgorithm.cc
 *
 *  @brief  Implementation of the sliding cone pfo mop up algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingConeFitResult.h"

#include "larpandoracontent/LArThreeDReco/LArPfoMopUp/SlidingConePfoMopUpAlgorithm.h"

using namespace pandora;

namespace lar_content
{

SlidingConePfoMopUpAlgorithm::SlidingConePfoMopUpAlgorithm() :
    m_useVertex(true),
    m_maxIterations(1000),
    m_maxHitsToConsider3DTrack(100),
    m_minHitsToConsider3DShower(20),
    m_halfWindowLayers(20),
    m_nConeFitLayers(20),
    m_nConeFits(5),
    m_coneLengthMultiplier(7.f),
    m_maxConeLength(126.f),
    m_coneTanHalfAngle1(0.5f),
    m_coneBoundedFraction1(0.5f),
    m_coneTanHalfAngle2(0.75f),
    m_coneBoundedFraction2(0.75f),
    m_minVertexLongitudinalDistance(-2.5f),
    m_maxVertexTransverseDistance(3.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SlidingConePfoMopUpAlgorithm::Run()
{
    const Vertex *pVertex(nullptr);
    this->GetInteractionVertex(pVertex);

    if (m_useVertex && !pVertex)
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "SlidingConePfoMopUpAlgorithm - interaction vertex not available for use." << std::endl;
        return STATUS_CODE_SUCCESS;
    }

    unsigned int nIterations(0);

    while (nIterations++ < m_maxIterations)
    {
        ClusterVector clusters3D;
        ClusterToPfoMap clusterToPfoMap;
        this->GetThreeDClusters(clusters3D, clusterToPfoMap);

        ClusterMergeMap clusterMergeMap;
        this->GetClusterMergeMap(pVertex, clusters3D, clusterToPfoMap, clusterMergeMap);

        if (!this->MakePfoMerges(clusterToPfoMap, clusterMergeMap))
            break;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SlidingConePfoMopUpAlgorithm::GetInteractionVertex(const Vertex *&pVertex) const
{
    if (!m_useVertex)
        return;

    const VertexList *pVertexList = nullptr;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetCurrentList(*this, pVertexList));

    pVertex =
        ((pVertexList && (pVertexList->size() == 1) && (VERTEX_3D == (*(pVertexList->begin()))->GetVertexType())) ? *(pVertexList->begin()) : nullptr);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SlidingConePfoMopUpAlgorithm::GetThreeDClusters(ClusterVector &clusters3D, ClusterToPfoMap &clusterToPfoMap) const
{
    for (const std::string &pfoListName : m_inputPfoListNames)
    {
        const PfoList *pPfoList(nullptr);

        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, pfoListName, pPfoList))
            continue;

        for (const Pfo *const pPfo : *pPfoList)
        {
            ClusterList pfoClusters3D;
            LArPfoHelper::GetThreeDClusterList(pPfo, pfoClusters3D);

            for (const Cluster *const pCluster3D : pfoClusters3D)
            {
                if (LArPfoHelper::IsTrack(pPfo) && (pCluster3D->GetNCaloHits() > m_maxHitsToConsider3DTrack))
                    continue;

                if (!clusterToPfoMap.insert(ClusterToPfoMap::value_type(pCluster3D, pPfo)).second)
                    throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);

                clusters3D.push_back(pCluster3D);
            }
        }
    }

    std::sort(clusters3D.begin(), clusters3D.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SlidingConePfoMopUpAlgorithm::GetClusterMergeMap(const Vertex *const pVertex, const ClusterVector &clusters3D,
    const ClusterToPfoMap &clusterToPfoMap, ClusterMergeMap &clusterMergeMap) const
{
    VertexAssociationMap vertexAssociationMap;
    const float pitchU{LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_U)};
    const float pitchV{LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_V)};
    const float pitchW{LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_W)};
    const float pitchMax{std::max({pitchU, pitchV, pitchW})};

    for (const Cluster *const pShowerCluster : clusters3D)
    {
        if ((pShowerCluster->GetNCaloHits() < m_minHitsToConsider3DShower) || !LArPfoHelper::IsShower(clusterToPfoMap.at(pShowerCluster)))
            continue;

        float coneLength(0.f);
        SimpleConeList simpleConeList;
        bool isShowerVertexAssociated(false);

        try
        {
            float layerPitch{0.f};
            const HitType view{LArClusterHelper::GetClusterHitType(pShowerCluster)};
            if (view == TPC_VIEW_U || view == TPC_VIEW_V || view == TPC_VIEW_W)
                layerPitch = LArGeometryHelper::GetWirePitch(this->GetPandora(), view);
            else
                layerPitch = pitchMax;

            const ThreeDSlidingConeFitResult slidingConeFitResult3D(pShowerCluster, m_halfWindowLayers, layerPitch);

            const CartesianVector &minLayerPosition(slidingConeFitResult3D.GetSlidingFitResult().GetGlobalMinLayerPosition());
            const CartesianVector &maxLayerPosition(slidingConeFitResult3D.GetSlidingFitResult().GetGlobalMaxLayerPosition());
            coneLength = std::min(m_coneLengthMultiplier * (maxLayerPosition - minLayerPosition).GetMagnitude(), m_maxConeLength);

            const float vertexToMinLayer(!pVertex ? 0.f : (pVertex->GetPosition() - minLayerPosition).GetMagnitude());
            const float vertexToMaxLayer(!pVertex ? 0.f : (pVertex->GetPosition() - maxLayerPosition).GetMagnitude());
            const ConeSelection coneSelection(
                !pVertex ? CONE_BOTH_DIRECTIONS : (vertexToMaxLayer > vertexToMinLayer) ? CONE_FORWARD_ONLY : CONE_BACKWARD_ONLY);

            slidingConeFitResult3D.GetSimpleConeList(m_nConeFitLayers, m_nConeFits, coneSelection, simpleConeList);
            isShowerVertexAssociated =
                this->IsVertexAssociated(pShowerCluster, pVertex, vertexAssociationMap, &(slidingConeFitResult3D.GetSlidingFitResult()));
        }
        catch (const StatusCodeException &)
        {
            continue;
        }

        for (const Cluster *const pNearbyCluster : clusters3D)
        {
            if (pNearbyCluster == pShowerCluster)
                continue;

            ClusterMerge bestClusterMerge(nullptr, 0.f, 0.f);

            for (const SimpleCone &simpleCone : simpleConeList)
            {
                const float boundedFraction1(simpleCone.GetBoundedHitFraction(pNearbyCluster, coneLength, m_coneTanHalfAngle1));
                const float boundedFraction2(simpleCone.GetBoundedHitFraction(pNearbyCluster, coneLength, m_coneTanHalfAngle2));
                const ClusterMerge clusterMerge(pShowerCluster, boundedFraction1, boundedFraction2);

                if (clusterMerge < bestClusterMerge)
                    bestClusterMerge = clusterMerge;
            }

            if (isShowerVertexAssociated && this->IsVertexAssociated(pNearbyCluster, pVertex, vertexAssociationMap))
                continue;

            if (bestClusterMerge.GetParentCluster() && (bestClusterMerge.GetBoundedFraction1() > m_coneBoundedFraction1) &&
                (bestClusterMerge.GetBoundedFraction2() > m_coneBoundedFraction2))
                clusterMergeMap[pNearbyCluster].push_back(bestClusterMerge);
        }
    }

    for (ClusterMergeMap::value_type &mapEntry : clusterMergeMap)
        std::sort(mapEntry.second.begin(), mapEntry.second.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SlidingConePfoMopUpAlgorithm::IsVertexAssociated(const Cluster *const pCluster, const Vertex *const pVertex,
    VertexAssociationMap &vertexAssociationMap, const ThreeDSlidingFitResult *const pSlidingFitResult) const
{
    if (!pVertex)
        return false;

    VertexAssociationMap::const_iterator iter = vertexAssociationMap.find(pCluster);

    if (vertexAssociationMap.end() != iter)
        return iter->second;

    const bool isVertexAssociated(this->IsVertexAssociated(pCluster, pVertex->GetPosition(), pSlidingFitResult));
    (void)vertexAssociationMap.insert(VertexAssociationMap::value_type(pCluster, isVertexAssociated));

    return isVertexAssociated;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SlidingConePfoMopUpAlgorithm::IsVertexAssociated(
    const Cluster *const pCluster, const CartesianVector &vertexPosition, const ThreeDSlidingFitResult *const pSlidingFitResult) const
{
    const float pitchU{LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_U)};
    const float pitchV{LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_V)};
    const float pitchW{LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_W)};
    const float pitchMax{std::max({pitchU, pitchV, pitchW})};

    try
    {
        float layerPitch{0.f};
        const HitType view{LArClusterHelper::GetClusterHitType(pCluster)};
        if (view == TPC_VIEW_U || view == TPC_VIEW_V || view == TPC_VIEW_W)
            layerPitch = LArGeometryHelper::GetWirePitch(this->GetPandora(), view);
        else
            layerPitch = pitchMax;

        const LArPointingCluster pointingCluster(
            pSlidingFitResult ? LArPointingCluster(*pSlidingFitResult) : LArPointingCluster(pCluster, m_halfWindowLayers, layerPitch));

        const bool useInner((pointingCluster.GetInnerVertex().GetPosition() - vertexPosition).GetMagnitudeSquared() <
                            (pointingCluster.GetOuterVertex().GetPosition() - vertexPosition).GetMagnitudeSquared());

        const LArPointingCluster::Vertex &daughterVertex(useInner ? pointingCluster.GetInnerVertex() : pointingCluster.GetOuterVertex());
        return LArPointingClusterHelper::IsNode(vertexPosition, daughterVertex, m_minVertexLongitudinalDistance, m_maxVertexTransverseDistance);
    }
    catch (const StatusCodeException &)
    {
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SlidingConePfoMopUpAlgorithm::MakePfoMerges(const ClusterToPfoMap &clusterToPfoMap, const ClusterMergeMap &clusterMergeMap) const
{
    ClusterVector daughterClusters;
    for (const ClusterMergeMap::value_type &mapEntry : clusterMergeMap)
        daughterClusters.push_back(mapEntry.first);
    std::sort(daughterClusters.begin(), daughterClusters.end(), LArClusterHelper::SortByNHits);

    bool pfosMerged(false);
    ClusterReplacementMap clusterReplacementMap;

    for (ClusterVector::const_reverse_iterator rIter = daughterClusters.rbegin(), rIterEnd = daughterClusters.rend(); rIter != rIterEnd; ++rIter)
    {
        const Cluster *const pDaughterCluster(*rIter);

        if (clusterReplacementMap.count(pDaughterCluster))
            throw StatusCodeException(STATUS_CODE_FAILURE);

        const Cluster *pParentCluster(clusterMergeMap.at(pDaughterCluster).at(0).GetParentCluster());

        if (clusterReplacementMap.count(pParentCluster))
            pParentCluster = clusterReplacementMap.at(pParentCluster);

        // ATTN Sign there was a reciprocal relationship in the cluster merge map (already actioned once)
        if (pDaughterCluster == pParentCluster)
            continue;

        // Key book-keeping on clusters and use cluster->pfo lookup
        const Pfo *const pDaughterPfo(clusterToPfoMap.at(pDaughterCluster));
        const Pfo *const pParentPfo(clusterToPfoMap.at(pParentCluster));
        this->MergeAndDeletePfos(pParentPfo, pDaughterPfo);
        pfosMerged = true;

        // Simple/placeholder book-keeping for reciprocal relationships and progressive merges
        clusterReplacementMap[pDaughterCluster] = pParentCluster;

        for (ClusterReplacementMap::value_type &mapEntry : clusterReplacementMap)
        {
            if (pDaughterCluster == mapEntry.second)
                mapEntry.second = pParentCluster;
        }
    }

    return pfosMerged;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

bool SlidingConePfoMopUpAlgorithm::ClusterMerge::operator<(const ClusterMerge &rhs) const
{
    if (!this->GetParentCluster() && !rhs.GetParentCluster())
        return false;

    if (this->GetParentCluster() && !rhs.GetParentCluster())
        return true;

    if (!this->GetParentCluster() && rhs.GetParentCluster())
        return false;

    if (std::fabs(this->GetBoundedFraction1() - rhs.GetBoundedFraction1()) > std::numeric_limits<float>::epsilon())
        return (this->GetBoundedFraction1() > rhs.GetBoundedFraction1());

    if (std::fabs(this->GetBoundedFraction2() - rhs.GetBoundedFraction2()) > std::numeric_limits<float>::epsilon())
        return (this->GetBoundedFraction2() > rhs.GetBoundedFraction2());

    return LArClusterHelper::SortByNHits(this->GetParentCluster(), rhs.GetParentCluster());
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SlidingConePfoMopUpAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputPfoListNames", m_inputPfoListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UseVertex", m_useVertex));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxIterations", m_maxIterations));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxHitsToConsider3DTrack", m_maxHitsToConsider3DTrack));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinHitsToConsider3DShower", m_minHitsToConsider3DShower));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingFitHalfWindow", m_halfWindowLayers));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NConeFitLayers", m_nConeFitLayers));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NConeFits", m_nConeFits));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ConeLengthMultiplier", m_coneLengthMultiplier));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxConeLength", m_maxConeLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ConeTanHalfAngle1", m_coneTanHalfAngle1));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ConeBoundedFraction1", m_coneBoundedFraction1));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ConeTanHalfAngle2", m_coneTanHalfAngle2));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ConeBoundedFraction2", m_coneBoundedFraction2));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinVertexLongitudinalDistance", m_minVertexLongitudinalDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxVertexTransverseDistance", m_maxVertexTransverseDistance));

    m_daughterListNames.insert(m_daughterListNames.end(), m_inputPfoListNames.begin(), m_inputPfoListNames.end());

    return PfoMopUpBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
