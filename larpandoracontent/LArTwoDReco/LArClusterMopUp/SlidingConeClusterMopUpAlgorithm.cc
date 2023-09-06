/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterMopUp/SlidingConeClusterMopUpAlgorithm.cc
 *
 *  @brief  Implementation of the sliding cone cluster mop up algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingConeFitResult.h"

#include "larpandoracontent/LArTwoDReco/LArClusterMopUp/SlidingConeClusterMopUpAlgorithm.h"

using namespace pandora;

namespace lar_content
{

SlidingConeClusterMopUpAlgorithm::SlidingConeClusterMopUpAlgorithm() :
    m_useVertex(true),
    m_maxIterations(1000),
    m_maxHitsToConsider3DTrack(100),
    m_minHitsToConsider3DShower(20),
    m_maxHitsToConsider2DCluster(50),
    m_halfWindowLayers(20),
    m_nConeFitLayers(40),
    m_nConeFits(5),
    m_coneLengthMultiplier(3.f),
    m_maxConeLength(126.f),
    m_coneTanHalfAngle(0.2f),
    m_coneBoundedFraction(0.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SlidingConeClusterMopUpAlgorithm::Run()
{
    const Vertex *pVertex(nullptr);
    this->GetInteractionVertex(pVertex);

    if (m_useVertex && !pVertex)
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "SlidingConeClusterMopUpAlgorithm - interaction vertex not available for use." << std::endl;
        return STATUS_CODE_SUCCESS;
    }

    ClusterVector clusters3D;
    ClusterToPfoMap clusterToPfoMap;
    this->GetThreeDClusters(clusters3D, clusterToPfoMap);

    ClusterVector availableClusters2D;
    this->GetAvailableTwoDClusters(availableClusters2D);

    ClusterMergeMap clusterMergeMap;
    this->GetClusterMergeMap(pVertex, clusters3D, availableClusters2D, clusterMergeMap);

    this->MakeClusterMerges(clusterToPfoMap, clusterMergeMap);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SlidingConeClusterMopUpAlgorithm::GetInteractionVertex(const Vertex *&pVertex) const
{
    if (!m_useVertex)
        return;

    const VertexList *pVertexList = nullptr;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetCurrentList(*this, pVertexList));

    pVertex =
        ((pVertexList && (pVertexList->size() == 1) && (VERTEX_3D == (*(pVertexList->begin()))->GetVertexType())) ? *(pVertexList->begin()) : nullptr);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SlidingConeClusterMopUpAlgorithm::GetThreeDClusters(ClusterVector &clusters3D, ClusterToPfoMap &clusterToPfoMap) const
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

                if (LArPfoHelper::IsShower(pPfo) && (pCluster3D->GetNCaloHits() < m_minHitsToConsider3DShower))
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

void SlidingConeClusterMopUpAlgorithm::GetAvailableTwoDClusters(ClusterVector &availableClusters2D) const
{
    for (const std::string &clusterListName : m_daughterListNames)
    {
        const ClusterList *pClusterList(nullptr);

        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, clusterListName, pClusterList))
            continue;

        for (const Cluster *const pCluster : *pClusterList)
        {
            if (!PandoraContentApi::IsAvailable(*this, pCluster))
                continue;

            if (pCluster->GetNCaloHits() > m_maxHitsToConsider2DCluster)
                continue;

            availableClusters2D.push_back(pCluster);
        }
    }

    std::sort(availableClusters2D.begin(), availableClusters2D.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SlidingConeClusterMopUpAlgorithm::GetClusterMergeMap(const Vertex *const pVertex, const ClusterVector &clusters3D,
    const ClusterVector &availableClusters2D, ClusterMergeMap &clusterMergeMap) const
{
    for (const Cluster *const pShowerCluster : clusters3D)
    {
        float coneLength3D(0.f);
        SimpleConeList simpleConeList3D;

        try
        {
            HitType view{LArClusterHelper::GetClusterHitType(pShowerCluster)};
            if (!(view == TPC_VIEW_U || view == TPC_VIEW_V))
                view = TPC_VIEW_W;
            const float layerPitch(LArGeometryHelper::GetWirePitch(this->GetPandora(), view));
            const ThreeDSlidingConeFitResult slidingConeFitResult3D(pShowerCluster, m_halfWindowLayers, layerPitch);

            const CartesianVector &minLayerPosition(slidingConeFitResult3D.GetSlidingFitResult().GetGlobalMinLayerPosition());
            const CartesianVector &maxLayerPosition(slidingConeFitResult3D.GetSlidingFitResult().GetGlobalMaxLayerPosition());
            coneLength3D = std::min(m_coneLengthMultiplier * (maxLayerPosition - minLayerPosition).GetMagnitude(), m_maxConeLength);

            const float vertexToMinLayer(!pVertex ? 0.f : (pVertex->GetPosition() - minLayerPosition).GetMagnitude());
            const float vertexToMaxLayer(!pVertex ? 0.f : (pVertex->GetPosition() - maxLayerPosition).GetMagnitude());
            const ConeSelection coneSelection(
                !pVertex ? CONE_BOTH_DIRECTIONS : (vertexToMaxLayer > vertexToMinLayer) ? CONE_FORWARD_ONLY : CONE_BACKWARD_ONLY);

            slidingConeFitResult3D.GetSimpleConeList(m_nConeFitLayers, m_nConeFits, coneSelection, simpleConeList3D);
        }
        catch (const StatusCodeException &)
        {
            continue;
        }

        for (const Cluster *const pNearbyCluster2D : availableClusters2D)
        {
            ClusterMerge bestClusterMerge(nullptr, 0.f, 0.f);
            const HitType hitType(LArClusterHelper::GetClusterHitType(pNearbyCluster2D));

            for (const SimpleCone &simpleCone3D : simpleConeList3D)
            {
                const CartesianVector coneBaseCentre3D(simpleCone3D.GetConeApex() + simpleCone3D.GetConeDirection() * coneLength3D);
                const CartesianVector coneApex2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), simpleCone3D.GetConeApex(), hitType));
                const CartesianVector coneBaseCentre2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), coneBaseCentre3D, hitType));

                const CartesianVector apexToBase2D(coneBaseCentre2D - coneApex2D);
                const SimpleCone simpleCone2D(coneApex2D, apexToBase2D.GetUnitVector(), apexToBase2D.GetMagnitude(), m_coneTanHalfAngle);
                const ClusterMerge clusterMerge(
                    pShowerCluster, simpleCone2D.GetBoundedHitFraction(pNearbyCluster2D), simpleCone2D.GetMeanRT(pNearbyCluster2D));

                if (clusterMerge < bestClusterMerge)
                    bestClusterMerge = clusterMerge;
            }

            if (bestClusterMerge.GetParentCluster() && (bestClusterMerge.GetBoundedFraction() > m_coneBoundedFraction))
                clusterMergeMap[pNearbyCluster2D].push_back(bestClusterMerge);
        }
    }

    for (ClusterMergeMap::value_type &mapEntry : clusterMergeMap)
        std::sort(mapEntry.second.begin(), mapEntry.second.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SlidingConeClusterMopUpAlgorithm::MakeClusterMerges(const ClusterToPfoMap &clusterToPfoMap, const ClusterMergeMap &clusterMergeMap) const
{
    ClusterVector daughterClusters;
    for (const ClusterMergeMap::value_type &mapEntry : clusterMergeMap)
        daughterClusters.push_back(mapEntry.first);
    std::sort(daughterClusters.begin(), daughterClusters.end(), LArClusterHelper::SortByNHits);

    for (ClusterVector::const_reverse_iterator rIter = daughterClusters.rbegin(), rIterEnd = daughterClusters.rend(); rIter != rIterEnd; ++rIter)
    {
        const Cluster *const pDaughterCluster(*rIter);
        const HitType daughterHitType(LArClusterHelper::GetClusterHitType(pDaughterCluster));
        const Cluster *const pParentCluster3D(clusterMergeMap.at(pDaughterCluster).at(0).GetParentCluster());

        const Pfo *const pParentPfo(clusterToPfoMap.at(pParentCluster3D));
        const Cluster *const pParentCluster(this->GetParentCluster(pParentPfo->GetClusterList(), daughterHitType));

        if (pParentCluster)
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                PandoraContentApi::MergeAndDeleteClusters(
                    *this, pParentCluster, pDaughterCluster, this->GetListName(pParentCluster), this->GetListName(pDaughterCluster)));
        }
        else
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pParentPfo, pDaughterCluster));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

bool SlidingConeClusterMopUpAlgorithm::ClusterMerge::operator<(const ClusterMerge &rhs) const
{
    if (!this->GetParentCluster() && !rhs.GetParentCluster())
        return false;

    if (this->GetParentCluster() && !rhs.GetParentCluster())
        return true;

    if (!this->GetParentCluster() && rhs.GetParentCluster())
        return false;

    if (std::fabs(this->GetBoundedFraction() - rhs.GetBoundedFraction()) > std::numeric_limits<float>::epsilon())
        return (this->GetBoundedFraction() > rhs.GetBoundedFraction());

    if (std::fabs(this->GetMeanRT() - rhs.GetMeanRT()) > std::numeric_limits<float>::epsilon())
        return (this->GetMeanRT() < rhs.GetMeanRT());

    return LArClusterHelper::SortByNHits(this->GetParentCluster(), rhs.GetParentCluster());
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SlidingConeClusterMopUpAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputPfoListNames", m_inputPfoListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UseVertex", m_useVertex));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxIterations", m_maxIterations));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxHitsToConsider3DTrack", m_maxHitsToConsider3DTrack));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinHitsToConsider3DShower", m_minHitsToConsider3DShower));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxHitsToConsider2DCluster", m_maxHitsToConsider2DCluster));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingFitHalfWindow", m_halfWindowLayers));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NConeFitLayers", m_nConeFitLayers));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NConeFits", m_nConeFits));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ConeLengthMultiplier", m_coneLengthMultiplier));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxConeLength", m_maxConeLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ConeTanHalfAngle", m_coneTanHalfAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ConeBoundedFraction", m_coneBoundedFraction));

    return PfoMopUpBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
