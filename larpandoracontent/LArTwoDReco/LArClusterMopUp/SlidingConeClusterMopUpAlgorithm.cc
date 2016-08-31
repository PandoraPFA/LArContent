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

// TODO Lots of refactoring!

SlidingConeClusterMopUpAlgorithm::SlidingConeClusterMopUpAlgorithm() :
    m_useVertex(true),
    m_maxIterations(1000),
    m_maxHitsToConsider3DTrack(100),
    m_minHitsToConsider3DShower(20),
    m_halfWindowLayers(20),
    m_nConeFitLayers(40),
    m_nConeFits(5),
    m_coneLengthMultiplier(3.f),
    m_maxConeLength(126.f),
    m_coneTanHalfAngle1(0.f),
    m_coneBoundedFraction1(-1.f),
    m_coneTanHalfAngle2(0.2f),
    m_coneBoundedFraction2(0.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SlidingConeClusterMopUpAlgorithm::Run()
{
    const Vertex *pVertex(nullptr);
    this->GetInteractionVertex(pVertex);

    if (m_useVertex && !pVertex)
    {
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

    pVertex = ((pVertexList && (pVertexList->size() == 1) && (VERTEX_3D == (*(pVertexList->begin()))->GetVertexType())) ? *(pVertexList->begin()) : nullptr);
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

            availableClusters2D.push_back(pCluster);
        }
    }

    std::sort(availableClusters2D.begin(), availableClusters2D.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SlidingConeClusterMopUpAlgorithm::GetClusterMergeMap(const Vertex *const pVertex, const ClusterVector &clusters3D,
    const ClusterVector &availableClusters2D, ClusterMergeMap &clusterMergeMap) const
{
    const float layerPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));

    for (const Cluster *const pShowerCluster : clusters3D)
    {
        float coneLength(0.f);
        SimpleConeList simpleConeList;

        try
        {
            const ThreeDSlidingConeFitResult slidingConeFitResult3D(pShowerCluster, m_halfWindowLayers, layerPitch);

            const CartesianVector &minLayerPosition(slidingConeFitResult3D.GetSlidingFitResult().GetGlobalMinLayerPosition());
            const CartesianVector &maxLayerPosition(slidingConeFitResult3D.GetSlidingFitResult().GetGlobalMaxLayerPosition());
            coneLength = std::min(m_coneLengthMultiplier * (maxLayerPosition - minLayerPosition).GetMagnitude(), m_maxConeLength);

            const float vertexToMinLayer(!pVertex ? 0.f : (pVertex->GetPosition() - minLayerPosition).GetMagnitude());
            const float vertexToMaxLayer(!pVertex ? 0.f : (pVertex->GetPosition() - maxLayerPosition).GetMagnitude());
            const ConeSelection coneSelection(!pVertex ? CONE_BOTH_DIRECTIONS : (vertexToMaxLayer > vertexToMinLayer) ? CONE_FORWARD_ONLY : CONE_BACKWARD_ONLY);

            slidingConeFitResult3D.GetSimpleConeList(m_nConeFitLayers, m_nConeFits, coneSelection, simpleConeList);
        }
        catch (const StatusCodeException &) {continue;}

        // TODO - refactor
        for (const Cluster *const pNearbyCluster : availableClusters2D)
        {
            const HitType hitType(LArClusterHelper::GetClusterHitType(pNearbyCluster));
            ClusterMerge bestClusterMerge(nullptr, 0.f, 0.f, 0.f);

            for (const SimpleCone &simpleCone : simpleConeList)
            {
                // TODO quick way of avoiding calculation
                const CartesianVector &coneApex3D(simpleCone.GetConeApex());
                const CartesianVector coneBaseCentre3D(simpleCone.GetConeApex() + simpleCone.GetConeDirection() * coneLength);

                const CartesianVector coneApex2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), coneApex3D, hitType));
                const CartesianVector coneBaseCentre2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), coneBaseCentre3D, hitType));

                const CartesianVector coneDirection2D((coneBaseCentre2D - coneApex2D).GetUnitVector());
                const float coneLength2D((coneBaseCentre2D - coneApex2D).GetMagnitude());

                float rTSum(0.f);
                unsigned int nMatchedHits1(0), nMatchedHits2(0);
                const unsigned int nClusterHits(pNearbyCluster->GetNCaloHits());

                const OrderedCaloHitList &orderedCaloHitList(pNearbyCluster->GetOrderedCaloHitList());

                for (const auto &mapEntry : orderedCaloHitList)
                {
                    for (const CaloHit *pCaloHit : *(mapEntry.second))
                    {
                        if (hitType != pCaloHit->GetHitType())
                            throw StatusCodeException(STATUS_CODE_FAILURE);

                        const CartesianVector displacement(pCaloHit->GetPositionVector() - coneApex2D);
                        const float rL(displacement.GetDotProduct(coneDirection2D));
                        const float rT(displacement.GetCrossProduct(coneDirection2D).GetMagnitude());
                        rTSum += rT;

                        if ((rL < 0.f) || (rL > coneLength2D))
                            continue;

                        if (rL * m_coneTanHalfAngle1 > rT)
                            ++nMatchedHits1;

                        if (rL * m_coneTanHalfAngle2 > rT)
                            ++nMatchedHits2;
                    }
                }

                const float boundedFraction1((nClusterHits > 0) ? static_cast<float>(nMatchedHits1) / static_cast<float>(nClusterHits) : 0.f);
                const float boundedFraction2((nClusterHits > 0) ? static_cast<float>(nMatchedHits2) / static_cast<float>(nClusterHits) : 0.f);
                const float meanRT((nClusterHits > 0) ? rTSum / static_cast<float>(nClusterHits) : 0.f);

                const ClusterMerge clusterMerge(pShowerCluster, boundedFraction1, boundedFraction2, meanRT);

                if (clusterMerge < bestClusterMerge)
                    bestClusterMerge = clusterMerge;
            }

            if (bestClusterMerge.GetParentCluster() && (bestClusterMerge.GetBoundedFraction1() > m_coneBoundedFraction1) && (bestClusterMerge.GetBoundedFraction2() > m_coneBoundedFraction2))
                clusterMergeMap[pNearbyCluster].push_back(bestClusterMerge);
        }
    }

    for (ClusterMergeMap::value_type &mapEntry : clusterMergeMap)
        std::sort(mapEntry.second.begin(), mapEntry.second.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SlidingConeClusterMopUpAlgorithm::MakeClusterMerges(const ClusterToPfoMap &clusterToPfoMap, const ClusterMergeMap &clusterMergeMap) const
{
    ClusterVector daughterClusters;
    for (const ClusterMergeMap::value_type &mapEntry : clusterMergeMap) daughterClusters.push_back(mapEntry.first);
    std::sort(daughterClusters.begin(), daughterClusters.end(), LArClusterHelper::SortByNHits);

//typedef std::map<const Cluster*, ClusterList> ClusterToClusterListMap; ClusterToClusterListMap parentTo2DDaughters;
//typedef std::map<const Cluster*, const Cluster*> ClusterToClusterMap; ClusterToClusterMap parent2DToParent3D;
//for (ClusterVector::const_reverse_iterator rIter = daughterClusters.rbegin(), rIterEnd = daughterClusters.rend(); rIter != rIterEnd; ++rIter)
//{
//    const Cluster *const pDaughterCluster(*rIter);
//    const HitType daughterHitType(LArClusterHelper::GetClusterHitType(pDaughterCluster));
//    const Cluster *pParentCluster3D(clusterMergeMap.at(pDaughterCluster).at(0).GetParentCluster());
//    const Cluster *pParentCluster(this->GetParentCluster(clusterToPfoMap.at(pParentCluster3D)->GetClusterList(), daughterHitType));
//    parentTo2DDaughters[pParentCluster].insert(pDaughterCluster);
//    parent2DToParent3D[pParentCluster] = pParentCluster3D;
//}
//for (const auto &mapEntry : parentTo2DDaughters)
//{
//    const Cluster *const pParentCluster3d(parent2DToParent3D.at(mapEntry.first));
//    const float layerPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
//    const ThreeDSlidingConeFitResult slidingConeFitResult3D(pParentCluster3d, m_halfWindowLayers, layerPitch);
//
//    const CartesianVector &minLayerPosition(slidingConeFitResult3D.GetSlidingFitResult().GetGlobalMinLayerPosition());
//    const CartesianVector &maxLayerPosition(slidingConeFitResult3D.GetSlidingFitResult().GetGlobalMaxLayerPosition());
//    const float coneLength = std::min(m_coneLengthMultiplier * (maxLayerPosition - minLayerPosition).GetMagnitude(), m_maxConeLength);
//
//    const Vertex *pVertex(nullptr);
//    this->GetInteractionVertex(pVertex);
//    const float vertexToMinLayer(!pVertex ? 0.f : (pVertex->GetPosition() - minLayerPosition).GetMagnitude());
//    const float vertexToMaxLayer(!pVertex ? 0.f : (pVertex->GetPosition() - maxLayerPosition).GetMagnitude());
//    const ConeSelection coneSelection(!pVertex ? CONE_BOTH_DIRECTIONS : (vertexToMaxLayer > vertexToMinLayer) ? CONE_FORWARD_ONLY : CONE_BACKWARD_ONLY);
//
//    SimpleConeList simpleConeList;
//    slidingConeFitResult3D.GetSimpleConeList(m_nConeFitLayers, m_nConeFits, coneSelection, simpleConeList);
//
//    const TwoDSlidingFitResult &fitResult1(slidingConeFitResult3D.GetSlidingFitResult().GetFirstFitResult());
//    const TwoDSlidingFitResult &fitResult2(slidingConeFitResult3D.GetSlidingFitResult().GetSecondFitResult());
//
//    const LayerFitContributionMap &contributionMap1(fitResult1.GetLayerFitContributionMap());
//    const LayerFitContributionMap &contributionMap2(fitResult2.GetLayerFitContributionMap());
//
//    const HitType hitType(mapEntry.first ? LArClusterHelper::GetClusterHitType(mapEntry.first) : LArClusterHelper::GetClusterHitType(*(mapEntry.second.begin())));
//
//    ClusterList temp3; if (mapEntry.first) temp3.insert(mapEntry.first);
//    ClusterList temp4; temp4.insert(mapEntry.second.begin(), mapEntry.second.end());
//    ClusterList temp5; temp5.insert(pParentCluster3d);
//
//    for (const SimpleCone &simpleCone : simpleConeList)
//    {
//break;
//        const CartesianVector &coneApex3D(simpleCone.GetConeApex());
//        const CartesianVector coneBaseCentre3D(simpleCone.GetConeApex() + simpleCone.GetConeDirection() * coneLength);
//        const CartesianVector baseDir3D(simpleCone.GetConeDirection().GetCrossProduct(CartesianVector(0.f, 1.f, 0.f)).GetUnitVector());
//        const float baseDisp3D(coneLength * m_coneTanHalfAngle2);
//        const CartesianVector baseMarkerA3D(coneBaseCentre3D + (baseDir3D * baseDisp3D));
//        const CartesianVector baseMarkerB3D(coneBaseCentre3D + (baseDir3D * (-1.f * baseDisp3D)));
//
//        const CartesianVector coneApex2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), coneApex3D, hitType));
//        const CartesianVector coneBaseCentre2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), coneBaseCentre3D, hitType));
//
//        const CartesianVector coneDirection2D((coneBaseCentre2D - coneApex2D).GetUnitVector());
//        const float coneLength2D((coneBaseCentre2D - coneApex2D).GetMagnitude());
//        const CartesianVector baseDir2D(coneDirection2D.GetCrossProduct(CartesianVector(0.f, 1.f, 0.f)).GetUnitVector());
//        const float baseDisp2D(coneLength2D * m_coneTanHalfAngle2);
//        const CartesianVector baseMarkerA2D(coneBaseCentre2D + (baseDir2D * baseDisp2D));
//        const CartesianVector baseMarkerB2D(coneBaseCentre2D + (baseDir2D * (-1.f * baseDisp2D)));
//
//        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &minLayerPosition, "minLayerPosition", GRAY, 1);
//        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &maxLayerPosition, "maxLayerPosition", GRAY, 1);
//        const int nSteps(static_cast<int>((maxLayerPosition - minLayerPosition).GetMagnitude() / layerPitch));
//        for (int iStep = 0; iStep < nSteps; ++iStep)
//        {
//            const float rL((static_cast<float>(iStep) + 0.5f) * layerPitch);
//
//            CartesianVector fitPosition3D(0.f, 0.f, 0.f);
//            if (STATUS_CODE_SUCCESS != slidingConeFitResult3D.GetSlidingFitResult().GetGlobalFitPosition(rL, fitPosition3D))
//                 continue;
//
//            if (!contributionMap1.count(fitResult1.GetLayer(rL)) && !contributionMap2.count(fitResult2.GetLayer(rL)))
//                continue;
//
//            const CartesianVector axis(slidingConeFitResult3D.GetSlidingFitResult().GetAxisIntercept() + slidingConeFitResult3D.GetSlidingFitResult().GetAxisDirection() * rL);
//            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &axis, "axis3D", GRAY, 1);
//            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &fitPosition3D, "fitPosition3D", ORANGE, 1);
//        }
//
//        std::cout << "MakeClusterMerges, nParents " << parentTo2DDaughters.size() << ", thisHitType " << (mapEntry.first ? LArClusterHelper::GetClusterHitType(mapEntry.first) : ECAL) << std::endl;
//        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &temp3, "ParentCluster2D", RED);
//        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &temp5, "ParentCluster3D", GREEN);
//        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &temp4, "DaughterClusters", BLUE);
//        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &coneApex3D, "coneApex3D", CYAN, 2);
//        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &coneBaseCentre3D, "coneBaseCentre3D", CYAN, 2);
//        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &baseMarkerA3D, "baseMarkerA3D", CYAN, 1);
//        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &baseMarkerB3D, "baseMarkerB3D", CYAN, 1);
//        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &coneApex2D, "coneApex2D", MAGENTA, 1);
//        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &coneBaseCentre2D, "coneBaseCentre2D", MAGENTA, 1);
//        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &baseMarkerA2D, "baseMarkerA2D", MAGENTA, 1);
//        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &baseMarkerB2D, "baseMarkerB2D", MAGENTA, 1);
//        PandoraMonitoringApi::ViewEvent(this->GetPandora());
//    }
////continue;
//    std::cout << "MakeClusterMerges, nParents " << parentTo2DDaughters.size() << ", thisHitType " << (mapEntry.first ? LArClusterHelper::GetClusterHitType(mapEntry.first) : ECAL) << std::endl;
//    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &temp3, "ParentCluster2D", RED);
//    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &temp5, "ParentCluster3D", GREEN);
//    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &temp4, "DaughterClusters", BLUE);
//    PandoraMonitoringApi::ViewEvent(this->GetPandora());
//}

    for (ClusterVector::const_reverse_iterator rIter = daughterClusters.rbegin(), rIterEnd = daughterClusters.rend(); rIter != rIterEnd; ++rIter)
    {
        const Cluster *const pDaughterCluster(*rIter);
        const HitType daughterHitType(LArClusterHelper::GetClusterHitType(pDaughterCluster));
        const Cluster *const pParentCluster3D(clusterMergeMap.at(pDaughterCluster).at(0).GetParentCluster());

        const Pfo *const pParentPfo(clusterToPfoMap.at(pParentCluster3D));
        const Cluster *const pParentCluster(this->GetParentCluster(pParentPfo->GetClusterList(), daughterHitType));

        if (pParentCluster)
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pParentCluster, pDaughterCluster,
                this->GetListName(pParentCluster), this->GetListName(pDaughterCluster)));
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

    if (std::fabs(this->GetBoundedFraction1() - rhs.GetBoundedFraction1()) > std::numeric_limits<float>::epsilon())
        return (this->GetBoundedFraction1() > rhs.GetBoundedFraction1());

    if (std::fabs(this->GetBoundedFraction2() - rhs.GetBoundedFraction2()) > std::numeric_limits<float>::epsilon())
        return (this->GetBoundedFraction2() > rhs.GetBoundedFraction2());

    if (std::fabs(this->GetMeanRT() - rhs.GetMeanRT()) > std::numeric_limits<float>::epsilon())
        return (this->GetMeanRT() < rhs.GetMeanRT());

    return LArClusterHelper::SortByNHits(this->GetParentCluster(), rhs.GetParentCluster());
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SlidingConeClusterMopUpAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "InputPfoListNames", m_inputPfoListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseVertex", m_useVertex));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxIterations", m_maxIterations));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxHitsToConsider3DTrack", m_maxHitsToConsider3DTrack));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinHitsToConsider3DShower", m_minHitsToConsider3DShower));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitHalfWindow", m_halfWindowLayers));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NConeFitLayers", m_nConeFitLayers));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NConeFits", m_nConeFits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ConeLengthMultiplier", m_coneLengthMultiplier));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxConeLength", m_maxConeLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ConeTanHalfAngle1", m_coneTanHalfAngle1));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ConeBoundedFraction1", m_coneBoundedFraction1));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ConeTanHalfAngle2", m_coneTanHalfAngle2));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ConeBoundedFraction2", m_coneBoundedFraction2));

    return PfoMopUpBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
