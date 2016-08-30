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
    m_nConeFitLayers(20),
    m_nConeFits(5),
    m_coneLengthMultiplier(7.f),
    m_maxConeLength(126.f),
    m_coneTanHalfAngle1(0.5f),
    m_coneBoundedFraction1(0.5f),
    m_coneTanHalfAngle2(0.75f),
    m_coneBoundedFraction2(0.75f)
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

    (void) this->MakeClusterMerges(clusterToPfoMap, clusterMergeMap);

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

                        if ((rL < 0.f) || (rL > coneLength))
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

bool SlidingConeClusterMopUpAlgorithm::MakeClusterMerges(const ClusterToPfoMap &clusterToPfoMap, const ClusterMergeMap &clusterMergeMap) const
{
    ClusterVector daughterClusters;
    for (const ClusterMergeMap::value_type &mapEntry : clusterMergeMap) daughterClusters.push_back(mapEntry.first);
    std::sort(daughterClusters.begin(), daughterClusters.end(), LArClusterHelper::SortByNHits);
typedef std::map<const Cluster*, ClusterList> HackMap;
HackMap hackMap;
for (ClusterVector::const_reverse_iterator rIter = daughterClusters.rbegin(), rIterEnd = daughterClusters.rend(); rIter != rIterEnd; ++rIter)
{
    const Cluster *const pDaughterCluster(*rIter);
    const HitType daughterHitType(LArClusterHelper::GetClusterHitType(pDaughterCluster));
    const Cluster *pParentCluster3D(clusterMergeMap.at(pDaughterCluster).at(0).GetParentCluster());
    const Cluster *pParentCluster(this->GetParentCluster(clusterToPfoMap.at(pParentCluster3D)->GetClusterList(), daughterHitType));
    hackMap[pParentCluster].insert(pDaughterCluster);
}
for (const auto &hackMapEntry : hackMap)
{
    std::cout << "MakeClusterMerges, nParents " << hackMap.size() << ", thisHitType " << (hackMapEntry.first ? LArClusterHelper::GetClusterHitType(hackMapEntry.first) : ECAL) << std::endl;
    ClusterList temp3; if (hackMapEntry.first) temp3.insert(hackMapEntry.first);
    ClusterList temp4; temp4.insert(hackMapEntry.second.begin(), hackMapEntry.second.end());
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &temp3, "ParentCluster", RED);
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &temp4, "DaughterClusters", BLUE);
    PandoraMonitoringApi::ViewEvent(this->GetPandora());
}
    bool clustersMerged(false);

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
            std::cout << "Adding missing cluster! " << std::endl;
PfoList pfoList; pfoList.insert(pParentPfo);
PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &pfoList, "pParentPfo", RED);
ClusterList clusterList; clusterList.insert(pDaughterCluster);
PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &clusterList, "pDaughterCluster", BLUE);
PandoraMonitoringApi::ViewEvent(this->GetPandora());
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pParentPfo, pDaughterCluster));
PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &pfoList, "pParentPfo", RED);
PandoraMonitoringApi::ViewEvent(this->GetPandora());
        }

        clustersMerged = true;
    }

    return clustersMerged;
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
