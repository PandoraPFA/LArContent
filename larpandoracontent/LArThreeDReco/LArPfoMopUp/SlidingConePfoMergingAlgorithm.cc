/**
 *  @file   larpandoracontent/LArThreeDReco/LArPfoMopUp/SlidingConePfoMergingAlgorithm.cc
 * 
 *  @brief  Implementation of the sliding cone pfo merging algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingConeFitResult.h"

#include "larpandoracontent/LArThreeDReco/LArPfoMopUp/SlidingConePfoMergingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

SlidingConePfoMergingAlgorithm::SlidingConePfoMergingAlgorithm() :
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

StatusCode SlidingConePfoMergingAlgorithm::Run()
{
    const Vertex *pVertex(nullptr);
    this->GetInteractionVertex(pVertex);

    if (m_useVertex && !pVertex)
    {
        std::cout << "SlidingConePfoMergingAlgorithm - interaction vertex not available for use." << std::endl;
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

void SlidingConePfoMergingAlgorithm::GetInteractionVertex(const Vertex *&pVertex) const
{
    if (!m_useVertex)
        return;

    const VertexList *pVertexList = nullptr;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetCurrentList(*this, pVertexList));

    pVertex = ((pVertexList && (pVertexList->size() == 1) && (VERTEX_3D == (*(pVertexList->begin()))->GetVertexType())) ? *(pVertexList->begin()) : nullptr);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SlidingConePfoMergingAlgorithm::GetThreeDClusters(ClusterVector &clusters3D, ClusterToPfoMap &clusterToPfoMap) const
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

void SlidingConePfoMergingAlgorithm::GetClusterMergeMap(const Vertex *const pVertex, const ClusterVector &clusters3D,
    const ClusterToPfoMap &clusterToPfoMap, ClusterMergeMap &clusterMergeMap) const
{
    const float layerPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));

    for (const Cluster *const pShowerCluster : clusters3D)
    {
        if (pShowerCluster->GetNCaloHits() < m_minHitsToConsider3DShower)
            continue;

        if (!LArPfoHelper::IsShower(clusterToPfoMap.at(pShowerCluster)))
            continue;

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
//if (boundedFraction1 > m_coneBoundedFraction1 && boundedFraction2 > m_coneBoundedFraction2)
//{
//std::cout << "ClosestDistance " << LArClusterHelper::GetClosestDistance(pShowerCluster, pNearbyCluster) << ", coneLength " << simpleCone.GetConeLength() << std::endl;
//std::cout << " boundedFraction1 " << boundedFraction1 << " boundedFraction2 " << boundedFraction2 << std::endl;
//ClusterList temp3; temp3.insert(pShowerCluster);
//ClusterList temp4; temp4.insert(pNearbyCluster);
//PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &temp3, "pShowerCluster", RED);
//PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &temp4, "pNearbyCluster", BLUE);
//PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &simpleCone.GetConeApex(), "coneApex", GRAY, 1);
//const CartesianVector maxBaseCentre(simpleCone.GetConeApex() + simpleCone.GetConeDirection() * coneLength);
//PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &maxBaseCentre, "coneBaseCentre", CYAN, 1);
//PandoraMonitoringApi::ViewEvent(this->GetPandora());
//}
                if (clusterMerge < bestClusterMerge)
                    bestClusterMerge = clusterMerge;
            }

            if ((bestClusterMerge.GetBoundedFraction1() > m_coneBoundedFraction1) && (bestClusterMerge.GetBoundedFraction2() > m_coneBoundedFraction2))
                clusterMergeMap[pNearbyCluster].push_back(bestClusterMerge);
        }
    }

    for (ClusterMergeMap::value_type &mapEntry : clusterMergeMap)
        std::sort(mapEntry.second.begin(), mapEntry.second.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SlidingConePfoMergingAlgorithm::MakePfoMerges(const ClusterToPfoMap &clusterToPfoMap, const ClusterMergeMap &clusterMergeMap) const
{
    ClusterVector daughterClusters;

    for (const ClusterMergeMap::value_type &mapEntry : clusterMergeMap)
        daughterClusters.push_back(mapEntry.first);

    std::sort(daughterClusters.begin(), daughterClusters.end(), LArClusterHelper::SortByNHits);

    bool pfosMerged(false);
    //ClusterReplacementMap clusterReplacementMap;

    for (ClusterVector::const_reverse_iterator rIter = daughterClusters.rbegin(), rIterEnd = daughterClusters.rend(); rIter != rIterEnd; ++rIter)
    {
std::cout << "FOR CONSIDERATION " << std::endl;
        const Cluster *const pDaughterCluster(*rIter);
        const ClusterMerge &bestClusterMerge(clusterMergeMap.at(pDaughterCluster).at(0));
        const Cluster *const pParentCluster(bestClusterMerge.GetParentCluster());

        const Pfo *const pDaughterPfo(clusterToPfoMap.at(pDaughterCluster));
        const Pfo *const pParentPfo(clusterToPfoMap.at(pParentCluster));
std::cout << " WILL MERGE boundedFraction1 " << bestClusterMerge.GetBoundedFraction1() << " boundedFraction2 " << bestClusterMerge.GetBoundedFraction2() << std::endl;
std::cout << "To Enlarge Cluster " << pParentCluster << std::endl;
std::cout << "To Delete Cluster " << pDaughterCluster << std::endl;
std::cout << "To Enlarge Pfo " << pParentPfo << std::endl;
std::cout << "To Delete Pfo " << pDaughterPfo << std::endl;
ClusterList temp3; temp3.insert(pParentCluster);
ClusterList temp4; temp4.insert(pDaughterCluster);
PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &temp3, "pParentCluster", RED);
PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &temp4, "pDaughterCluster", BLUE);
PfoList temp5; temp5.insert(pParentPfo);
PfoList temp6; temp6.insert(pDaughterPfo);
PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &temp5, "pParentPfo", RED);
PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &temp6, "pDaughterPfo", BLUE);
PandoraMonitoringApi::ViewEvent(this->GetPandora());

        // TODO
        this->MergeAndDeletePfos(pParentPfo, pDaughterPfo);
        pfosMerged = true;
std::cout << "MERGED " << std::endl;
std::cout << "Enlarged Cluster " << pParentCluster << std::endl;
std::cout << "Deleted Cluster " << pDaughterCluster << std::endl;
std::cout << "Enlarged Pfo " << pParentPfo << std::endl;
std::cout << "Deleted Pfo " << pDaughterPfo << std::endl;
return true;
    }

    return pfosMerged;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

bool SlidingConePfoMergingAlgorithm::ClusterMerge::operator<(const ClusterMerge &rhs) const
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

StatusCode SlidingConePfoMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
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

    m_daughterListNames.insert(m_daughterListNames.end(), m_inputPfoListNames.begin(), m_inputPfoListNames.end());

    return PfoMergingBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
