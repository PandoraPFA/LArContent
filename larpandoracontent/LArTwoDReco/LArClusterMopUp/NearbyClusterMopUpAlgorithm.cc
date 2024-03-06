/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterMopUp/NearbyClusterMopUpAlgorithm.cc
 *
 *  @brief  Implementation of the nearby cluster mop up algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterMopUp/NearbyClusterMopUpAlgorithm.h"

using namespace pandora;

namespace lar_content
{

NearbyClusterMopUpAlgorithm::NearbyClusterMopUpAlgorithm() :
    m_minHitsInCluster(5), m_vertexProximity(5.f), m_minClusterSeparation(2.5f), m_touchingDistance(0.001f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NearbyClusterMopUpAlgorithm::ClusterMopUp(const ClusterList &pfoClusters, const ClusterList &remnantClusters) const
{
    ClusterAssociationMap clusterAssociationMap;

    const VertexList *pVertexList(NULL);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pVertexList));
    const Vertex *const pVertex(
        ((pVertexList->size() == 1) && (VERTEX_3D == (*(pVertexList->begin()))->GetVertexType())) ? *(pVertexList->begin()) : NULL);

    ClusterVector sortedPfoClusters(pfoClusters.begin(), pfoClusters.end());
    std::sort(sortedPfoClusters.begin(), sortedPfoClusters.end(), LArClusterHelper::SortByNHits);

    ClusterVector sortedRemnantClusters(remnantClusters.begin(), remnantClusters.end());
    std::sort(sortedRemnantClusters.begin(), sortedRemnantClusters.end(), LArClusterHelper::SortByNHits);

    for (const Cluster *const pClusterP : sortedPfoClusters)
    {
        const HitType hitType(LArClusterHelper::GetClusterHitType(pClusterP));
        const CartesianVector vertexPosition2D(
            !pVertex ? CartesianVector(0.f, 0.f, 0.f) : LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), hitType));

        const float innerPV((vertexPosition2D - pClusterP->GetCentroid(pClusterP->GetInnerPseudoLayer())).GetMagnitude());
        const float outerPV((vertexPosition2D - pClusterP->GetCentroid(pClusterP->GetOuterPseudoLayer())).GetMagnitude());

        for (const Cluster *const pClusterR : sortedRemnantClusters)
        {
            if (pClusterR->GetNCaloHits() < m_minHitsInCluster)
                continue;

            const float innerRV((vertexPosition2D - pClusterR->GetCentroid(pClusterR->GetInnerPseudoLayer())).GetMagnitude());
            const float outerRV((vertexPosition2D - pClusterR->GetCentroid(pClusterR->GetOuterPseudoLayer())).GetMagnitude());

            // ATTN Could use pointing clusters here, for consistency with other vertex association mechanics
            if (pVertex && (((innerPV < m_vertexProximity) || (outerPV < m_vertexProximity)) &&
                               ((innerRV < m_vertexProximity) || (outerRV < m_vertexProximity))))
                continue;

            const float innerRP(LArClusterHelper::GetClosestDistance(pClusterR->GetCentroid(pClusterR->GetInnerPseudoLayer()), pClusterP));
            const float outerRP(LArClusterHelper::GetClosestDistance(pClusterR->GetCentroid(pClusterR->GetOuterPseudoLayer()), pClusterP));

            const float minSeparation(std::min(innerRP, outerRP));

            if (minSeparation > m_minClusterSeparation)
                continue;

            // ATTN Use of ClusterMopUp base algorithm assumes bigger figure of merit means better association
            if (m_touchingDistance < std::numeric_limits<float>::epsilon())
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

            AssociationDetails &associationDetails(clusterAssociationMap[pClusterR]);
            const float figureOfMerit((minSeparation < m_touchingDistance) ? 1.f / m_touchingDistance : 1.f / minSeparation);

            if (!associationDetails.insert(AssociationDetails::value_type(pClusterP, figureOfMerit)).second)
                throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
        }
    }

    this->MakeClusterMerges(clusterAssociationMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode NearbyClusterMopUpAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinHitsInCluster", m_minHitsInCluster));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexProximity", m_vertexProximity));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterSeparation", m_minClusterSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TouchingDistance", m_touchingDistance));

    return ClusterMopUpBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
