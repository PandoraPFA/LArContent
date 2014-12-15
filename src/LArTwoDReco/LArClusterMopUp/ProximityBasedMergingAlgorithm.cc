/**
 *  @file   LArContent/src/LArTwoDReco/LArClusterMopUp/ProximityBasedMergingAlgorithm.cc
 * 
 *  @brief  Implementation of the proximity based merging algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"

#include "LArTwoDReco/LArClusterMopUp/ProximityBasedMergingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ProximityBasedMergingAlgorithm::ProximityBasedMergingAlgorithm() :
    m_minHitsInCluster(5),
    m_vertexProximity(5.f),
    m_minClusterSeparation(2.5f),
    m_touchingDistance(0.001f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ProximityBasedMergingAlgorithm::ClusterMopUp(const ClusterList &pfoClusters, const ClusterList &remnantClusters,
    const ClusterToListNameMap &clusterToListNameMap) const
{
    ClusterAssociationMap clusterAssociationMap;

    const VertexList *pVertexList(NULL);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pVertexList));
    const Vertex *pVertex(((pVertexList->size() == 1) && (VERTEX_3D == (*(pVertexList->begin()))->GetVertexType())) ? *(pVertexList->begin()) : NULL);

    for (ClusterList::const_iterator pIter = pfoClusters.begin(), pIterEnd = pfoClusters.end(); pIter != pIterEnd; ++pIter)
    {
        Cluster *const pClusterP(*pIter);

        const HitType hitType(LArClusterHelper::GetClusterHitType(pClusterP));
        const CartesianVector vertexPosition2D(!pVertex ? CartesianVector(0.f, 0.f, 0.f) :
            LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), hitType));

        const float innerPV((vertexPosition2D - pClusterP->GetCentroid(pClusterP->GetInnerPseudoLayer())).GetMagnitude());
        const float outerPV((vertexPosition2D - pClusterP->GetCentroid(pClusterP->GetOuterPseudoLayer())).GetMagnitude());

        for (ClusterList::const_iterator rIter = remnantClusters.begin(), rIterEnd = remnantClusters.end(); rIter != rIterEnd; ++rIter)
        {
            Cluster *const pClusterR(*rIter);

            if (pClusterR->GetNCaloHits() < m_minHitsInCluster)
                continue;

            const float innerRV((vertexPosition2D - pClusterR->GetCentroid(pClusterR->GetInnerPseudoLayer())).GetMagnitude());
            const float outerRV((vertexPosition2D - pClusterR->GetCentroid(pClusterR->GetOuterPseudoLayer())).GetMagnitude());

            // ATTN Could use pointing clusters here, for consistency with other vertex association mechanics
            if (pVertex && (((innerPV < m_vertexProximity) || (outerPV < m_vertexProximity)) && ((innerRV < m_vertexProximity) || (outerRV < m_vertexProximity))))
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

    this->MakeClusterMerges(clusterAssociationMap, clusterToListNameMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ProximityBasedMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinHitsInCluster", m_minHitsInCluster));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VertexProximity", m_vertexProximity));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterSeparation", m_minClusterSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TouchingDistance", m_touchingDistance));

    return ClusterMopUpAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
