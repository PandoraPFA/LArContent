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

ProximityBasedMergingAlgorithm::ProximityBasedMergingAlgorithm()
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

            if (pClusterR->GetNCaloHits() < 5)
                continue;

            const float innerRV((vertexPosition2D - pClusterR->GetCentroid(pClusterR->GetInnerPseudoLayer())).GetMagnitude());
            const float outerRV((vertexPosition2D - pClusterR->GetCentroid(pClusterR->GetOuterPseudoLayer())).GetMagnitude());

            if (pVertex && (((innerPV < 5.f) || (outerPV < 5.f)) && ((innerRV < 5.f) || (outerRV < 5.f)))) // TODO use pointing clusters
                continue;

            const float innerRP(LArClusterHelper::GetClosestDistance(pClusterR->GetCentroid(pClusterR->GetInnerPseudoLayer()), pClusterP));
            const float outerRP(LArClusterHelper::GetClosestDistance(pClusterR->GetCentroid(pClusterR->GetOuterPseudoLayer()), pClusterP));

            const float minSeparation(std::min(innerRP, outerRP));

            if (minSeparation > 2.5f) // TODO
                continue;

            AssociationDetails &associationDetails(clusterAssociationMap[pClusterR]);
            const float figureOfMerit((minSeparation < 0.01f) ? std::numeric_limits<float>::max() : 1.f / minSeparation);

            if (!associationDetails.insert(AssociationDetails::value_type(pClusterP, figureOfMerit)).second)
                throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
        }
    }

    this->MakeClusterMerges(clusterAssociationMap, clusterToListNameMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ProximityBasedMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    return ClusterMopUpAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
