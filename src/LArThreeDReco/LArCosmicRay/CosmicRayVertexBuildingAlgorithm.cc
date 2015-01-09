/**
 *  @file   LArContent/src/LArThreeDReco/LArCosmicRay/CosmicRayVertexBuildingAlgorithm.cc
 *
 *  @brief  Implementation of the cosmic-ray vertex building algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArGeometryHelper.h"
#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArPfoHelper.h"
#include "LArPlugins/LArTransformationPlugin.h"

#include "LArThreeDReco/LArCosmicRay/CosmicRayVertexBuildingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CosmicRayVertexBuildingAlgorithm::CosmicRayVertexBuildingAlgorithm() :
    m_useParentShowerVertex(false),
    m_halfWindowLayers(30)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayVertexBuildingAlgorithm::Run()
{
    const PfoList *pPfoList = NULL;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_parentPfoListName,
        pPfoList));

    if (NULL == pPfoList)
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "CosmicRayVertexBuildingAlgorithm: pfo list unavailable." << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    PfoVector pfoVector;
    LArPointingClusterMap pointingClusterMap;

    this->GetCosmicPfos(pPfoList, pfoVector);
    this->BuildPointingClusterMap(pfoVector, pointingClusterMap);
    this->BuildCosmicRayParticles(pointingClusterMap, pfoVector);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayVertexBuildingAlgorithm::GetCosmicPfos(const PfoList *const pPfoList, PfoVector &pfoVector) const
{
    PfoList outputList;

    for (PfoList::const_iterator pIter = pPfoList->begin(), pIterEnd = pPfoList->end(); pIter != pIterEnd; ++pIter)
    {
        if (!LArPfoHelper::IsFinalState(*pIter))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        LArPfoHelper::GetAllDownstreamPfos(*pIter, outputList);
    }

    for (PfoList::const_iterator pIter = outputList.begin(), pIterEnd = outputList.end(); pIter != pIterEnd; ++pIter)
    {
        ClusterList clusterList;
        LArPfoHelper::GetClusters(*pIter, TPC_3D, clusterList);

        if (clusterList.empty())
            continue;

        pfoVector.push_back(*pIter);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayVertexBuildingAlgorithm::BuildPointingClusterMap(const PfoVector &pfoVector, LArPointingClusterMap &pointingClusterMap) const
{
    const float slidingFitPitch(LArGeometryHelper::GetLArTransformationPlugin(this->GetPandora())->GetWireZPitch());

    for (PfoVector::const_iterator pIter = pfoVector.begin(), pIterEnd = pfoVector.end(); pIter != pIterEnd; ++pIter)
    {
        const ParticleFlowObject *const pPfo = *pIter;

        if (!LArPfoHelper::IsTrack(pPfo))
            continue;

        ClusterList clusterList;
        LArPfoHelper::GetClusters(pPfo, TPC_3D, clusterList);

        for (ClusterList::const_iterator cIter = clusterList.begin(), cIterEnd = clusterList.end(); cIter != cIterEnd; ++cIter)
        {
            const Cluster *const pCluster = *cIter;

            try
            {
                const LArPointingCluster pointingCluster(pCluster, m_halfWindowLayers, slidingFitPitch);

                if (!pointingClusterMap.insert(LArPointingClusterMap::value_type(pCluster, pointingCluster)).second)
                    throw StatusCodeException(STATUS_CODE_FAILURE);
            }
            catch (StatusCodeException &statusCodeException)
            {
                if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                    throw statusCodeException;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayVertexBuildingAlgorithm::BuildCosmicRayParticles(const LArPointingClusterMap &pointingClusterMap, const PfoVector &pfoVector) const
{
    for (PfoVector::const_iterator iter = pfoVector.begin(), iterEnd = pfoVector.end(); iter != iterEnd; ++iter)
    {
        ParticleFlowObject *const pPfo = *iter;

        if (LArPfoHelper::IsFinalState(pPfo))
        {
            this->BuildCosmicRayParent(pointingClusterMap, pPfo);
        }
        else
        {
            this->BuildCosmicRayDaughter(pPfo);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayVertexBuildingAlgorithm::BuildCosmicRayParent(const LArPointingClusterMap &pointingClusterMap, ParticleFlowObject *const pPfo) const
{
    ClusterList clusterList;
    LArPfoHelper::GetClusters(pPfo, TPC_3D, clusterList);

    if (clusterList.empty())
        return;

    // Take highest point as vertex of parent Pfos (TODO: do something more sophisticated for horizontal events)
    bool foundVtx(false);
    CartesianVector vtxPosition(0.f, 0.f, 0.f);
    CartesianVector vtxDirection(0.f, 0.f, 0.f);

    bool foundEnd(false);
    CartesianVector endPosition(0.f, 0.f, 0.f);
    CartesianVector endDirection(0.f, 0.f, 0.f);

    for (ClusterList::const_iterator cIter1 = clusterList.begin(), cIterEnd1 = clusterList.end(); cIter1 != cIterEnd1; ++cIter1)
    {
        const Cluster *const pCluster = *cIter1;

        try
        {
            CartesianVector minPosition(0.f, 0.f, 0.f), maxPosition(0.f,0.f,0.f);
            CartesianVector minDirection(0.f, 0.f, 0.f), maxDirection(0.f,0.f,0.f);

            LArPointingClusterMap::const_iterator cIter2 = pointingClusterMap.find(pCluster);

            if (pointingClusterMap.end() != cIter2)
            {
                const LArPointingCluster &pointingCluster(cIter2->second);

                minPosition = pointingCluster.GetInnerVertex().GetPosition();
                maxPosition = pointingCluster.GetOuterVertex().GetPosition();
                minDirection = pointingCluster.GetInnerVertex().GetDirection();
                maxDirection = pointingCluster.GetOuterVertex().GetDirection();
            }
            else
            {
                LArClusterHelper::GetExtremalCoordinates(pCluster, minPosition, maxPosition);
                minDirection = (maxPosition - minPosition).GetUnitVector();
                maxDirection = (minPosition - maxPosition).GetUnitVector();
            }

            if ((maxPosition - minPosition).GetMagnitudeSquared() < std::numeric_limits<float>::epsilon())
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);

            if (!foundVtx || (minPosition.GetY() > std::max(maxPosition.GetY(), vtxPosition.GetY())))
            {
                foundVtx = true;
                vtxPosition = minPosition;
                vtxDirection = minDirection;
            }

            if (!foundVtx || (maxPosition.GetY() > std::max(minPosition.GetY(), vtxPosition.GetY())))
            {
                foundVtx = true;
                vtxPosition = maxPosition;
                vtxDirection = maxDirection;
            }

            if (!foundEnd || (minPosition.GetY() < std::min(maxPosition.GetY(), endPosition.GetY())))
            {
                foundEnd = true;
                endPosition = minPosition;
                endDirection = minDirection;
            }

            if (!foundEnd || (maxPosition.GetY() < std::min(minPosition.GetY(), endPosition.GetY())))
            {
                foundEnd = true;
                endPosition = maxPosition;
                endDirection = maxDirection;
            }
        }
        catch(StatusCodeException &statusCodeException)
        {
            if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                throw statusCodeException;

            continue;
        }
    }

    if (!(foundVtx && foundEnd))
        return;

    this->SetParticleParameters(vtxPosition, vtxDirection, pPfo);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayVertexBuildingAlgorithm::BuildCosmicRayDaughter(ParticleFlowObject *const pDaughterPfo) const
{
    if (pDaughterPfo->GetParentPfoList().size() != 1)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const ParticleFlowObject *const pParentPfo = *(pDaughterPfo->GetParentPfoList().begin());

    ClusterList parentList, daughterList;
    LArPfoHelper::GetClusters(pParentPfo, TPC_3D, parentList);
    LArPfoHelper::GetClusters(pDaughterPfo, TPC_3D, daughterList);

    if (daughterList.empty() || parentList.empty())
        return;

    bool foundVtx(false);
    float vtxDistanceSquared(0.f);
    CartesianVector vtxPosition(0.f, 0.f, 0.f);

    for (ClusterList::const_iterator dIter = daughterList.begin(), dIterEnd = daughterList.end(); dIter != dIterEnd; ++dIter)
    {
        const Cluster *const pDaughterCluster = *dIter;

        for (ClusterList::const_iterator pIter = parentList.begin(), pIterEnd = parentList.end(); pIter != pIterEnd; ++pIter)
        {
            const Cluster *const pParentCluster = *pIter;

            CartesianVector closestDaughterPosition(0.f, 0.f, 0.f), closestParentPosition(0.f, 0.f, 0.f);
            LArClusterHelper::GetClosestPositions(pDaughterCluster, pParentCluster, closestDaughterPosition, closestParentPosition);

            const float closestDistanceSquared((closestDaughterPosition - closestParentPosition).GetMagnitudeSquared());

            if (!foundVtx || closestDistanceSquared < vtxDistanceSquared)
            {
                foundVtx = true;
                vtxDistanceSquared = closestDistanceSquared;
                vtxPosition = (m_useParentShowerVertex ? closestParentPosition : closestDaughterPosition);
            }
        }
    }

    if (!foundVtx)
        return;

    this->SetParticleParameters(vtxPosition, CartesianVector(0.f, 0.f, 0.f), pDaughterPfo);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayVertexBuildingAlgorithm::SetParticleParameters(const CartesianVector &vtxPosition, const CartesianVector &vtxDirection,
    ParticleFlowObject *const pPfo) const
{
    if (!pPfo->GetVertexList().empty())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    pPfo->SetMomentum(vtxDirection);

    const VertexList *pVertexList = NULL; std::string vertexListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, vertexListName));

    PandoraContentApi::Vertex::Parameters parameters;
    parameters.m_position = vtxPosition;
    parameters.m_vertexType = VERTEX_3D;

    Vertex *pVertex(NULL);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertex));

    if (!pVertexList->empty())
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_vertexListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo<Vertex>(*this, pPfo, pVertex));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayVertexBuildingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputPfoListName", m_parentPfoListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputVertexListName", m_vertexListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseParentForShowerVertex", m_useParentShowerVertex));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitHalfWindow", m_halfWindowLayers));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
