/**
 *  @file   larpandoracontent/LArThreeDReco/LArEventBuilding/NeutrinoDaughterVerticesAlgorithm.cc
 *
 *  @brief  Implementation of the neutrino daughter vertices algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArThreeDReco/LArEventBuilding/NeutrinoDaughterVerticesAlgorithm.h"

using namespace pandora;

namespace lar_content
{

NeutrinoDaughterVerticesAlgorithm::NeutrinoDaughterVerticesAlgorithm() : m_useParentShowerVertex(false), m_halfWindowLayers(20), m_cheatedPfoListName("")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode NeutrinoDaughterVerticesAlgorithm::Run()
{
    const PfoList *pPfoList = NULL;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_neutrinoListName, pPfoList));

    if (!pPfoList || pPfoList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "NeutrinoDaughterVerticesAlgorithm: unable to find pfo list " << m_neutrinoListName << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    const VertexList *pGammaVertexList = NULL;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, "GammaVertices", pGammaVertexList));

    PfoVector pfoVector;
    LArPointingClusterMap pointingClusterMap;

    this->GetDaughterPfos(pPfoList, pfoVector);
    this->BuildPointingClusterMap(pfoVector, pointingClusterMap);
    this->BuildDaughterParticles(pointingClusterMap, pfoVector, pGammaVertexList);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoDaughterVerticesAlgorithm::GetDaughterPfos(const PfoList *const pPfoList, PfoVector &pfoVector) const
{
    PfoList outputList;

    for (PfoList::const_iterator pIter = pPfoList->begin(), pIterEnd = pPfoList->end(); pIter != pIterEnd; ++pIter)
    {
        if (!LArPfoHelper::IsNeutrino(*pIter) && (*pIter)->GetVertexList().size() != 1)
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        LArPfoHelper::GetAllDownstreamPfos(*pIter, outputList);
    }

    const PfoList *pCheatedPfoList(nullptr);

    if (!m_cheatedPfoListName.empty())
      PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_cheatedPfoListName, pCheatedPfoList));

    for (PfoList::const_iterator pIter = outputList.begin(), pIterEnd = outputList.end(); pIter != pIterEnd; ++pIter)
    {
        ClusterList clusterList;
        LArPfoHelper::GetClusters(*pIter, TPC_3D, clusterList);

        if (clusterList.empty())
            continue;

	if (pCheatedPfoList && (std::find(pCheatedPfoList->begin(), pCheatedPfoList->end(), *pIter) != pCheatedPfoList->end()))
	  continue;

        pfoVector.push_back(*pIter);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoDaughterVerticesAlgorithm::BuildPointingClusterMap(const PfoVector &pfoList, LArPointingClusterMap &pointingClusterMap) const
{
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));

    for (PfoVector::const_iterator pIter = pfoList.begin(), pIterEnd = pfoList.end(); pIter != pIterEnd; ++pIter)
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

void NeutrinoDaughterVerticesAlgorithm::BuildDaughterParticles(const LArPointingClusterMap &pointingClusterMap, const PfoVector &pfoVector, 
    const VertexList *const pGammaVertexList) const
{
    for (PfoVector::const_iterator iter = pfoVector.begin(), iterEnd = pfoVector.end(); iter != iterEnd; ++iter)
    {
        const ParticleFlowObject *const pPfo = *iter;

        if (LArPfoHelper::IsTrack(pPfo))
        {
            this->BuildDaughterTrack(pointingClusterMap, pPfo, pGammaVertexList);
        }
        else
        {
            this->BuildDaughterShower(pPfo, pGammaVertexList);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoDaughterVerticesAlgorithm::BuildDaughterTrack(const LArPointingClusterMap &pointingClusterMap, const ParticleFlowObject *const pDaughterPfo, 
    const VertexList *const pGammaVertexList) const
{
    if (pDaughterPfo->GetParentPfoList().size() != 1)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const ParticleFlowObject *const pParentPfo = *(pDaughterPfo->GetParentPfoList().begin());

    ClusterList parentList, daughterList;
    LArPfoHelper::GetClusters(pParentPfo, TPC_3D, parentList);
    LArPfoHelper::GetClusters(pDaughterPfo, TPC_3D, daughterList);

    if (parentList.empty() && pParentPfo->GetVertexList().empty())
        return;

    bool foundVtx(false);
    float vtxDistance(0.f);
    CartesianVector vtxPosition(0.f, 0.f, 0.f);
    CartesianVector vtxDirection(0.f, 0.f, 0.f);

    for (ClusterList::const_iterator dIter = daughterList.begin(), dIterEnd = daughterList.end(); dIter != dIterEnd; ++dIter)
    {
        const Cluster *const pDaughterCluster = *dIter;

        CartesianVector minPosition(0.f, 0.f, 0.f), maxPosition(0.f, 0.f, 0.f);
        CartesianVector minDirection(0.f, 0.f, 0.f), maxDirection(0.f, 0.f, 0.f);
        bool foundDirection(false);

        LArPointingClusterMap::const_iterator cIter = pointingClusterMap.find(pDaughterCluster);

        if (pointingClusterMap.end() != cIter)
        {
            const LArPointingCluster &pointingCluster(cIter->second);

            minPosition = pointingCluster.GetInnerVertex().GetPosition();
            maxPosition = pointingCluster.GetOuterVertex().GetPosition();
            minDirection = pointingCluster.GetInnerVertex().GetDirection();
            maxDirection = pointingCluster.GetOuterVertex().GetDirection();
            foundDirection = true;
        }
        else
        {
            LArClusterHelper::GetExtremalCoordinates(pDaughterCluster, minPosition, maxPosition);
        }

        if ((maxPosition - minPosition).GetMagnitudeSquared() < std::numeric_limits<float>::epsilon())
            continue;

        if (!foundDirection)
        {
            minDirection = (maxPosition - minPosition).GetUnitVector();
            maxDirection = (minPosition - maxPosition).GetUnitVector();
        }

        float minDistance(std::numeric_limits<float>::max());
        float maxDistance(std::numeric_limits<float>::max());

        for (ClusterList::const_iterator pIter = parentList.begin(), pIterEnd = parentList.end(); pIter != pIterEnd; ++pIter)
        {
            const Cluster *const pParentCluster = *pIter;
            minDistance = std::min(minDistance, (LArClusterHelper::GetClosestDistance(minPosition, pParentCluster)));
            maxDistance = std::min(maxDistance, (LArClusterHelper::GetClosestDistance(maxPosition, pParentCluster)));
        }

        if (LArPfoHelper::IsNeutrino(pParentPfo) && !pParentPfo->GetVertexList().empty())
        {
            const Vertex *const pVertex = *(pParentPfo->GetVertexList().begin());
            minDistance = std::min(minDistance, (pVertex->GetPosition() - minPosition).GetMagnitude());
            maxDistance = std::min(maxDistance, (pVertex->GetPosition() - maxPosition).GetMagnitude());
        }

        if (!foundVtx || (minDistance < vtxDistance))
        {
            foundVtx = true;
            vtxDistance = minDistance;
            vtxPosition = minPosition;
            vtxDirection = minDirection;
        }

        if (!foundVtx || (maxDistance < vtxDistance))
        {
            foundVtx = true;
            vtxDistance = maxDistance;
            vtxPosition = maxPosition;
            vtxDirection = maxDirection;
        }
    }

    if (!foundVtx)
        return;

    this->SetParticleParameters(vtxPosition, vtxDirection, pDaughterPfo, pGammaVertexList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoDaughterVerticesAlgorithm::BuildDaughterShower(const ParticleFlowObject *const pDaughterPfo, const VertexList *const pGammaVertexList) const
{
    if (pDaughterPfo->GetParentPfoList().size() != 1)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const ParticleFlowObject *const pParentPfo = *(pDaughterPfo->GetParentPfoList().begin());

    if (LArPfoHelper::IsNeutrino(pParentPfo) && pParentPfo->GetVertexList().size() != 1)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    ClusterList parentList, daughterList;
    LArPfoHelper::GetClusters(pParentPfo, TPC_3D, parentList);
    LArPfoHelper::GetClusters(pDaughterPfo, TPC_3D, daughterList);

    if (daughterList.empty())
        return;

    if (LArPfoHelper::IsNeutrino(pParentPfo))
    {
        const Vertex *const pVertex = *(pParentPfo->GetVertexList().begin());
        const CartesianVector vtxPosition(
            m_useParentShowerVertex ? pVertex->GetPosition() : LArClusterHelper::GetClosestPosition(pVertex->GetPosition(), daughterList));

        return this->SetParticleParameters(vtxPosition, CartesianVector(0.f, 0.f, 0.f), pDaughterPfo, pGammaVertexList);
    }

    if (parentList.empty())
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

    this->SetParticleParameters(vtxPosition, CartesianVector(0.f, 0.f, 0.f), pDaughterPfo, pGammaVertexList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoDaughterVerticesAlgorithm::SetParticleParameters(const CartesianVector &vtxPosition, const CartesianVector &vtxDirection, 
    const ParticleFlowObject *const pPfo, const VertexList *const pGammaVertexList) const
{
    PandoraContentApi::ParticleFlowObject::Metadata metadata;
    metadata.m_momentum = vtxDirection;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, metadata));

    if (!pPfo->GetVertexList().empty())
    {
        if (std::find(pGammaVertexList->begin(), pGammaVertexList->end(), pPfo->GetVertexList().front()) != pGammaVertexList->end())
            return;

        throw StatusCodeException(STATUS_CODE_FAILURE);
    }

    const VertexList *pVertexList = NULL;
    std::string vertexListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, vertexListName));

    PandoraContentApi::Vertex::Parameters parameters;
    parameters.m_position = vtxPosition;
    parameters.m_vertexLabel = VERTEX_INTERACTION;
    parameters.m_vertexType = VERTEX_3D;

    const Vertex *pVertex(NULL);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertex));

    if (!pVertexList->empty())
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_vertexListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo<Vertex>(*this, pPfo, pVertex));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode NeutrinoDaughterVerticesAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NeutrinoPfoListName", m_neutrinoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CheatedPfoListName", m_cheatedPfoListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputVertexListName", m_vertexListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "UseParentForShowerVertex", m_useParentShowerVertex));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingFitHalfWindow", m_halfWindowLayers));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
