/**
 *  @file   LArContent/src/LArThreeDReco/LArEventBuilding/CosmicRayBuildingAlgorithm.cc
 *
 *  @brief  Implementation of the cosmic-ray building algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArGeometryHelper.h"
#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArPfoHelper.h"
#include "LArPlugins/LArTransformationPlugin.h"

#include "LArThreeDReco/LArEventBuilding/CosmicRayBuildingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CosmicRayBuildingAlgorithm::CosmicRayBuildingAlgorithm() :
    m_halfWindowLayers(30)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayBuildingAlgorithm::Run()
{
    try
    {
        PfoList pfoList;
        this->GetInputPfoList(pfoList);

        ThreeDSlidingFitResultMap slidingFitResultMap;
        this->BuildSlidingFitResultMap(pfoList, slidingFitResultMap);
        this->BuildCosmicRayParticles(slidingFitResultMap, pfoList);
    }
    catch (StatusCodeException &)
    {
        std::cout << "CosmicRayBuildingAlgorithm: unable to build cosmic rays." << std::endl;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayBuildingAlgorithm::GetInputPfoList(PfoList &pfoList) const
{
    for (StringVector::const_iterator iter = m_pfoListNames.begin(), iterEnd = m_pfoListNames.end(); iter != iterEnd; ++iter)
    {
        const PfoList *pPfoList(NULL);

        if (STATUS_CODE_SUCCESS == PandoraContentApi::GetList(*this, *iter, pPfoList))
        {
            pfoList.insert(pPfoList->begin(), pPfoList->end());
        }
        else
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "CosmicRayBuildingAlgorithm: pfo list " << *iter << " unavailable." << std::endl;
        }
    }

    if (pfoList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayBuildingAlgorithm::BuildSlidingFitResultMap(const PfoList &pfoList, ThreeDSlidingFitResultMap &slidingFitResultMap) const
{
    const float slidingFitPitch(LArGeometryHelper::GetLArTransformationPlugin(this->GetPandora())->GetWireZPitch());

    for (PfoList::const_iterator pIter = pfoList.begin(), pIterEnd = pfoList.end(); pIter != pIterEnd; ++pIter)
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
                const ThreeDSlidingFitResult slidingFitResult(pCluster, m_halfWindowLayers, slidingFitPitch);

                if (!slidingFitResultMap.insert(ThreeDSlidingFitResultMap::value_type(pCluster, slidingFitResult)).second)
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

void CosmicRayBuildingAlgorithm::BuildCosmicRayParticles(const ThreeDSlidingFitResultMap &slidingFitResultMap, const PfoList &pfoList) const
{
    for (PfoList::const_iterator iter = pfoList.begin(), iterEnd = pfoList.end(); iter != iterEnd; ++iter)
    {
        ParticleFlowObject *const pPfo = *iter;

        if (LArPfoHelper::IsFinalState(pPfo))
        {
            this->BuildCosmicRayParent(slidingFitResultMap, pPfo);
        }
        else
        {
            this->BuildCosmicRayDaughter(pPfo);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayBuildingAlgorithm::BuildCosmicRayParent(const ThreeDSlidingFitResultMap &slidingFitResultMap, ParticleFlowObject *const pPfo) const
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

            ThreeDSlidingFitResultMap::const_iterator cIter2 = slidingFitResultMap.find(pCluster);

            if (slidingFitResultMap.end() != cIter2)
            {
                const ThreeDSlidingFitResult &slidingFitResult(cIter2->second);

                minPosition = slidingFitResult.GetGlobalMinLayerPosition();
                maxPosition = slidingFitResult.GetGlobalMaxLayerPosition();
                minDirection = slidingFitResult.GetGlobalMinLayerDirection();
                maxDirection = slidingFitResult.GetGlobalMaxLayerDirection();
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

    const float scaleFactor((vtxDirection.GetDotProduct(endPosition - vtxPosition) > 0.f) ? +1.f : -1.f);
    this->SetParticleParameters(vtxPosition, vtxDirection * scaleFactor, pPfo);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayBuildingAlgorithm::BuildCosmicRayDaughter(ParticleFlowObject *const pDaughterPfo) const
{
    if (pDaughterPfo->GetParentPfoList().size() != 1)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const ParticleFlowObject *const pParentPfo = *(pDaughterPfo->GetParentPfoList().begin());

    ClusterList parentList, daughterList;
    LArPfoHelper::GetClusters(pParentPfo, TPC_3D, parentList);
    LArPfoHelper::GetClusters(pDaughterPfo, TPC_3D, daughterList);

    // Take closest extremal position to parent Pfo as vertex of daughter Pfos (TODO: Make this more sophisticated)
    bool foundVtx(false);
    float vtxDistanceSquared(0.f);
    CartesianVector vtxPosition(0.f, 0.f, 0.f);
    CartesianVector vtxDirection(0.f, 0.f, 0.f);

    for (ClusterList::const_iterator dIter = daughterList.begin(), dIterEnd = daughterList.end(); dIter != dIterEnd; ++dIter)
    {
        const Cluster *const pDaughterCluster = *dIter;

        try
        {
            CartesianVector minDaughterPosition(0.f, 0.f, 0.f), maxDaughterPosition(0.f, 0.f, 0.f);
            LArClusterHelper::GetExtremalCoordinates(pDaughterCluster, minDaughterPosition, maxDaughterPosition);

            for (ClusterList::const_iterator pIter = parentList.begin(), pIterEnd = parentList.end(); pIter != pIterEnd; ++pIter)
            {
                const Cluster *const pParentCluster = *pIter;

                const CartesianVector minParentPosition(LArClusterHelper::GetClosestPosition(minDaughterPosition, pParentCluster));
                const CartesianVector maxParentPosition(LArClusterHelper::GetClosestPosition(maxDaughterPosition, pParentCluster));
                const CartesianVector minParentDirection((maxDaughterPosition - minParentPosition).GetUnitVector());
                const CartesianVector maxParentDirection((minDaughterPosition - maxParentPosition).GetUnitVector());

                const float minDistanceSquared((minParentPosition - minDaughterPosition).GetMagnitudeSquared());
                const float maxDistanceSquared((maxParentPosition - maxDaughterPosition).GetMagnitudeSquared());

                if (!foundVtx || (minDistanceSquared < vtxDistanceSquared))
                {
                    foundVtx = true;
                    vtxDistanceSquared = minDistanceSquared;
                    vtxPosition = minParentPosition;
                    vtxDirection = minParentDirection;
                }

                if (!foundVtx || (maxDistanceSquared < vtxDistanceSquared))
                {
                    foundVtx = true;
                    vtxDistanceSquared = maxDistanceSquared;
                    vtxPosition = maxParentPosition;
                    vtxDirection = maxParentDirection;
                }
            }
        }
        catch(StatusCodeException &statusCodeException)
        {
            if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                throw statusCodeException;

            continue;
        }
    }

    if (!foundVtx)
        return;

    this->SetParticleParameters(vtxPosition, vtxDirection, pDaughterPfo);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayBuildingAlgorithm::SetParticleParameters(const CartesianVector &vtxPosition, const CartesianVector &vtxDirection,
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

StatusCode CosmicRayBuildingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputPfoListNames", m_pfoListNames));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputVertexListName", m_vertexListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitHalfWindow", m_halfWindowLayers));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
