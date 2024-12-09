/**
 *  @file   larpandoradlcontent/LArThreeDReco/LArEventBuilding/MLPPrimaryHierarchyTool.cc
 *
 *  @brief  Implementation of the MLP primary hierarchy tool
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/StatusCodes.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/LArHierarchyPfo.h"
#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/MLPPrimaryHierarchyTool.h"

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

MLPPrimaryHierarchyTool::MLPPrimaryHierarchyTool() :
    m_pfoListNames({"TrackParticles3D", "ShowerParticles3D"}),
    m_nSpacepointsMin(0.f),
    m_nSpacepointsMax(2000.f),
    m_nuSeparationMin(-50.f),
    m_nuSeparationMax(500.f),
    m_vertexRegionRadius(5.f),
    m_vertexRegionNHitsMin(-10.f),
    m_vertexRegionNHitsMax(100.f),
    m_vertexRegionNParticlesMin(-1.f),
    m_vertexRegionNParticlesMax(8.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MLPPrimaryHierarchyTool::Run(const Algorithm *const pAlgorithm, const HierarchyPfo &hierarchyPfo, const ParticleFlowObject *const pNeutrinoPfo, 
    const bool useUpstream)
{
    // Pick out the correct pfo vertex/direction
    const CartesianVector particleVertex(useUpstream ? hierarchyPfo.GetUpstreamVertex() : hierarchyPfo.GetDownstreamVertex());
    const CartesianVector particleDirection(useUpstream ? hierarchyPfo.GetUpstreamDirection() : hierarchyPfo.GetDownstreamDirection());

    // Check that we could actually calculate pfo vertex/direction, return if not
    if ((particleVertex.GetZ() < -990.f) || (particleVertex.GetZ() < -990.f))
        return STATUS_CODE_NOT_INITIALIZED;

    if ((particleDirection.GetZ() < -990.f) || (particleDirection.GetZ() < -990.f))
        return STATUS_CODE_NOT_INITIALIZED;

    MLPPrimaryNetworkParams primaryNetworkParams;

    // Set primaryNetworkParams
    this->SetNSpacepoints(hierarchyPfo.GetPfo(), primaryNetworkParams);
    this->SetNuVertexSep(particleVertex, pNeutrinoPfo, primaryNetworkParams);
    this->SetVertexRegionParams(pAlgorithm, hierarchyPfo.GetPfo(), particleVertex, primaryNetworkParams);
    /////////////////////////////////////////
    std::cout << "BEFORE NORMALISATION" << std::endl;
    std::cout << "NSpacepoints: " << primaryNetworkParams.m_nSpacepoints << std::endl;
    std::cout << "NuVertexSep: " << primaryNetworkParams.m_nuSeparation << std::endl;
    std::cout << "StartRegionNHits:" << primaryNetworkParams.m_vertexRegionNHits << std::endl;
    std::cout << "StartRegionNParticles:" << primaryNetworkParams.m_vertexRegionNParticles << std::endl;
    /////////////////////////////////////////

    // Normalise
    this->NormaliseNetworkParams(primaryNetworkParams);

    /////////////////////////////////////////
    // std::cout << "AFTER NORMALISATION" << std::endl;
    // std::cout << "NSpacepoints: " << primaryNetworkParams.m_nSpacepoints << std::endl;
    // std::cout << "NuVertexSep: " << primaryNetworkParams.m_nuSeparation << std::endl;
    /////////////////////////////////////////

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPPrimaryHierarchyTool::SetNSpacepoints(const ParticleFlowObject *const pPfo, MLPPrimaryNetworkParams &primaryNetworkParams) const
{
    ClusterList clusterList3D;
    LArPfoHelper::GetThreeDClusterList(pPfo, clusterList3D);

    int total3DHits(0);

    for (const Cluster *const pCluster3D : clusterList3D)
        total3DHits += pCluster3D->GetNCaloHits();

    primaryNetworkParams.m_nSpacepoints = total3DHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPPrimaryHierarchyTool::SetNuVertexSep(const CartesianVector &particleVertex, const ParticleFlowObject *const pNeutrinoPfo,
    MLPPrimaryNetworkParams &primaryNetworkParams) const
{ 
    if (pNeutrinoPfo->GetVertexList().empty())
        return;
 
    // Get the neutrino vertex
    const Vertex *const pNeutrinoVertex(LArPfoHelper::GetVertex(pNeutrinoPfo));
    const CartesianVector nuVertex(pNeutrinoVertex->GetPosition().GetX(), pNeutrinoVertex->GetPosition().GetY(), pNeutrinoVertex->GetPosition().GetZ());

    // Work out separation
    const float separation((particleVertex - nuVertex).GetMagnitude());

    primaryNetworkParams.m_nuSeparation = separation;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPPrimaryHierarchyTool::SetVertexRegionParams(const Algorithm *const pAlgorithm, const ParticleFlowObject *const pPfo, 
    const CartesianVector &particleVertex, MLPPrimaryNetworkParams &primaryNetworkParams) const
{ 
    int hitCount = 0;
    int particleCount = 0;    

    for (const std::string &pfoListName : m_pfoListNames)
    {
        const PfoList *pPfoList(nullptr);
        if (PandoraContentApi::GetList(*pAlgorithm, pfoListName, pPfoList) != STATUS_CODE_SUCCESS)
            continue;

        for (const ParticleFlowObject *const pOtherPfo : *pPfoList)
        {
            if (pPfo == pOtherPfo)
                continue;

            std::cout << "YO" << std::endl;

            bool isClose(false);

            CartesianPointVector otherPfoPositions3D;
            LArPfoHelper::GetCoordinateVector(pOtherPfo, TPC_3D, otherPfoPositions3D);

            for (const CartesianVector &otherPfoPosition : otherPfoPositions3D)
            {
                const double sepSq = (otherPfoPosition - particleVertex).GetMagnitudeSquared();

                if (sepSq < (m_vertexRegionRadius * m_vertexRegionRadius))
                {
                    isClose = true;
                    ++hitCount;
                }
            }

            if (isClose)
                ++particleCount;
        }
    }

    primaryNetworkParams.m_vertexRegionNHits = hitCount;
    primaryNetworkParams.m_vertexRegionNParticles = particleCount;
}


//------------------------------------------------------------------------------------------------------------------------------------------

void MLPPrimaryHierarchyTool::NormaliseNetworkParams(MLPPrimaryNetworkParams &primaryNetworkParams) const
{
    this->NormaliseNetworkParam(m_nSpacepointsMin, m_nSpacepointsMax, primaryNetworkParams.m_nSpacepoints);
    this->NormaliseNetworkParam(m_nuSeparationMin, m_nuSeparationMax, primaryNetworkParams.m_nuSeparation);
    this->NormaliseNetworkParam(m_vertexRegionNHitsMin, m_vertexRegionNHitsMax, primaryNetworkParams.m_vertexRegionNHits);
    this->NormaliseNetworkParam(m_vertexRegionNParticlesMin, m_vertexRegionNParticlesMax, primaryNetworkParams.m_vertexRegionNParticles);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLPPrimaryHierarchyTool::NormaliseNetworkParam(const float minLimit, const float maxLimit, float &primaryNetworkParam) const
{
    const float interval(minLimit + maxLimit);

    if (primaryNetworkParam < minLimit)
        primaryNetworkParam = minLimit;

    if (primaryNetworkParam > maxLimit)
        primaryNetworkParam = maxLimit;

    primaryNetworkParam /= interval;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MLPPrimaryHierarchyTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "PfoListNames", m_pfoListNames));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NSpacepointsMin", m_nSpacepointsMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NSpacepointsMax", m_nSpacepointsMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NuSeparationMin", m_nuSeparationMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NuSeparationMax", m_nuSeparationMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexRegionRadius", m_vertexRegionRadius));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexRegionNHitsMin", m_vertexRegionNHitsMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexRegionNHitsMax", m_vertexRegionNHitsMax));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexRegionNParticlesMin", m_vertexRegionNParticlesMin));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexRegionNParticlesMax", m_vertexRegionNParticlesMax));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
