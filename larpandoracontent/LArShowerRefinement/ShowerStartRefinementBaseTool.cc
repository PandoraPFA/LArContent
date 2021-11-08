/**
 *  @file   larpandoracontent/LArShowerRefinement/ShowerStartRefinementBaseTool.cc
 *
 *  @brief  Implementation of the shower start refinement base algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArShowerRefinement/ShowerStartRefinementBaseTool.h"

using namespace pandora;

namespace lar_content
{

ShowerStartRefinementBaseTool::ShowerStartRefinementBaseTool() : m_maxDistanceForConnection(5.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ShowerStartRefinementBaseTool::HasPathToNuVertex(const ParticleFlowObject *const pShowerPfo, const CartesianVector &neutrinoVertex) const
{
    // this algorithm has to run after we have 3D hits but not after shower vertex creation
    // therefore, find the closest space point to the neutrino vertex

    ClusterList clusterList3D;
    LArPfoHelper::GetClusters(pShowerPfo, TPC_3D, clusterList3D);

    for (const Cluster *const pCluster3D : clusterList3D)
    {
        CaloHitList caloHitList3D;
        pCluster3D->GetOrderedCaloHitList().FillCaloHitList(caloHitList3D);

        for (const CaloHit *const caloHit3D : caloHitList3D)
        {
            const CartesianVector &hitPosition(caloHit3D->GetPositionVector());
            const float separationSquared((neutrinoVertex - hitPosition).GetMagnitudeSquared());

            if (separationSquared < (m_maxDistanceForConnection * m_maxDistanceForConnection))
                return true;
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerStartRefinementBaseTool::BuildProtoShowers(const ParticleFlowObject *const pShowerPfo, ProtoShowerVector &protoShowerVector) const
{
    // cheat this to find all shower starts of significant MC particles in the pfo
    this->FindShowerCores(pShowerPfo, protoShowerVector);

    if (protoShowerVector.empty())
    {
        std::cout << "NO SHOWER CORES FOUND" << std::endl;
        return;
    }

    // bail if it can't be done? or take out of vector?
    this->FindShowerStartPositions(pShowerPfo, protoShowerVector);

    this->FindConnectionPathways(pShowerPfo, protoShowerVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerStartRefinementBaseTool::FindShowerCores(const ParticleFlowObject *const pShowerPfo, ProtoShowerVector &protoShowerVector) const
{
    // do fancy logic
    std::cout << pShowerPfo << protoShowerVector.size() << std::endl;
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerStartRefinementBaseTool::FindShowerStartPositions(const ParticleFlowObject *const pShowerPfo, ProtoShowerVector &protoShowerVector) const
{
    // do fancy logic
    std::cout << pShowerPfo << protoShowerVector.size() << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerStartRefinementBaseTool::FindConnectionPathways(const ParticleFlowObject *const pShowerPfo, ProtoShowerVector &protoShowerVector) const
{
    // do fancy logic
    std::cout << pShowerPfo << protoShowerVector.size() << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ShowerStartRefinementBaseTool::IsElectronPathway(const ProtoShower &protoShower)
{
    std::cout << protoShower.m_showerCore.m_startPosition.GetX() << std::endl;
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerStartRefinementBaseTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxDistanceForConnection", m_maxDistanceForConnection));

    return STATUS_CODE_SUCCESS;
}


} // namespace lar_content
