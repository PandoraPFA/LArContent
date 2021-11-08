/**
 *  @file   larpandoracontent/LArShowerRefinement/GammaStartRefinementTool.cc
 *
 *  @brief  Implementation of the gamma start refinement tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArShowerRefinement/GammaStartRefinementTool.h"

using namespace pandora;

namespace lar_content
{

GammaStartRefinementTool::GammaStartRefinementTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool GammaStartRefinementTool::Run(ShowerStartRefinementAlgorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo, const CartesianVector &nuVertexPosition)
{
    std::cout << "AAAAAAAAA" << std::endl;

    // only apply gamma refinement algorithm to showers
    if (!LArPfoHelper::IsShower(pShowerPfo))
        return false;

    std::cout << "BBBBBBBBBB" << std::endl;

    // ISOBEL - SHOULD I ALSO HAVE A HIT CUT?
    if (!this->HasPathToNuVertex(pShowerPfo, nuVertexPosition))
        return false;

    std::cout << "CCCCCCCCCCCC" << std::endl;

    ProtoShowerVector protoShowerVector;
    this->BuildProtoShowers(pShowerPfo, protoShowerVector);

    std::cout << "DDDDDDDDDDD" << std::endl;

    if (protoShowerVector.empty())
        return false;

    std::cout << "EEEEEEEEEEEE" << std::endl;

    // pfo splitting alg? (can be quite simple if we are able to build the protoShowers sufficiently...)
    // set ProtoShower parent pfo address to match!
    if (protoShowerVector.size() != 1)
        return false;

    for (const ProtoShower &protoShower : protoShowerVector)
    {
        if (this->IsElectronPathway(protoShower))
            continue;

        this->RemoveConnectionPathway(protoShower);

        // change metadata to say that we think that this gamma is a gamma (so that we don't extend it in the future tool)
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::RemoveConnectionPathway(const ProtoShower &protoShower)
{
    std::cout << protoShower.m_showerCore.m_startPosition.GetX() << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode GammaStartRefinementTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
