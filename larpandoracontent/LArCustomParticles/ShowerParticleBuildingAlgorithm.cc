/**
 *  @file   larpandoracontent/LArCustomParticles/ShowerParticleBuildingAlgorithm.cc
 *
 *  @brief  Implementation of the 3D shower building algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArObjects/LArShowerPfo.h"

#include "larpandoracontent/LArCustomParticles/ShowerParticleBuildingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ShowerParticleBuildingAlgorithm::ShowerParticleBuildingAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerParticleBuildingAlgorithm::CreatePfo(const ParticleFlowObject *const /*pInputPfo*/, const ParticleFlowObject *&/*pOutputPfo*/) const
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerParticleBuildingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    return CustomParticleCreationAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
