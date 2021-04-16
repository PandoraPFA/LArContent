/**
 *  @file   larpandoracontent/LArCheating/CheatingClusterCharacterisationAlgorithm.cc
 *
 *  @brief  Implementation of the cheating cluster characterisation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArCheating/CheatingClusterCharacterisationAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

using namespace pandora;

namespace lar_content
{

bool CheatingClusterCharacterisationAlgorithm::IsClearTrack(const Cluster *const pCluster) const
{
    try
    {
        // ATTN Slightly curious definition of a clear track, but this is most-likely what is needed for shower-growing
        const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCluster));

        if ((PHOTON != pMCParticle->GetParticleId()) && (E_MINUS != std::abs(pMCParticle->GetParticleId())))
            return true;
    }
    catch (const StatusCodeException &)
    {
    }

    return false;
}

} // namespace lar_content
