/**
 *  @file   LArContent/src/LArHelpers/LArMCParticleHelper.cc
 * 
 *  @brief  Implementation of the lar monte carlo particle helper class.
 * 
 *  $Log: $
 */

#include "Objects/MCParticle.h"

#include "LArHelpers/LArMCParticleHelper.h"

namespace lar
{

using namespace pandora;

bool LArMCParticleHelper::IsNeutrinoInduced(const MCParticle *const pMCParticle)
{
    return ((pMCParticle->GetParentList().size() == 1) && (LArMCParticleHelper::IsNeutrino(*(pMCParticle->GetParentList().begin()))));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsNeutrino(const MCParticle *const pMCParticle)
{
    const int absoluteParticleId(std::abs(pMCParticle->GetParticleId()));

    if ((NU_E == absoluteParticleId) || (NU_MU == absoluteParticleId) || (NU_TAU == absoluteParticleId))
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const MCParticle *LArMCParticleHelper::GetParentMCParticle(const MCParticle *const pMCParticle)
{
    const MCParticle *pParentMCParticle = pMCParticle;

    while (pParentMCParticle->GetParentList().empty() == false)
    {
        pParentMCParticle = *(pParentMCParticle->GetParentList().begin());
    }

    return pParentMCParticle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

int LArMCParticleHelper::GetPrimaryNeutrino(const MCParticle *const pMCParticle)
{
    const MCParticle *pParentMCParticle = LArMCParticleHelper::GetParentMCParticle(pMCParticle);

    if (LArMCParticleHelper::IsNeutrino(pParentMCParticle))
        return pParentMCParticle->GetParticleId();

    return 0;
}

} // namespace lar
