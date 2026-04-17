/**
 *  @file   larpandoracontent/LArMetrics/BaseValidationTool.cc
 *
 *  @brief  Implementation of the base validation tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMetrics/BaseValidationTool.h"

using namespace pandora;

namespace lar_content
{

void BaseValidationTool::GetHitsOfType(const CaloHitList &inputList, const HitType hitType, CaloHitVector &outputVector, float &totalEnergy)
{
    totalEnergy = 0.f;

    for (const CaloHit *const pCaloHit : inputList)
    {
        if (pCaloHit->GetHitType() == hitType)
        {
            totalEnergy += pCaloHit->GetElectromagneticEnergy();
            outputVector.push_back(pCaloHit);
        }
    }
}

} // namespace lar_content
