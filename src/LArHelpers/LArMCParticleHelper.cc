/**
 *  @file   LArContent/src/LArHelpers/LArMCParticleHelper.cc
 * 
 *  @brief  Implementation of the lar monte carlo particle helper class.
 * 
 *  $Log: $
 */

#include "Objects/MCParticle.h"
#include "Objects/CaloHit.h"
#include "Objects/Cluster.h"

#include "Pandora/PdgTable.h"

#include "LArHelpers/LArMCParticleHelper.h"

#include <cstdlib>

namespace lar_content
{

using namespace pandora;

bool LArMCParticleHelper::IsNeutrinoFinalState(const MCParticle *const pMCParticle)
{
    return ((pMCParticle->GetParentList().size() == 1) && (LArMCParticleHelper::IsNeutrino(*(pMCParticle->GetParentList().begin()))));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsNeutrinoInduced(const MCParticle *const pMCParticle)
{
    return (LArMCParticleHelper::GetParentNeutrinoId(pMCParticle) != 0);
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

const MCParticle *LArMCParticleHelper::GetParentNeutrino(const MCParticle *const pMCParticle)
{
    const MCParticle *const pParentMCParticle = LArMCParticleHelper::GetParentMCParticle(pMCParticle);  

    if(!LArMCParticleHelper::IsNeutrino(pParentMCParticle))
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return pParentMCParticle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

int LArMCParticleHelper::GetParentNeutrinoId(const MCParticle *const pMCParticle)
{
    try
    {
        const MCParticle *const pParentMCParticle = LArMCParticleHelper::GetParentNeutrino(pMCParticle);
        return pParentMCParticle->GetParticleId();
    }
    catch (const StatusCodeException &)
    {
        return 0;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsNeutrinoInduced(const Cluster *const pCluster, const float minWeight)
{
    return (LArMCParticleHelper::GetNeutrinoWeight(pCluster) > minWeight);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsNeutrinoInduced(const CaloHit *const pCaloHit, const float minWeight)
{
    return (LArMCParticleHelper::GetNeutrinoWeight(pCaloHit) > minWeight);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArMCParticleHelper::GetNeutrinoWeight(const Cluster *const pCluster)
{
    float neutrinoWeight(0.f);
    float totalWeight(0.f);

    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            const CaloHit *const pCaloHit = *hitIter;

            try
            {
                // note: order is important here
                neutrinoWeight += LArMCParticleHelper::GetNeutrinoWeight(pCaloHit);
                totalWeight += 1.f;
            }
            catch (const StatusCodeException &)
            {
            }
        }
    }

    if (totalWeight > std::numeric_limits<float>::epsilon())
        return (neutrinoWeight/totalWeight);

    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArMCParticleHelper::GetNeutrinoWeight(const CaloHit *const pCaloHit)
{
    const MCParticleWeightMap &weights = pCaloHit->GetMCParticleWeightMap();

    if (weights.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    float neutrinoWeight(0.f);
    float totalWeight(0.f);

    for (MCParticleWeightMap::const_iterator iter = weights.begin(), iterEnd = weights.end(); iter != iterEnd; ++iter)
    {
        const MCParticle *const pMCParticle = iter->first;
        const float weight = iter->second;

        if (LArMCParticleHelper::IsNeutrinoInduced(pMCParticle))
            neutrinoWeight += weight;

        totalWeight += weight;
    }

    if (totalWeight > std::numeric_limits<float>::epsilon())
        return (neutrinoWeight/totalWeight);

    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMCParticleHelper::IsPrimary(const pandora::MCParticle *const pMCParticle)
{
    // TODO This will require careful thought and modification as the use-case is clarified
    return (LArMCParticleHelper::IsNeutrinoFinalState(pMCParticle) ||
        (pMCParticle->GetParentList().empty() && (PI_ZERO != pMCParticle->GetParticleId())) ||
        ((1 == pMCParticle->GetParentList().size()) && (PI_ZERO == (*(pMCParticle->GetParentList().begin()))->GetParticleId())) );
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMCParticleHelper::GetMCPrimaryMap(const MCParticleList *const pMCParticleList, MCRelationMap &mcPrimaryMap)
{
    for (MCParticleList::const_iterator iter = pMCParticleList->begin(), iterEnd = pMCParticleList->end(); iter != iterEnd; ++iter)
    {
        const MCParticle *const pMCParticle = *iter;

        if (LArMCParticleHelper::IsPrimary(pMCParticle))
        {
            if (!LArMCParticleHelper::IsNeutrino(pMCParticle))
                mcPrimaryMap[pMCParticle] = pMCParticle;
        }
        else if (pMCParticle->GetParentList().empty())
        {
            // ATTN Protect against case where a parentless particle is denied primary status
            continue;
        }
        else
        {
            const MCParticle *pParentMCParticle = pMCParticle;

            while (true)
            {
                pParentMCParticle = *(pParentMCParticle->GetParentList().begin());

                if (mcPrimaryMap.find(pParentMCParticle) != mcPrimaryMap.end())
                {
                    mcPrimaryMap[pMCParticle] = mcPrimaryMap[pParentMCParticle];
                    break;
                }
                else if (LArMCParticleHelper::IsPrimary(pParentMCParticle))
                {
                    mcPrimaryMap[pMCParticle] = pParentMCParticle;
                    break;
                }
            }
        }
    }
}

} // namespace lar_content
