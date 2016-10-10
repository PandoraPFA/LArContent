/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/DeltaRayShowerHitsTool.cc
 * 
 *  @brief  Implementation of the delta ray shower hits tool.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArThreeDReco/LArHitCreation/DeltaRayShowerHitsTool.h"
#include "larpandoracontent/LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

void DeltaRayShowerHitsTool::Run(ThreeDHitCreationAlgorithm *const pAlgorithm, const ParticleFlowObject *const pPfo, const CaloHitVector &inputTwoDHits,
    CaloHitVector &newThreeDHits)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    try
    {
        if (!LArPfoHelper::IsShower(pPfo))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        if (pPfo->GetParentPfoList().size() != 1)
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);

        const ParticleFlowObject *const pParentPfo = *(pPfo->GetParentPfoList().begin());

        CaloHitList caloHitList3D;
        LArPfoHelper::GetCaloHits(pParentPfo, TPC_3D, caloHitList3D);

        if (caloHitList3D.empty())
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);

        CaloHitVector caloHitVector3D(caloHitList3D.begin(), caloHitList3D.end());
        std::sort(caloHitVector3D.begin(), caloHitVector3D.end(), LArClusterHelper::SortHitsByPosition);

        this->CreateThreeDHits(pAlgorithm, inputTwoDHits, caloHitVector3D, newThreeDHits);
    }
    catch (StatusCodeException &)
    {
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayShowerHitsTool::CreateThreeDHits(ThreeDHitCreationAlgorithm *const pAlgorithm, const CaloHitVector &inputTwoDHits, const CaloHitVector &caloHitList3D, 
    CaloHitVector &newThreeDHits) const
{
    for (const CaloHit *const pCaloHit2D : inputTwoDHits)
    {
        try
        {
            const HitType hitType(pCaloHit2D->GetHitType());
            const HitType hitType1((TPC_VIEW_U == hitType) ? TPC_VIEW_V : (TPC_VIEW_V == hitType) ? TPC_VIEW_W : TPC_VIEW_U);
            const HitType hitType2((TPC_VIEW_U == hitType) ? TPC_VIEW_W : (TPC_VIEW_V == hitType) ? TPC_VIEW_U : TPC_VIEW_V);

            bool foundClosestPosition(false);
            float closestDistanceSquared(std::numeric_limits<float>::max());
            CartesianVector closestPosition3D(0.f, 0.f, 0.f);

            for (const CaloHit *const pCaloHit3D : caloHitList3D)
            {
                const CartesianVector thisPosition3D(pCaloHit3D->GetPositionVector());
                const CartesianVector thisPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), thisPosition3D, hitType));
                const float thisDistanceSquared((pCaloHit2D->GetPositionVector() - thisPosition2D).GetMagnitudeSquared());

                if (thisDistanceSquared <  closestDistanceSquared)
                {
                    foundClosestPosition = true;
                    closestDistanceSquared = thisDistanceSquared;
                    closestPosition3D = thisPosition3D;
                }
            }

            if (!foundClosestPosition)
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);

            const CartesianVector position1(LArGeometryHelper::ProjectPosition(this->GetPandora(), closestPosition3D, hitType1));
            const CartesianVector position2(LArGeometryHelper::ProjectPosition(this->GetPandora(), closestPosition3D, hitType2));

            CartesianVector position3D(0.f, 0.f, 0.f);
            float chiSquared(std::numeric_limits<float>::max());
            this->GetPosition3D(pCaloHit2D, hitType1, hitType2, position1, position2, position3D, chiSquared);

            const CaloHit *pCaloHit3D(NULL);
            pAlgorithm->CreateThreeDHit(pCaloHit2D, position3D, pCaloHit3D);
            newThreeDHits.push_back(pCaloHit3D);
        }
        catch (StatusCodeException &)
        {
        }
    }
}

} // namespace lar_content
