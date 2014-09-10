/**
 *  @file   LArContent/src/LArThreeDReco/LArHitCreation/DeltaRayShowerHitsTool.cc
 * 
 *  @brief  Implementation of the delta ray shower hits tool.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArGeometryHelper.h"
#include "LArHelpers/LArPfoHelper.h"

#include "LArThreeDReco/LArHitCreation/DeltaRayShowerHitsTool.h"
#include "LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

void DeltaRayShowerHitsTool::Run(ThreeDHitCreationAlgorithm *pAlgorithm, const ParticleFlowObject *const pPfo, const CaloHitList &inputTwoDHits,
    CaloHitList &newThreeDHits)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << this->GetType() << std::endl;

    try
    {
        if (!LArPfoHelper::IsShower(pPfo))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        if (pPfo->GetParentPfoList().size() != 1)
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);

        const ParticleFlowObject *pParentPfo = *(pPfo->GetParentPfoList().begin());

        CaloHitList caloHitList3D;
        LArPfoHelper::GetCaloHits(pParentPfo, TPC_3D, caloHitList3D);

        if (caloHitList3D.empty())
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);

        this->CreateThreeDHits(pAlgorithm, inputTwoDHits, caloHitList3D, newThreeDHits);
    }
    catch (StatusCodeException &)
    {
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayShowerHitsTool::CreateThreeDHits(ThreeDHitCreationAlgorithm *pAlgorithm, const CaloHitList &inputTwoDHits, const CaloHitList &caloHitList3D, 
    CaloHitList &newThreeDHits) const
{
    for (CaloHitList::const_iterator iter1 = inputTwoDHits.begin(), iterEnd1 = inputTwoDHits.end(); iter1 != iterEnd1; ++iter1)
    {
        try
        {
            CaloHit *pCaloHit2D(*iter1);

            const HitType hitType(pCaloHit2D->GetHitType());
            const HitType hitType1((TPC_VIEW_U == hitType) ? TPC_VIEW_V : (TPC_VIEW_V == hitType) ? TPC_VIEW_W : TPC_VIEW_U);
            const HitType hitType2((TPC_VIEW_U == hitType) ? TPC_VIEW_W : (TPC_VIEW_V == hitType) ? TPC_VIEW_U : TPC_VIEW_V);

            bool foundClosestPosition(false);
            float closestDistanceSquared(std::numeric_limits<float>::max());
            CartesianVector closestPosition3D(0.f, 0.f, 0.f);

            for (CaloHitList::const_iterator iter2 = caloHitList3D.begin(), iterEnd2 = caloHitList3D.end(); iter2 != iterEnd2; ++iter2)
            {
                const CartesianVector thisPosition3D((*iter2)->GetPositionVector());
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

            CaloHit *pCaloHit3D(NULL);
            pAlgorithm->CreateThreeDHit(pCaloHit2D, position3D, pCaloHit3D);
            newThreeDHits.insert(pCaloHit3D);
        }
        catch (StatusCodeException &)
        {
        }
    }
}

} // namespace lar_content
