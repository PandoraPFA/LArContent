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

void DeltaRayShowerHitsTool::Run(ThreeDHitCreationAlgorithm *const pAlgorithm, const ParticleFlowObject *const pPfo,
    const CaloHitVector &inputTwoDHits, ProtoHitVector &protoHitVector)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    try
    {
        if (!LArPfoHelper::IsShower(pPfo) || (1 != pPfo->GetParentPfoList().size()))
            return;

        const ParticleFlowObject *const pParentPfo(pPfo->GetParentPfoList().front());

        CaloHitList parentHitList3D;
        LArPfoHelper::GetCaloHits(pParentPfo, TPC_3D, parentHitList3D);

        if (parentHitList3D.empty())
            return;

        CaloHitVector parentHitVector3D(parentHitList3D.begin(), parentHitList3D.end());
        std::sort(parentHitVector3D.begin(), parentHitVector3D.end(), LArClusterHelper::SortHitsByPosition);

        this->CreateDeltaRayShowerHits3D(inputTwoDHits, parentHitVector3D, protoHitVector);
    }
    catch (StatusCodeException &)
    {
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayShowerHitsTool::CreateDeltaRayShowerHits3D(
    const CaloHitVector &inputTwoDHits, const CaloHitVector &parentHits3D, ProtoHitVector &protoHitVector) const
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

            for (const CaloHit *const pCaloHit3D : parentHits3D)
            {
                const CartesianVector thisPosition3D(pCaloHit3D->GetPositionVector());
                const CartesianVector thisPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), thisPosition3D, hitType));
                const float thisDistanceSquared((pCaloHit2D->GetPositionVector() - thisPosition2D).GetMagnitudeSquared());

                if (thisDistanceSquared < closestDistanceSquared)
                {
                    foundClosestPosition = true;
                    closestDistanceSquared = thisDistanceSquared;
                    closestPosition3D = thisPosition3D;
                }
            }

            if (!foundClosestPosition)
                continue;

            const CartesianVector position1(LArGeometryHelper::ProjectPosition(this->GetPandora(), closestPosition3D, hitType1));
            const CartesianVector position2(LArGeometryHelper::ProjectPosition(this->GetPandora(), closestPosition3D, hitType2));

            ProtoHit protoHit(pCaloHit2D);
            this->GetBestPosition3D(hitType1, hitType2, position1, position2, protoHit);

            if (protoHit.IsPositionSet())
                protoHitVector.emplace_back(protoHit);
        }
        catch (StatusCodeException &)
        {
        }
    }
}

} // namespace lar_content
