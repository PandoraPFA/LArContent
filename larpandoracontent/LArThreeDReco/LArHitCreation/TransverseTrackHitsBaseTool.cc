/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/TransverseTrackHitsBaseTool.cc
 *
 *  @brief  Implementation of the transverse track hits base tool.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArHitCreation/TransverseTrackHitsBaseTool.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

using namespace pandora;

namespace lar_content
{

void TransverseTrackHitsBaseTool::GetTrackHits3D(
    const CaloHitVector &inputTwoDHits, const MatchedSlidingFitMap &matchedSlidingFitMap, ProtoHitVector &protoHitVector) const
{
    for (const CaloHit *const pCaloHit2D : inputTwoDHits)
    {
        try
        {
            ProtoHit protoHit(pCaloHit2D);
            this->GetTransverseTrackHit3D(matchedSlidingFitMap, protoHit);
            this->AddTransverseChi2(matchedSlidingFitMap, protoHit);

            if (protoHit.IsPositionSet() && (protoHit.GetChi2() < m_chiSquaredCut))
                protoHitVector.emplace_back(protoHit);
        }
        catch (StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TransverseTrackHitsBaseTool::AddTransverseChi2(const MatchedSlidingFitMap &matchedSlidingFitMap, ProtoHit &protoHit) const
{
    // TODO Develop a proper treatment of the |dz/dx| * sigmaX uncertainty
    double chiSquared(protoHit.GetChi2());
    const HitType inputHitType(protoHit.GetParentCaloHit2D()->GetHitType());
    const CartesianVector &inputPosition3D(protoHit.GetPosition3D());

    for (const MatchedSlidingFitMap::value_type &mapEntry : matchedSlidingFitMap)
    {
        if (mapEntry.first == inputHitType)
            continue;

        const CartesianVector inputPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), inputPosition3D, mapEntry.first));
        chiSquared += this->GetTransverseChi2(inputPosition2D, mapEntry.second);
    }

    protoHit.SetPosition3D(inputPosition3D, chiSquared);
}

//------------------------------------------------------------------------------------------------------------------------------------------

double TransverseTrackHitsBaseTool::GetTransverseChi2(const CartesianVector &position2D, const TwoDSlidingFitResult &fitResult) const
{
    const float minLayerX(fitResult.GetGlobalMinLayerPosition().GetX());
    const float maxLayerX(fitResult.GetGlobalMaxLayerPosition().GetX());

    const float minX(std::min(minLayerX, maxLayerX));
    const float maxX(std::max(minLayerX, maxLayerX));

    if (((position2D.GetX() - minX) > -std::numeric_limits<float>::epsilon()) && ((position2D.GetX() - maxX) < +std::numeric_limits<float>::epsilon()))
        return 0.;

    const float minLayerDistanceSquared((position2D - fitResult.GetGlobalMinLayerPosition()).GetMagnitudeSquared());
    const float maxLayerDistanceSquared((position2D - fitResult.GetGlobalMaxLayerPosition()).GetMagnitudeSquared());

    double px(0.), pz(0.);

    if (minLayerDistanceSquared < maxLayerDistanceSquared)
    {
        px = fitResult.GetGlobalMinLayerDirection().GetX();
        pz = fitResult.GetGlobalMinLayerDirection().GetZ();
    }
    else if (maxLayerDistanceSquared < minLayerDistanceSquared)
    {
        px = fitResult.GetGlobalMaxLayerDirection().GetX();
        pz = fitResult.GetGlobalMaxLayerDirection().GetZ();
    }

    if (std::fabs(px) > std::numeric_limits<double>::epsilon())
    {
        return (0.5 * (pz * pz) / (px * px)); // TODO: Figure out the appropriate scale factor
    }

    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

} // namespace lar_content
