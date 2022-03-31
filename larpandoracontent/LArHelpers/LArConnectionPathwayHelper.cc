/**
 *  @file   larpandoracontent/LArHelpers/LArConnectionPathwayHelper.cc
 *
 *  @brief  Implementation of the connection pathway helper class.
 *
 *  $Log: $
 */

#include "Objects/CartesianVector.h"
#include "Objects/CaloHit.h"

#include "Pandora/PandoraInternal.h"
#include "Pandora/Pandora.h"
#include "Pandora/Algorithm.h"


#include "larpandoracontent/LArShowerRefinement/LArProtoShower.h"

#include "larpandoracontent/LArHelpers/LArConnectionPathwayHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include <limits>

using namespace pandora;

namespace lar_content
{

//------------------------------------------------------------------------------------------------------------------------------------------
// Signal Definitions
//------------------------------------------------------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------------------------------------------------------
// Determining Whether ProtoShowers Are Consistent
//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::AreShowerStartsConsistent(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU, 
    const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const float maxXSeparation, const float maxSeparation)
{
    float metric(0.0);

    return LArConnectionPathwayHelper::AreShowerStartsConsistent(pAlgorithm, protoShowerU, protoShowerV, protoShowerW, maxXSeparation, 
        maxSeparation, metric);
}

//------------------------------------------------------------------------------------------------------------------------------------------
// also check start direction??
bool LArConnectionPathwayHelper::AreShowerStartsConsistent(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU, 
    const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const float maxXSeparation, const float maxSeparation, float &metric)
{
    const CartesianVector &showerStartU(protoShowerU.m_showerCore.m_startPosition);
    const CartesianVector &showerStartV(protoShowerV.m_showerCore.m_startPosition);
    const CartesianVector &showerStartW(protoShowerW.m_showerCore.m_startPosition);

    const float xSeparationUV(std::fabs(showerStartU.GetX() - showerStartV.GetX()));
    const float xSeparationUW(std::fabs(showerStartU.GetX() - showerStartW.GetX()));
    const float xSeparationVW(std::fabs(showerStartV.GetX() - showerStartW.GetX()));

    if ((xSeparationUV > maxXSeparation) || (xSeparationUW > maxXSeparation) || (xSeparationVW > maxXSeparation))
        return false;

    float chi2(0.f);
    CartesianVector projectionU(0.f, 0.f, 0.f), projectionV(0.f, 0.f, 0.f), projectionW(0.f, 0.f, 0.f);

    LArGeometryHelper::MergeTwoPositions(pAlgorithm->GetPandora(), TPC_VIEW_V, TPC_VIEW_W, showerStartV, showerStartW, projectionU, chi2);
    LArGeometryHelper::MergeTwoPositions(pAlgorithm->GetPandora(), TPC_VIEW_W, TPC_VIEW_U, showerStartW, showerStartU, projectionV, chi2);
    LArGeometryHelper::MergeTwoPositions(pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, showerStartU, showerStartV, projectionW, chi2);

    const float separationU((projectionU - showerStartU).GetMagnitude());
    const float separationV((projectionV - showerStartV).GetMagnitude());
    const float separationW((projectionW - showerStartW).GetMagnitude());

    /////////////////////////////////
    /*
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &projectionU, "U PROJECTION", RED, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &projectionV, "V PROJECTION", RED, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &projectionW, "W PROJECTION", RED, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerStartU, "U SHOWER START", BLACK, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerStartV, "V SHOWER START", BLACK, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerStartW, "W SHOWER START", BLACK, 2);
    std::cout << "separationU: " << separationU << std::endl;
    std::cout << "separationV: " << separationV << std::endl;
    std::cout << "separationW: " << separationW << std::endl;
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    */
    /////////////////////////////////

    metric = (separationU + separationV + separationW) / 3.f;

    return (metric < maxSeparation);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::AreDirectionsConsistent(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU, 
    const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const float maxOpeningAngle)
{
    float metric(0.0);

    return LArConnectionPathwayHelper::AreDirectionsConsistent(pAlgorithm, protoShowerU, protoShowerV, protoShowerW, maxOpeningAngle, metric);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::AreDirectionsConsistent(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU, 
    const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const float maxOpeningAngle, float &metric)
{
    CartesianVector directionU(protoShowerU.m_connectionPathway.m_startDirection);
    CartesianVector directionV(protoShowerV.m_connectionPathway.m_startDirection);
    CartesianVector directionW(protoShowerW.m_connectionPathway.m_startDirection);

    const CartesianVector wireAxis(0.f, 0.f, 1.f);

    float wireDeviationU(directionU.GetOpeningAngle(wireAxis));
    wireDeviationU = std::min(wireDeviationU, static_cast<float>(M_PI - wireDeviationU));

    float wireDeviationV(directionV.GetOpeningAngle(wireAxis));
    wireDeviationV = std::min(wireDeviationV, static_cast<float>(M_PI - wireDeviationV));

    float wireDeviationW(directionW.GetOpeningAngle(wireAxis));
    wireDeviationW = std::min(wireDeviationW, static_cast<float>(M_PI - wireDeviationW));

    bool isIsochronous((wireDeviationU < 2.f) && (wireDeviationV < 2.f) && (wireDeviationW < 2.f));

    if (isIsochronous)
    {
        int positiveCount(0);
        positiveCount += directionU.GetX() > 0.f ? 0 : 1;
        positiveCount += directionV.GetX() > 0.f ? 0 : 1;
        positiveCount += directionW.GetX() > 0.f ? 0 : 1;

        if (positiveCount >= 2)
        {
            directionU = CartesianVector(std::fabs(directionU.GetX()), 0.f, directionU.GetZ());
            directionV = CartesianVector(std::fabs(directionV.GetX()), 0.f, directionV.GetZ());
            directionW = CartesianVector(std::fabs(directionW.GetX()), 0.f, directionW.GetZ());
        }
    }

    if (directionU.GetX() * directionV.GetX() < 0.f)
        return false;

    if (directionU.GetX() * directionW.GetX() < 0.f)
        return false;

    if (directionV.GetX() * directionW.GetX() < 0.f)
        return false;
    
    const CartesianVector projectionU(LArGeometryHelper::MergeTwoDirections(pAlgorithm->GetPandora(), TPC_VIEW_V, TPC_VIEW_W, directionV, directionW));
    const CartesianVector projectionV(LArGeometryHelper::MergeTwoDirections(pAlgorithm->GetPandora(), TPC_VIEW_W, TPC_VIEW_U, directionW, directionU));
    const CartesianVector projectionW(LArGeometryHelper::MergeTwoDirections(pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, directionU, directionV));

    float openingAngleU(directionU.GetOpeningAngle(projectionU) * 180.f / M_PI);
    float openingAngleV(directionV.GetOpeningAngle(projectionV) * 180.f / M_PI);
    float openingAngleW(directionW.GetOpeningAngle(projectionW) * 180.f / M_PI);

    if (isIsochronous)
    {
        openingAngleU = std::min(openingAngleU, 180.f - openingAngleU);
        openingAngleV = std::min(openingAngleV, 180.f - openingAngleV);
        openingAngleW = std::min(openingAngleW, 180.f - openingAngleW);
    }

    /////////////////////////////////
    /*
    const CartesianVector &nuVertexU(protoShowerU.m_connectionPathway.m_startPosition);;
    const CartesianVector &nuVertexV(protoShowerV.m_connectionPathway.m_startPosition);
    const CartesianVector &nuVertexW(protoShowerW.m_connectionPathway.m_startPosition);
    const CartesianVector endU(nuVertexU + (directionU * 10.f));
    const CartesianVector endV(nuVertexV + (directionV * 10.f));
    const CartesianVector endW(nuVertexW + (directionW * 10.f));
    const CartesianVector projectionEndU(nuVertexU + (projectionU * 10.f));
    const CartesianVector projectionEndV(nuVertexV + (projectionV * 10.f));
    const CartesianVector projectionEndW(nuVertexW + (projectionW * 10.f));
    PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &nuVertexU, &projectionEndU, "PROJECTION U", RED, 2, 2);
    PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &nuVertexV, &projectionEndV, "PROJECTION V", RED, 2, 2);
    PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &nuVertexW, &projectionEndW, "PROJECTION W", RED, 2, 2);
    PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &nuVertexU, &endU, "DIRECTION U", BLACK, 2, 2);
    PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &nuVertexV, &endV, "DIRECTION V", BLACK, 2, 2);
    PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &nuVertexW, &endW, "DIRECTION W", BLACK, 2, 2);
    std::cout << "angularDeviationU: " << openingAngleU << std::endl;
    std::cout << "angularDeviationV: " << openingAngleV << std::endl;
    std::cout << "angularDeviationW: " << openingAngleW << std::endl;
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    */
    /////////////////////////////////

    metric = (openingAngleU + openingAngleV + openingAngleW) / 3.f;

    return (metric < maxOpeningAngle);
}

//------------------------------------------------------------------------------------------------------------------------------------------
// Finding 3D Shower Start Position
//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FindShowerVertexFromPosition(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU,
    const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, CartesianVector &showerStart3D)
{
    const CartesianVector showerStartU(protoShowerU.m_showerCore.m_startPosition);
    const CartesianVector showerStartV(protoShowerV.m_showerCore.m_startPosition);
    const CartesianVector showerStartW(protoShowerW.m_showerCore.m_startPosition);

    float chi2(0.0);
    LArGeometryHelper::MergeThreePositions3D(pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, showerStartU, showerStartV, showerStartW, showerStart3D, chi2);

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FindShowerVertexFromDirection(pandora::Algorithm *const pAlgorithm, const CartesianVector nuVertexPosition, 
    const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, CartesianVector &showerStart3D)
{
    const CartesianVector &showerStartU(protoShowerU.m_showerCore.m_startPosition);
    const CartesianVector &showerStartV(protoShowerV.m_showerCore.m_startPosition);
    const CartesianVector &showerStartW(protoShowerW.m_showerCore.m_startPosition);

    const CartesianVector &projectedNuVertexU(protoShowerU.m_connectionPathway.m_startPosition);
    const CartesianVector &projectedNuVertexV(protoShowerV.m_connectionPathway.m_startPosition);
    const CartesianVector &projectedNuVertexW(protoShowerW.m_connectionPathway.m_startPosition);

    const CartesianVector &peakDirectionU(protoShowerU.m_connectionPathway.m_startDirection);
    const CartesianVector &peakDirectionV(protoShowerV.m_connectionPathway.m_startDirection);
    const CartesianVector &peakDirectionW(protoShowerW.m_connectionPathway.m_startDirection);

    const CartesianVector &displacementU(showerStartU - projectedNuVertexU);
    const CartesianVector &displacementV(showerStartV - projectedNuVertexV);
    const CartesianVector &displacementW(showerStartW - projectedNuVertexW);

    const float transverseU(peakDirectionU.GetCrossProduct(displacementU).GetMagnitude());
    const float transverseV(peakDirectionV.GetCrossProduct(displacementV).GetMagnitude());
    const float transverseW(peakDirectionW.GetCrossProduct(displacementW).GetMagnitude());

    if ((transverseU > 1.f) || (transverseV > 1.f) || (transverseW > 1.f))
        return false;

    const CartesianVector xAxis(1.f, 0.f, 0.f);
    const float cosThetaU(peakDirectionU.GetCosOpeningAngle(xAxis));
    const float cosThetaV(peakDirectionV.GetCosOpeningAngle(xAxis));
    const float cosThetaW(peakDirectionW.GetCosOpeningAngle(xAxis));

    const float x1(showerStartU.GetX());
    const CartesianVector showerStartV1(projectedNuVertexV + (peakDirectionV * ((x1 - projectedNuVertexV.GetX()) / cosThetaV)));
    const CartesianVector showerStartW1(projectedNuVertexW + (peakDirectionW * ((x1 - projectedNuVertexW.GetX()) / cosThetaW)));

    const float x2(showerStartV.GetX());
    const CartesianVector showerStartU2(projectedNuVertexU + (peakDirectionU * ((x2 - projectedNuVertexU.GetX()) / cosThetaU)));
    const CartesianVector showerStartW2(projectedNuVertexW + (peakDirectionW * ((x2 - projectedNuVertexW.GetX()) / cosThetaW)));

    const float x3(showerStartW.GetX());
    const CartesianVector showerStartU3(projectedNuVertexU + (peakDirectionU * ((x3 - projectedNuVertexU.GetX()) / cosThetaU)));
    const CartesianVector showerStartV3(projectedNuVertexV + (peakDirectionV * ((x3 - projectedNuVertexV.GetX()) / cosThetaV)));

    float chi2(0.0);
    CartesianVector showerStart1(0.f, 0.f, 0.f), showerStart2(0.f, 0.f, 0.f), showerStart3(0.f, 0.f, 0.f);

    LArGeometryHelper::MergeThreePositions3D(pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, showerStartU, showerStartV1, showerStartW1, showerStart1, chi2);
    LArGeometryHelper::MergeThreePositions3D(pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, showerStartU2, showerStartV, showerStartW2, showerStart2, chi2);
    LArGeometryHelper::MergeThreePositions3D(pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, showerStartU3, showerStartV3, showerStartW, showerStart3, chi2);

    float minSeparation(std::numeric_limits<float>::max());

    for (const CartesianVector showerStart : {showerStart1, showerStart2, showerStart3})
    {
        const float separation((showerStart - nuVertexPosition).GetMagnitude());

        if (separation < minSeparation)
            showerStart3D = showerStart;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FindShowerVertexFromXProjection(pandora::Algorithm *const pAlgorithm, const CartesianVector nuVertexPosition, 
    const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const float maxSeparation, CartesianVector &showerStart3D)
{
    bool found(false);
    float minSeparation(std::numeric_limits<float>::max());

    for (const ProtoShower &protoShower : {protoShowerU, protoShowerV, protoShowerW})
    {
        const CartesianVector showerStart(protoShower.m_showerCore.m_startPosition);
        CartesianVector showerStart1(0.f, 0.f, 0.f), showerStart2(0.f, 0.f, 0.f);

        const HitType hitType(protoShower.m_connectionPathway.m_pathwayHitList.front()->GetHitType());
        const HitType hitType1(hitType == TPC_VIEW_U ? TPC_VIEW_V : hitType == TPC_VIEW_V ? TPC_VIEW_W : TPC_VIEW_U);
        const HitType hitType2(hitType == TPC_VIEW_U ? TPC_VIEW_W : hitType == TPC_VIEW_V ? TPC_VIEW_U : TPC_VIEW_V);

        const ProtoShower &protoShower1(hitType == TPC_VIEW_U ? protoShowerV : hitType == TPC_VIEW_V ? protoShowerW : protoShowerU);
        const ProtoShower &protoShower2(hitType == TPC_VIEW_U ? protoShowerW : hitType == TPC_VIEW_V ? protoShowerU : protoShowerV);

        bool found1(false);
        float lowestL1(std::numeric_limits<float>::max());

        for (const CaloHit *const pCaloHit : protoShower1.m_connectionPathway.m_pathwayHitList)
        {
            if (std::fabs(pCaloHit->GetPositionVector().GetX() - showerStart.GetX()) > 0.5f)
                continue;

            float lVertex(protoShower1.m_connectionPathway.m_startDirection.GetDotProduct(pCaloHit->GetPositionVector() - protoShower1.m_connectionPathway.m_startPosition));

            if ((lVertex > 0.f) && (lVertex < lowestL1))
            {
                lowestL1 = lVertex;
                showerStart1 = pCaloHit->GetPositionVector();
                found1 = true;
            }
        }

        if (!found1)
            continue;

        bool found2(false);
        float lowestL2(std::numeric_limits<float>::max());

        for (const CaloHit *const pCaloHit : protoShower2.m_connectionPathway.m_pathwayHitList)
        {
            if (std::fabs(pCaloHit->GetPositionVector().GetX() - showerStart.GetX()) > 0.5f)
                continue;

            float lVertex(protoShower2.m_connectionPathway.m_startDirection.GetDotProduct(pCaloHit->GetPositionVector() - protoShower2.m_connectionPathway.m_startPosition));

            if ((lVertex > 0.f) && (lVertex < lowestL2))
            {
                lowestL2 = lVertex;
                showerStart2 = pCaloHit->GetPositionVector();
                found2 = true;
            }
        }

        if (!found2)
            continue;

        // Now make sure that they agree...

        float chi2(0.f);
        CartesianVector projection(0.f, 0.f, 0.f), projection1(0.f, 0.f, 0.f), projection2(0.f, 0.f, 0.f);

        LArGeometryHelper::MergeTwoPositions(pAlgorithm->GetPandora(), hitType1, hitType2, showerStart1, showerStart2, projection, chi2);
        const float separationSquared((projection - showerStart).GetMagnitudeSquared());

        LArGeometryHelper::MergeTwoPositions(pAlgorithm->GetPandora(), hitType, hitType2, showerStart, showerStart2, projection1, chi2);
        const float separationSquared1((projection1 - showerStart1).GetMagnitudeSquared());

        LArGeometryHelper::MergeTwoPositions(pAlgorithm->GetPandora(), hitType, hitType1, showerStart, showerStart1, projection2, chi2);
        const float separationSquared2((projection2 - showerStart2).GetMagnitudeSquared());

    //////////////////
        /*
    std::cout << "SETTING THE NEW VERTEX" << std::endl;
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &projectionU, "U PROJECTION", RED, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &projectionV, "V PROJECTION", RED, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &projectionW, "W PROJECTION", RED, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerStartU, "U SHOWER START", BLACK, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerStartV, "V SHOWER START", BLACK, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerStartW, "W SHOWER START", BLACK, 2);
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
        */
    //////////////////
    

        if ((separationSquared > maxSeparation * maxSeparation) || (separationSquared1 > maxSeparation * maxSeparation) || 
            (separationSquared2 > maxSeparation * maxSeparation))
        {
            continue;
        }

        CartesianVector mergedShowerStart(0.f, 0.f, 0.f);
        LArGeometryHelper::MergeThreePositions3D(pAlgorithm->GetPandora(), hitType, hitType1, hitType2, showerStart, showerStart1, showerStart2, mergedShowerStart, chi2);

        const float separation((mergedShowerStart - nuVertexPosition).GetMagnitude());

        if (separation < minSeparation)
        {
            found = true;
            showerStart3D = mergedShowerStart;
        }
    }

    return found;
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content
