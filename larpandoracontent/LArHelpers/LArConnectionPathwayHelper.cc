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

#include "PandoraMonitoringApi.h"

#include "larpandoracontent/LArShowerRefinement/LArProtoShower.h"

#include "larpandoracontent/LArHelpers/LArConnectionPathwayHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"


#include <limits>

using namespace pandora;

namespace lar_content
{



//------------------------------------------------------------------------------------------------------------------------------------------
// Electron Tree
//------------------------------------------------------------------------------------------------------------------------------------------

LArConnectionPathwayHelper::ElectronTreeVariables::ElectronTreeVariables() :
    m_pathwayLength(-10.f),
    m_pathwayShowerStartDelta(-10.f),
    m_pathwayMaxScatteringAngle(-10.f),
    m_pathwayEnergyMeanU(-10.f),
    m_pathwayEnergyMeanV(-10.f),
    m_pathwayEnergyMeanW(-10.f),
    m_pathwayEnergySigmaU(-10.f),
    m_pathwayEnergySigmaV(-10.f),
    m_pathwayEnergySigmaW(-10.f),
    m_pathwayMinEnergyMean(-10.f),
    m_pathwayMiddleEnergyMean(-10.f),
    m_pathwayMaxEnergyMean(-10.f),
    m_pathwayMinEnergySigma(-10.f),
    m_pathwayMiddleEnergySigma(-10.f),
    m_pathwayMaxEnergySigma(-10.f),
    m_postShowerStartLength(-10.f),
    m_postShowerStartNHits(-10),
    m_postShowerStartScatterAngle(-10.f),
    m_postShowerStartMeanTransverseAngle(-10.f),
    m_postShowerStartMeanRadialDistance(-10.f),
    m_postShowerStartRadialDistanceSigma(-10.f),
    m_postShowerStartEnergyWeightedMeanRadialDistance(-10.f),
    m_postShowerStartEstimatedMoliereRadius(-10.f),
    m_postShowerStartInitialGapSize(-10.f),
    m_postShowerStartMaxGapSize(-10.f),
    m_initialRegionDistanceToNuVertex(-10.f),
    m_initialRegionDistanceInGaps(-10.f),
    m_initialRegionMaxGapSize(-10.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArConnectionPathwayHelper::FillElectronTreeVariables(pandora::Algorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo, 
    const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const CartesianVector &nuVertexPosition, 
    const LArConnectionPathwayHelper::Consistency &consistency, LArConnectionPathwayHelper::ElectronTreeVariables &electronTreeVariables)
{
    CartesianVector minShowerStart3D(0.f, 0.f, 0.f), middleShowerStart3D(0.f, 0.f, 0.f), maxShowerStart3D(0.f, 0.f, 0.f);

    if (!LArConnectionPathwayHelper::FindShowerStarts3D(pAlgorithm, protoShowerU, protoShowerV, protoShowerW, consistency, nuVertexPosition, 
        minShowerStart3D, middleShowerStart3D, maxShowerStart3D))
    {
        std::cout << "middleShowerStart3D: " << middleShowerStart3D << std::endl;

        std::cout << "COULD NOT FIND ANY SHOWER VERTICES!" << std::endl;
        return;
    }
   /////////////////////////////////

    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &minShowerStart3D, "MIN", BLACK, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &middleShowerStart3D, "MIDDLE", RED, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &maxShowerStart3D, "MAX", BLUE, 2);
    PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &protoShowerU.m_spineHitList, "U SPINE", GREEN);
    PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &protoShowerV.m_spineHitList, "V SPINE", GREEN);
    PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &protoShowerW.m_spineHitList, "W SPINE", GREEN);
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());

    /////////////////////////////////

    LArConnectionPathwayHelper::FillConnectionPathwayVariables(pAlgorithm, protoShowerU, protoShowerV, protoShowerW, nuVertexPosition, minShowerStart3D,
        middleShowerStart3D, maxShowerStart3D, electronTreeVariables);

    std::cout << "middleShowerStart3D: " << middleShowerStart3D << std::endl;

    LArConnectionPathwayHelper::FillPostShowerStartVariables(pAlgorithm, pShowerPfo, nuVertexPosition, protoShowerU, protoShowerV, protoShowerW, middleShowerStart3D, electronTreeVariables);

    LArConnectionPathwayHelper::FillInitialRegionVariables(pAlgorithm, nuVertexPosition, protoShowerU, protoShowerV, protoShowerW, 
        middleShowerStart3D, electronTreeVariables);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FillAmbiguousHitVariables(pandora::Algorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo, 
    const ElectronProtoShower &protoShowerU, const ElectronProtoShower &protoShowerV, const ElectronProtoShower &protoShowerW, const CartesianVector &nuVertexPosition,
    CartesianVector &middleShowerStart3D, const TwoDSlidingFitResult &spineFitU, const TwoDSlidingFitResult &spineFitV, const TwoDSlidingFitResult &spineFitW,
    const CaloHitList *const pCaloHitListU, const CaloHitList *const pCaloHitListV, const CaloHitList *const pCaloHitListW)
{
    int nViewsWithAmbiguousHits(0);

    const int nAmbiguousHitsU(protoShowerU.m_ambiguousHitList.size());
    nViewsWithAmbiguousHits += (nAmbiguousHitsU == 0) ? 0 : 1;

    const int nAmbiguousHitsV(protoShowerV.m_ambiguousHitList.size());
    nViewsWithAmbiguousHits += (nAmbiguousHitsV == 0) ? 0 : 1;

    const int nAmbiguousHitsW(protoShowerW.m_ambiguousHitList.size());
    nViewsWithAmbiguousHits += (nAmbiguousHitsW == 0) ? 0 : 1;

    const int nAmbiguousHits(nAmbiguousHitsU + nAmbiguousHitsV + nAmbiguousHitsW);
    std::cout << "nAmbiguousHits: " << nAmbiguousHits << std::endl;

    const int minNAmbiguousHits(std::min(std::min(nAmbiguousHitsU, nAmbiguousHitsV), nAmbiguousHitsW));
    std::cout << "minNAmbiguousHits: " << minNAmbiguousHits << std::endl;

    const int maxNAmbiguousHits(std::max(std::max(nAmbiguousHitsU, nAmbiguousHitsV), nAmbiguousHitsW));
    std::cout << "maxNAmbiguousHits: " << maxNAmbiguousHits << std::endl;

    const bool isDownstream(middleShowerStart3D.GetZ() > nuVertexPosition.GetZ());
    const CartesianVector &startDirectionU(isDownstream ? spineFitU.GetGlobalMinLayerDirection() : spineFitU.GetGlobalMaxLayerDirection() * (-1.f));
    const CartesianVector &startDirectionV(isDownstream ? spineFitV.GetGlobalMinLayerDirection() : spineFitV.GetGlobalMaxLayerDirection() * (-1.f));
    const CartesianVector &startDirectionW(isDownstream ? spineFitW.GetGlobalMinLayerDirection() : spineFitW.GetGlobalMaxLayerDirection() * (-1.f));

    CartesianVector connectionPathwayDirection3D(0.f, 0.f, 0.f);

    if (!LArConnectionPathwayHelper::FindDirection3D(pAlgorithm, isDownstream, protoShowerU.m_connectionPathway.m_startPosition, protoShowerV.m_connectionPathway.m_startPosition,
        protoShowerW.m_connectionPathway.m_startPosition, nuVertexPosition, startDirectionU, startDirectionV, startDirectionW, connectionPathwayDirection3D))
    {
        std::cout << "CANT FIND A DIRECTION FOR THE START OF THE CONNECTION PATHWAY" << std::endl;
        return false;
    }

    float unaccountedHitEnergyU(-10.f), showerEnergyRatioU(-10.f);
    LArConnectionPathwayHelper::GetViewAmbiguousHitVariables(pAlgorithm, protoShowerU, TPC_VIEW_U, nuVertexPosition, connectionPathwayDirection3D, 
        pCaloHitListU, unaccountedHitEnergyU, showerEnergyRatioU);

    float unaccountedHitEnergyV(-10.f), showerEnergyRatioV(-10.f);
    LArConnectionPathwayHelper::GetViewAmbiguousHitVariables(pAlgorithm, protoShowerV, TPC_VIEW_V, nuVertexPosition, connectionPathwayDirection3D, 
        pCaloHitListV, unaccountedHitEnergyV, showerEnergyRatioV);

    float unaccountedHitEnergyW(-10.f), showerEnergyRatioW(-10.f);
    LArConnectionPathwayHelper::GetViewAmbiguousHitVariables(pAlgorithm, protoShowerW, TPC_VIEW_W, nuVertexPosition, connectionPathwayDirection3D, 
        pCaloHitListW, unaccountedHitEnergyW, showerEnergyRatioW);

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArConnectionPathwayHelper::GetViewAmbiguousHitVariables(pandora::Algorithm *const pAlgorithm, const ElectronProtoShower &protoShower, 
    const HitType hitType, const CartesianVector &nuVertexPosition, const CartesianVector &connectionPathwayDirection3D,
    const CaloHitList *const pCaloHitList, float &unaccountedHitEnergy, float &showerEnergyRatio)
{
    const CartesianVector projectedNuVertex(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), nuVertexPosition, hitType));

    // Find the owner connection pathways of the ambiguous hits
    std::map<int, CaloHitList> ambiguousHitSpinesTemp;
    CaloHitList hitsToExcludeInEnergyCalcs; // to avoid double counting the energy

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        if (std::find(protoShower.m_ambiguousHitList.begin(), protoShower.m_ambiguousHitList.end(), pCaloHit) != protoShower.m_ambiguousHitList.end())
            continue;

        int count(0);

        // A hit can be in more than one spine
        for (unsigned int i = 0; i < protoShower.m_ambiguousDirectionVector.size(); ++i)
        {
            const CartesianVector &significantDirection(protoShower.m_ambiguousDirectionVector[i]);
            const CartesianVector displacement(pCaloHit->GetPositionVector() - projectedNuVertex);
            const float thisT(significantDirection.GetCrossProduct(displacement).GetMagnitude());
            const float thisL(significantDirection.GetDotProduct(displacement));

            if ((thisL > 0.f) && (thisT < 0.75f))
            {
                ++count;
                ambiguousHitSpinesTemp[i].push_back(pCaloHit);
            }

            if (count == 2)
                hitsToExcludeInEnergyCalcs.push_back(pCaloHit);
        }
    }

    if (ambiguousHitSpinesTemp.empty())
        return;

    // Find continuous pathways i.e. make sure other connection pathway objects are decent
    std::map<int, CaloHitList> ambiguousHitSpines;

    for (const auto &entry : ambiguousHitSpinesTemp)
    {
        CaloHitList continuousSpine(LArConnectionPathwayHelper::FindAmbiguousContinuousSpine(entry.second, protoShower.m_ambiguousHitList, projectedNuVertex));

        if (continuousSpine.size() > 0)
            ambiguousHitSpines[entry.first] = continuousSpine;
    }

    if (ambiguousHitSpines.empty())
        return;

    // We can't do a de/dx without a 3D direction :( 
    // Which is really hard to obtain for the other particles in the event... so we might have to do mean energy? 
    // But we can make sure that we only consider a certain length of true spine (lets say 5cm)
    // Pick this to avoid going into the shower i.e. hopefully shower dedx will be a bit constant

    const CartesianVector seedPoint3D(nuVertexPosition + (connectionPathwayDirection3D * 5.f));
    const CartesianVector seedPoint(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), seedPoint3D, hitType));
    const CartesianVector startDirection(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), connectionPathwayDirection3D, hitType));
    const float lRange(startDirection.GetDotProduct(seedPoint - projectedNuVertex));

    // Need to find max start position...
    float ambiguousHitEnergyMean(0.f);
    float startL(-std::numeric_limits<float>::max());

    for (const CaloHit *const pAmbiguousCaloHit : protoShower.m_ambiguousHitList)
    {
        const float thisT(startDirection.GetCrossProduct(pAmbiguousCaloHit->GetPositionVector() - projectedNuVertex).GetMagnitude());
        const float thisL(startDirection.GetDotProduct(pAmbiguousCaloHit->GetPositionVector() - projectedNuVertex));

        if ((thisL > startL) && (thisT < 0.75f))
            startL = thisL;

        ambiguousHitEnergyMean += pAmbiguousCaloHit->GetElectromagneticEnergy() * 1000.f;
    }

    if (startL < 0.f)
        return;

    ambiguousHitEnergyMean /= protoShower.m_ambiguousHitList.size();

    // Get mean energy of other pathways, avoiding the double counting hits (ambiguous to the ambiguous hits)
    float otherEnergyMeanSum(0.f);

    for (const auto &entry : ambiguousHitSpines)
    {
        int nOtherEnergyHits(0);
        float otherEnergyMean(0.f);

        for (const CaloHit *const pOtherCaloHit : entry.second)
        {
            if (std::find(hitsToExcludeInEnergyCalcs.begin(), hitsToExcludeInEnergyCalcs.end(), pOtherCaloHit) != hitsToExcludeInEnergyCalcs.end())
                continue;

            otherEnergyMean += pOtherCaloHit->GetElectromagneticEnergy() * 1000.f;
            ++nOtherEnergyHits;
        }

        if (nOtherEnergyHits == 0)
            continue;

        otherEnergyMean /= static_cast<float>(nOtherEnergyHits);
        otherEnergyMeanSum += otherEnergyMean;
    }

    // Get the spine mean energy

    float spineEnergyMean(0.f);
    int nSpineEnergyHits(0);

    for (const CaloHit *const pSpineCaloHit : protoShower.m_spineHitList)
    {
        if (std::find(protoShower.m_ambiguousHitList.begin(), protoShower.m_ambiguousHitList.end(), pSpineCaloHit) != protoShower.m_ambiguousHitList.end())
            continue;

        const float thisL(startDirection.GetDotProduct(pSpineCaloHit->GetPositionVector() - projectedNuVertex));

        if ((thisL > startL) && (thisL < (startL + lRange)))
        {
            spineEnergyMean += pSpineCaloHit->GetElectromagneticEnergy() * 1000.f;
            ++nSpineEnergyHits;
        }
    }

    if (nSpineEnergyHits == 0)
        return;

    spineEnergyMean /= static_cast<float>(nSpineEnergyHits);

    unaccountedHitEnergy = ambiguousHitEnergyMean - otherEnergyMeanSum - spineEnergyMean;
    std::cout << "unaccountedHitEnergy: " << unaccountedHitEnergy << std::endl;
    showerEnergyRatio = (ambiguousHitEnergyMean - otherEnergyMeanSum) / spineEnergyMean;
    std::cout << "showerEnergyRatio: " << showerEnergyRatio << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CaloHitList LArConnectionPathwayHelper::FindAmbiguousContinuousSpine(const CaloHitList &caloHitList, const CaloHitList &ambiguousHitList, 
    const CartesianVector &projectedNuVertexPosition)
{
    CaloHitList continuousHitList;

    CaloHitVector caloHitVector(caloHitList.begin(), caloHitList.end());
    std::sort(caloHitVector.begin(), caloHitVector.end(), LArConnectionPathwayHelper::SortByDistanceToPoint(projectedNuVertexPosition));

    for (unsigned int i = 0; i < caloHitVector.size(); ++i)
    {
        CaloHitList connectedHitList;
        connectedHitList.push_back(caloHitVector[i]);

        if (LArClusterHelper::GetClosestDistance(connectedHitList.front()->GetPositionVector(), ambiguousHitList) > 1.f)
            continue;

        bool found(true);

        while(found)
        {
            found = false;

            for (unsigned int j = (i + 1); j < caloHitVector.size(); ++j)
            {
                const CaloHit *const pCaloHit(caloHitVector[j]);

                if (std::find(connectedHitList.begin(), connectedHitList.end(), pCaloHit) != connectedHitList.end())
                    continue;

                if (LArClusterHelper::GetClosestDistance(pCaloHit->GetPositionVector(), connectedHitList) < 1.f)
                {
                    found = true;
                    connectedHitList.push_back(pCaloHit);
                    break;
                }
            }
        }

        if ((connectedHitList.size() >= 2) || (connectedHitList.size() == caloHitVector.size()))
        {
            continuousHitList.insert(continuousHitList.begin(), connectedHitList.begin(), connectedHitList.end());
            break;
        }
    }

    return continuousHitList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FillInitialRegionVariables(pandora::Algorithm *const pAlgorithm, const CartesianVector &nuVertexPosition, 
    const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, 
    CartesianVector &middleShowerStart3D, LArConnectionPathwayHelper::ElectronTreeVariables &electronTreeVariables)
{
    const bool isDownstream(middleShowerStart3D.GetZ() > nuVertexPosition.GetZ());

    CartesianPointVector spinePositionsU, spinePositionsV, spinePositionsW;

    for (const CaloHit *const pCaloHit : protoShowerU.m_spineHitList)
        spinePositionsU.push_back(pCaloHit->GetPositionVector());

    for (const CaloHit *const pCaloHit : protoShowerV.m_spineHitList)
        spinePositionsV.push_back(pCaloHit->GetPositionVector());

    for (const CaloHit *const pCaloHit : protoShowerW.m_spineHitList)
        spinePositionsW.push_back(pCaloHit->GetPositionVector());

    try
    {
        // try 10 so i can make it consistent?
        const TwoDSlidingFitResult spineFitU(&spinePositionsU, 10, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));
        const TwoDSlidingFitResult spineFitV(&spinePositionsV, 10, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));
        const TwoDSlidingFitResult spineFitW(&spinePositionsW, 10, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));

        const CartesianVector &startDirectionU(isDownstream ? spineFitU.GetGlobalMinLayerDirection() : spineFitU.GetGlobalMaxLayerDirection() * (-1.f));
        const CartesianVector &startDirectionV(isDownstream ? spineFitV.GetGlobalMinLayerDirection() : spineFitV.GetGlobalMaxLayerDirection() * (-1.f));
        const CartesianVector &startDirectionW(isDownstream ? spineFitW.GetGlobalMinLayerDirection() : spineFitW.GetGlobalMaxLayerDirection() * (-1.f));

        CartesianVector connectionPathwayDirection3D(0.f, 0.f, 0.f);

        if (!LArConnectionPathwayHelper::FindDirection3D(pAlgorithm, isDownstream, protoShowerU.m_connectionPathway.m_startPosition, protoShowerV.m_connectionPathway.m_startPosition,
            protoShowerW.m_connectionPathway.m_startPosition, nuVertexPosition, startDirectionU, startDirectionV, startDirectionW, connectionPathwayDirection3D))
        {
            std::cout << "CANT FIND A DIRECTION FOR THE START OF THE CONNECTION PATHWAY" << std::endl;
            return false;
        }

        const CartesianVector projectedNuVertexU(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), nuVertexPosition, TPC_VIEW_U));
        const CartesianVector projectedNuVertexV(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), nuVertexPosition, TPC_VIEW_V));
        const CartesianVector projectedNuVertexW(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), nuVertexPosition, TPC_VIEW_W));

        // The 3D hits of the connecting pathway might not exist -.- so we are going to have to try and make them...

        std::map<float, const CaloHit*> longitudinalProjectionsU, longitudinalProjectionsV, longitudinalProjectionsW;

        for (const CaloHit *const pCaloHit : protoShowerU.m_spineHitList)
            longitudinalProjectionsU[startDirectionU.GetDotProduct(pCaloHit->GetPositionVector() - projectedNuVertexU)] = pCaloHit;

        for (const CaloHit *const pCaloHit : protoShowerV.m_spineHitList)
            longitudinalProjectionsV[startDirectionV.GetDotProduct(pCaloHit->GetPositionVector() - projectedNuVertexV)] = pCaloHit;

        for (const CaloHit *const pCaloHit : protoShowerW.m_spineHitList)
            longitudinalProjectionsW[startDirectionW.GetDotProduct(pCaloHit->GetPositionVector() - projectedNuVertexW)] = pCaloHit;

        // maybe we want to limit this? to be like 10cm?
        const float length(std::min(10.f, (middleShowerStart3D - nuVertexPosition).GetMagnitude()));

        bool inPathway(true);
        float stepSize(0.5f);
        unsigned int count(0);

        // search for the closest point
        bool foundClosestPoint(false);
        CartesianVector closestPosition3D(0.f, 0.f, 0.f);
        const CaloHit *pBestHitU(nullptr), *pBestHitV(nullptr), *pBestHitW(nullptr);

        // search for any gaps
        float gapSize(0.f);
        float distanceInGaps(0.f);
        float largestGapSize(-std::numeric_limits<float>::max());

        while (inPathway)
        {
            const CartesianVector lowerBoundary(nuVertexPosition + (connectionPathwayDirection3D * count * stepSize));
            const CartesianVector lowerBoundaryU(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), lowerBoundary, TPC_VIEW_U));
            const CartesianVector lowerBoundaryV(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), lowerBoundary, TPC_VIEW_V));
            const CartesianVector lowerBoundaryW(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), lowerBoundary, TPC_VIEW_W));

            const float uLowerBoundaryProjectionL(startDirectionU.GetDotProduct(lowerBoundaryU - projectedNuVertexU));
            const float vLowerBoundaryProjectionL(startDirectionV.GetDotProduct(lowerBoundaryV - projectedNuVertexV));
            const float wLowerBoundaryProjectionL(startDirectionW.GetDotProduct(lowerBoundaryW - projectedNuVertexW));

            const CartesianVector upperBoundary(nuVertexPosition + (connectionPathwayDirection3D * ((count + 1) * stepSize)));
            const CartesianVector upperBoundaryU(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), upperBoundary, TPC_VIEW_U));
            const CartesianVector upperBoundaryV(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), upperBoundary, TPC_VIEW_V));
            const CartesianVector upperBoundaryW(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), upperBoundary, TPC_VIEW_W));

            const float uUpperBoundaryProjectionL(startDirectionU.GetDotProduct(upperBoundaryU - projectedNuVertexU));
            const float vUpperBoundaryProjectionL(startDirectionV.GetDotProduct(upperBoundaryV - projectedNuVertexV));
            const float wUpperBoundaryProjectionL(startDirectionW.GetDotProduct(upperBoundaryW - projectedNuVertexW));

            bool foundHit3D(false);
            float lowestL(std::numeric_limits<float>::max());

            for (const auto &entryU : longitudinalProjectionsU)
            {
                if ((entryU.first > uLowerBoundaryProjectionL) && (entryU.first < uUpperBoundaryProjectionL))
                {
                    for (const auto &entryV : longitudinalProjectionsV)
                    {
                        if ((entryV.first > vLowerBoundaryProjectionL) && (entryV.first < vUpperBoundaryProjectionL))
                        {
                            for (const auto &entryW : longitudinalProjectionsW)
                            {
                                if ((entryW.first > wLowerBoundaryProjectionL) && (entryW.first < wUpperBoundaryProjectionL))
                                 {
                                     float metric(0.f);
                                     if (LArConnectionPathwayHelper::AreShowerStartsConsistent(pAlgorithm, entryU.second->GetPositionVector(), 
                                         entryV.second->GetPositionVector(), entryW.second->GetPositionVector(), 5.f, 2.f, metric))
                                     {
                                         foundHit3D = true;

                                         CartesianVector position3D(0.f, 0.f, 0.f);

                                         LArGeometryHelper::MergeThreePositions3D(pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, entryU.second->GetPositionVector(), 
                                             entryV.second->GetPositionVector(), entryW.second->GetPositionVector(), position3D, metric);

                                         const float thisL(connectionPathwayDirection3D.GetDotProduct(position3D - nuVertexPosition));

                                         if (!foundClosestPoint && (thisL < lowestL))
                                         {
                                             lowestL = thisL;
                                             pBestHitU = entryU.second; pBestHitV = entryV.second; pBestHitW = entryW.second;
                                             closestPosition3D = position3D;
                                         }
                                     }
                                 }
                            }
                        }
                    }
                }
            }

            if (!foundHit3D)
            {
                gapSize += stepSize;
            }
            else
            {
                if (gapSize > 1.f)
                    distanceInGaps += gapSize;

                if (gapSize > largestGapSize)
                    largestGapSize = gapSize;

                gapSize = 0.f;
            }

            if (!foundClosestPoint && foundHit3D)
            {
                const CartesianVector &closestHitU(pBestHitU->GetPositionVector());
                const CartesianVector &closestHitV(pBestHitV->GetPositionVector());
                const CartesianVector &closestHitW(pBestHitW->GetPositionVector());

                PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &closestHitU, "closestHitU", BLUE, 2);
                PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &closestHitV, "closestHitV", BLUE, 2);
                PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &closestHitW, "closestHitW", BLUE, 2);
                PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &closestPosition3D, "closestPosition3D", RED, 2);
                PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
                foundClosestPoint = true;
            }

            if ((upperBoundary - nuVertexPosition).GetMagnitude() > length)
                inPathway = false;

            ++count;
        }

        if (foundClosestPoint)
        {
            const float distanceToNuVertex((closestPosition3D - nuVertexPosition).GetMagnitude());
            electronTreeVariables.m_initialRegionDistanceToNuVertex = distanceToNuVertex;
            std::cout << "electronTreeVariables.m_initialRegionDistanceToNuVertex: " << electronTreeVariables.m_initialRegionDistanceToNuVertex << std::endl;
        }

        electronTreeVariables.m_initialRegionDistanceInGaps = distanceInGaps;
        std::cout << "electronTreeVariables.m_initialRegionDistanceInGaps: " << electronTreeVariables.m_initialRegionDistanceInGaps << std::endl;
        electronTreeVariables.m_initialRegionMaxGapSize = largestGapSize;
        std::cout << "electronTreeVariables.m_initialRegionMaxGapSize: " << electronTreeVariables.m_initialRegionMaxGapSize << std::endl;
    }
    catch (...)
    {
        return false;
    }

    return true;
}


//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FillPostShowerStartVariables(pandora::Algorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo, 
    const CartesianVector &nuVertexPosition, const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, 
    CartesianVector &middleShowerStart3D, LArConnectionPathwayHelper::ElectronTreeVariables &electronTreeVariables)
{
    const bool isDownstream(middleShowerStart3D.GetZ() > nuVertexPosition.GetZ());

    CartesianPointVector spinePositionsU, spinePositionsV, spinePositionsW;

    for (const CaloHit *const pCaloHit : protoShowerU.m_spineHitList)
        spinePositionsU.push_back(pCaloHit->GetPositionVector());

    for (const CaloHit *const pCaloHit : protoShowerV.m_spineHitList)
        spinePositionsV.push_back(pCaloHit->GetPositionVector());

    for (const CaloHit *const pCaloHit : protoShowerW.m_spineHitList)
        spinePositionsW.push_back(pCaloHit->GetPositionVector());

    ///////////////////////////////////////////////////
    PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &protoShowerU.m_spineHitList, "U SPINE", GREEN);
    PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &protoShowerV.m_spineHitList, "V SPINE", GREEN);
    PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &protoShowerW.m_spineHitList, "W SPINE", GREEN);
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    ///////////////////////////////////////////////////

    try
    {
        // try 10 so i can make it consistent?
        const TwoDSlidingFitResult spineFitU(&spinePositionsU, 10, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));
        const TwoDSlidingFitResult spineFitV(&spinePositionsV, 10, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));
        const TwoDSlidingFitResult spineFitW(&spinePositionsW, 10, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));

        ///////////////////////////////////////////////////
        const LayerFitResultMap &layerFitResultMapU(spineFitU.GetLayerFitResultMap());
        const LayerFitResultMap &layerFitResultMapV(spineFitV.GetLayerFitResultMap());
        const LayerFitResultMap &layerFitResultMapW(spineFitW.GetLayerFitResultMap());

        const CartesianVector &startDirectionU(isDownstream ? spineFitU.GetGlobalMinLayerDirection() : spineFitU.GetGlobalMaxLayerDirection());
        const CartesianVector &startDirectionV(isDownstream ? spineFitV.GetGlobalMinLayerDirection() : spineFitV.GetGlobalMaxLayerDirection());
        const CartesianVector &startDirectionW(isDownstream ? spineFitW.GetGlobalMinLayerDirection() : spineFitW.GetGlobalMaxLayerDirection());

        CartesianVector connectionPathwayDirection3D(0.f, 0.f, 0.f);

        if (!LArConnectionPathwayHelper::FindDirection3D(pAlgorithm, isDownstream, protoShowerU.m_connectionPathway.m_startPosition, protoShowerV.m_connectionPathway.m_startPosition,
            protoShowerW.m_connectionPathway.m_startPosition, nuVertexPosition, startDirectionU, startDirectionV, startDirectionW, connectionPathwayDirection3D))
        {
            std::cout << "CANT FIND A DIRECTION FOR THE START OF THE CONNECTION PATHWAY" << std::endl;
            return false;
        }

        std::cout << "/////////////////////////////////" << std::endl;
        std::cout << "looking at the U fit..." << std::endl;

        for (auto &entry : layerFitResultMapU)
        {
            CartesianVector position(0.f, 0.f, 0.f);
            spineFitU.GetGlobalFitPosition(entry.second.GetL(), position);

            PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &position, "U spine fit layer", GREEN, 2);
        }

        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());

        std::cout << "/////////////////////////////////" << std::endl;
        std::cout << "looking at the V fit..." << std::endl;

        for (auto &entry : layerFitResultMapV)
        {
            CartesianVector position(0.f, 0.f, 0.f);
            spineFitV.GetGlobalFitPosition(entry.second.GetL(), position);

            PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &position, "V spine fit layer", GREEN, 2);
        }

        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());

        std::cout << "/////////////////////////////////" << std::endl;
        std::cout << "looking at the W fit..." << std::endl;

        for (auto &entry : layerFitResultMapW)
        {
            CartesianVector position(0.f, 0.f, 0.f);
            spineFitW.GetGlobalFitPosition(entry.second.GetL(), position);

            PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &position, "W spine fit layer", GREEN, 2);
        }

        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
        ///////////////////////////////////////////////////

        const CartesianVector projectedShowerStartU(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), middleShowerStart3D, TPC_VIEW_U));
        const CartesianVector projectedShowerStartV(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), middleShowerStart3D, TPC_VIEW_V));
        const CartesianVector projectedShowerStartW(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), middleShowerStart3D, TPC_VIEW_W));

        // attempt to get a 3D direction...
        float lShowerStartU(0.f), lShowerStartV(0.f), lShowerStartW(0.f);
        float tShowerStartU(0.f), tShowerStartV(0.f), tShowerStartW(0.f);

        spineFitU.GetLocalPosition(projectedShowerStartU, lShowerStartU, tShowerStartU);
        spineFitV.GetLocalPosition(projectedShowerStartV, lShowerStartV, tShowerStartV);
        spineFitW.GetLocalPosition(projectedShowerStartW, lShowerStartW, tShowerStartW);

        const int showerStartLayerU(spineFitU.GetLayer(lShowerStartU));
        const int showerStartLayerV(spineFitV.GetLayer(lShowerStartV));
        const int showerStartLayerW(spineFitW.GetLayer(lShowerStartW));

        if ((layerFitResultMapU.find(showerStartLayerU) == layerFitResultMapU.end()) || (layerFitResultMapV.find(showerStartLayerV) == layerFitResultMapV.end()) ||
            (layerFitResultMapW.find(showerStartLayerW) == layerFitResultMapW.end()))
        {
            std::cout << "shower start layer is not in layer fit result map" << std::endl;
            return false;
        }

        const float gradientU(layerFitResultMapU.at(showerStartLayerU).GetGradient());
        const float gradientV(layerFitResultMapV.at(showerStartLayerV).GetGradient());
        const float gradientW(layerFitResultMapW.at(showerStartLayerW).GetGradient());
        CartesianVector showerDirectionU(0.f, 0.f, 0.f), showerDirectionV(0.f, 0.f, 0.f), showerDirectionW(0.f, 0.f, 0.f);

        spineFitU.GetGlobalDirection(gradientU, showerDirectionU);
        spineFitV.GetGlobalDirection(gradientV, showerDirectionV);
        spineFitW.GetGlobalDirection(gradientW, showerDirectionW);

        ////////////////////////////////////
        const CartesianVector endU(projectedShowerStartU + (showerDirectionU * 100.f));
        PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &projectedShowerStartU, &endU, "PROJECTION U", RED, 2, 2);
        PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &protoShowerU.m_spineHitList, "U SPINE1", GREEN);
        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());

        const CartesianVector endV(projectedShowerStartV + (showerDirectionV * 100.f));
        PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &projectedShowerStartV, &endV, "PROJECTION V", RED, 2, 2);
        PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &protoShowerV.m_spineHitList, "V SPINE1", GREEN);
        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());

        const CartesianVector endW(projectedShowerStartW + (showerDirectionW * 100.f));
        PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &projectedShowerStartW, &endW, "PROJECTION W", RED, 2, 2);
        PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &protoShowerW.m_spineHitList, "W SPINE1", GREEN);
        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
        ////////////////////////////////////

        CartesianVector postShowerStartDirection3D(0.f, 0.f, 0.f);

        if (!LArConnectionPathwayHelper::FindDirection3D(pAlgorithm, isDownstream, projectedShowerStartU, projectedShowerStartV, projectedShowerStartW, middleShowerStart3D,
            showerDirectionU, showerDirectionV, showerDirectionW, postShowerStartDirection3D))
        {
            std::cout << "CANNOT FIND 3D POST SHOWER START DIRECTION..." << std::endl;
            return false;
        }

        CaloHitList caloHitList3D;
        LArPfoHelper::GetCaloHits(pShowerPfo, TPC_3D, caloHitList3D);

        CaloHitList postShowerHitList3D; // Do i use all hits?? or all hits in event?? o.O 
        CartesianPointVector postShowerPositions3D;

        for (const CaloHit *const pCaloHit : caloHitList3D)
        {
            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            const CartesianVector displacement(hitPosition - middleShowerStart3D);
            const float l(postShowerStartDirection3D.GetDotProduct(displacement));
            const float t(postShowerStartDirection3D.GetCrossProduct(displacement).GetMagnitude());

            if (((isDownstream && (l > 0.f)) || (!isDownstream && (l < 0.f))) && (t < 14.f))
            {
                postShowerHitList3D.push_back(pCaloHit);
                postShowerPositions3D.push_back(pCaloHit->GetPositionVector());
            }
        }

        PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &postShowerHitList3D, "postShowerHitList3D", RED);
        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());

        ThreeDSlidingFitResult slidingFitResult3D(&postShowerPositions3D, 1000, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));

        for (const CaloHit *const pCaloHit : caloHitList3D)
        {
            if (std::find(postShowerHitList3D.begin(), postShowerHitList3D.end(), pCaloHit) != postShowerHitList3D.end())
                continue;

            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            const CartesianVector displacement(hitPosition - middleShowerStart3D);
            const float l(slidingFitResult3D.GetGlobalMinLayerDirection().GetDotProduct(displacement));
            const float t(slidingFitResult3D.GetGlobalMinLayerDirection().GetCrossProduct(displacement).GetMagnitude());

            if (((isDownstream && (l > 0.f)) || (!isDownstream && (l < 0.f))) && (t < 14.f))
            {
                postShowerHitList3D.push_back(pCaloHit);
                postShowerPositions3D.push_back(pCaloHit->GetPositionVector());
            }
        }

        PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &postShowerHitList3D, "postShowerHitList3D", RED);
        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());

        ThreeDSlidingFitResult fullSlidingFitResult3D(&postShowerPositions3D, 1000, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));

        const float length((fullSlidingFitResult3D.GetGlobalMaxLayerPosition() - fullSlidingFitResult3D.GetGlobalMinLayerPosition()).GetMagnitude());
        electronTreeVariables.m_postShowerStartLength = length;
        std::cout << "electronTreeVariables.m_postShowerStartLength: " << electronTreeVariables.m_postShowerStartLength << std::endl;

        const int nPostShowerHits(postShowerHitList3D.size());
        electronTreeVariables.m_postShowerStartNHits = nPostShowerHits;
        std::cout << "electronTreeVariables.m_postShowerStartNHits: " << electronTreeVariables.m_postShowerStartNHits << std::endl;

        const float pathwayShowerAngle(fullSlidingFitResult3D.GetGlobalMinLayerDirection().GetOpeningAngle(connectionPathwayDirection3D) * 180.f / M_PI);
        electronTreeVariables.m_postShowerStartScatterAngle = pathwayShowerAngle;
        std::cout << "electronTreeVariables.m_postShowerStartScatterAngle: " << electronTreeVariables.m_postShowerStartScatterAngle << std::endl; 
 
        // Me checking the coordinate system of the 3D sliding fit
        std::cout << "axis direction from 3D: " << fullSlidingFitResult3D.GetAxisDirection() << std::endl;
        std::cout << "axis direction from first fit: " << fullSlidingFitResult3D.GetFirstFitResult().GetAxisDirection() << std::endl;
        std::cout << "ortho axis direction from first fit: " << fullSlidingFitResult3D.GetFirstFitResult().GetOrthoDirection() << std::endl;
        std::cout << "axis direction from second fit: " << fullSlidingFitResult3D.GetSecondFitResult().GetAxisDirection() << std::endl;
        std::cout << "ortho axis direction from second fit: " << fullSlidingFitResult3D.GetSecondFitResult().GetOrthoDirection() << std::endl;

        std::cout << "axis intercept: " << fullSlidingFitResult3D.GetAxisIntercept() << std::endl;
        std::cout << "first axis intercept: " << fullSlidingFitResult3D.GetFirstFitResult().GetAxisIntercept() << std::endl;
        std::cout << "second axis intercept: " << fullSlidingFitResult3D.GetSecondFitResult().GetAxisIntercept() << std::endl;

        std::cout << "3D axis with first ortho: " << fullSlidingFitResult3D.GetAxisDirection().GetDotProduct(fullSlidingFitResult3D.GetFirstFitResult().GetOrthoDirection()) << std::endl;
        std::cout << "3D axis with second ortho: " << fullSlidingFitResult3D.GetAxisDirection().GetDotProduct(fullSlidingFitResult3D.GetSecondFitResult().GetOrthoDirection()) << std::endl;

        std::cout << "RIGHT HANDED CHECK: " << fullSlidingFitResult3D.GetFirstFitResult().GetOrthoDirection().GetCrossProduct(fullSlidingFitResult3D.GetSecondFitResult().GetOrthoDirection()) 
                  << std::endl;

        const CartesianVector directionAxis(fullSlidingFitResult3D.GetAxisDirection());
        const CartesianVector &orthoAxis1(fullSlidingFitResult3D.GetFirstFitResult().GetOrthoDirection());
        const CartesianVector &orthoAxis2(fullSlidingFitResult3D.GetSecondFitResult().GetOrthoDirection());

        // Get angular asymmetry of hits (i don't think that there is any reason for this to be energy weighted)

        float meanTransverseAngle(0.f);

        for (const CaloHit *const pCaloHit : postShowerHitList3D)
        {
            const CartesianVector position(pCaloHit->GetPositionVector() - fullSlidingFitResult3D.GetAxisIntercept());
            const float angleToDirection(directionAxis.GetOpeningAngle(position));
            const float directionProjectionMagnitude(position.GetMagnitude() * std::cos(angleToDirection));
            const CartesianVector projectedPosition(position - (directionAxis * directionProjectionMagnitude));
            float angleToOrthoAxis1(orthoAxis1.GetOpeningAngle(projectedPosition));
            const float thisL(orthoAxis2.GetDotProduct(projectedPosition));

            if (thisL < 0.f)
                angleToOrthoAxis1 = M_PI + (M_PI - angleToOrthoAxis1);

            meanTransverseAngle += angleToOrthoAxis1;
        }

        // postShowerHitList3D.size() will not be zero as fit has succeeded
        meanTransverseAngle /= postShowerHitList3D.size();
        electronTreeVariables.m_postShowerStartMeanTransverseAngle = meanTransverseAngle;
        std::cout << electronTreeVariables.m_postShowerStartMeanTransverseAngle << std::endl;

        // Get transverse spread - will be larger for showers and not for tracks and will be affected by any long tracks before shower

        float meanRadialDistance(0.f);

        for (const CaloHit *const pCaloHit : postShowerHitList3D)
        {
            const CartesianVector position(pCaloHit->GetPositionVector() - fullSlidingFitResult3D.GetAxisIntercept());
            const float angleToDirection(directionAxis.GetOpeningAngle(position));
            const float directionProjectionMagnitude(position.GetMagnitude() * std::cos(angleToDirection));
            const CartesianVector projectedPosition(position - (directionAxis * directionProjectionMagnitude));
            float thisT(projectedPosition.GetMagnitude());
            const float thisL(orthoAxis2.GetDotProduct(projectedPosition));

            if (thisL < 0.f)
                thisT *= (-1.f);

            meanRadialDistance += thisT;
        }

        meanRadialDistance /= postShowerHitList3D.size();
        electronTreeVariables.m_postShowerStartMeanRadialDistance = meanRadialDistance;
        std::cout << "electronTreeVariables.m_postShowerStartMeanRadialDistance: " << electronTreeVariables.m_postShowerStartMeanRadialDistance << std::endl;

        float radialDistanceSigma(0.f);

        for (const CaloHit *const pCaloHit : postShowerHitList3D)
        {
            const CartesianVector position(pCaloHit->GetPositionVector() - fullSlidingFitResult3D.GetAxisIntercept());
            const float angleToDirection(directionAxis.GetOpeningAngle(position));
            const float directionProjectionMagnitude(position.GetMagnitude() * std::cos(angleToDirection));
            const CartesianVector projectedPosition(position - (directionAxis * directionProjectionMagnitude));
            float thisT(projectedPosition.GetMagnitude());
            const float thisL(orthoAxis2.GetDotProduct(projectedPosition));

            if (thisL < 0.f)
                thisT *= (-1.f);

            radialDistanceSigma += std::pow((thisT - meanRadialDistance), 2);
        }

        radialDistanceSigma = std::sqrt(radialDistanceSigma / postShowerHitList3D.size());
        electronTreeVariables.m_postShowerStartRadialDistanceSigma = radialDistanceSigma;
        std::cout << "electronTreeVariables.m_postShowerStartRadialDistanceSigma: " << electronTreeVariables.m_postShowerStartRadialDistanceSigma << std::endl;

        float energyWeightedMeanRadialDistance(0.f);

        for (const CaloHit *const pCaloHit : postShowerHitList3D)
        {
            const CartesianVector position(pCaloHit->GetPositionVector() - fullSlidingFitResult3D.GetAxisIntercept());
            const float angleToDirection(directionAxis.GetOpeningAngle(position));
            const float directionProjectionMagnitude(position.GetMagnitude() * std::cos(angleToDirection));
            const CartesianVector projectedPosition(position - (directionAxis * directionProjectionMagnitude));
            float thisT(projectedPosition.GetMagnitude());
            const float thisL(orthoAxis2.GetDotProduct(projectedPosition));

            if (thisL < 0.f)
                thisT *= (-1.f);

            energyWeightedMeanRadialDistance += thisT * pCaloHit->GetElectromagneticEnergy() * 1000.f;
        }

        energyWeightedMeanRadialDistance /= postShowerHitList3D.size();
        electronTreeVariables.m_postShowerStartEnergyWeightedMeanRadialDistance = energyWeightedMeanRadialDistance;
        std::cout << "electronTreeVariables.m_postShowerStartEnergyWeightedMeanRadialDistance: " << electronTreeVariables.m_postShowerStartEnergyWeightedMeanRadialDistance << std::endl;

        float energyWeightedRadialDistanceSigma(0.f);

        for (const CaloHit *const pCaloHit : postShowerHitList3D)
        {
            const CartesianVector position(pCaloHit->GetPositionVector() - fullSlidingFitResult3D.GetAxisIntercept());
            const float angleToDirection(directionAxis.GetOpeningAngle(position));
            const float directionProjectionMagnitude(position.GetMagnitude() * std::cos(angleToDirection));
            const CartesianVector projectedPosition(position - (directionAxis * directionProjectionMagnitude));
            float thisT(projectedPosition.GetMagnitude());
            const float thisL(orthoAxis2.GetDotProduct(projectedPosition));

            if (thisL < 0.f)
                thisT *= (-1.f);

            energyWeightedRadialDistanceSigma += std::pow(((thisT * pCaloHit->GetElectromagneticEnergy() * 1000.f) - energyWeightedMeanRadialDistance), 2);
        }

        energyWeightedRadialDistanceSigma = std::sqrt(energyWeightedRadialDistanceSigma / postShowerHitList3D.size());

        const float estimatedMoliereRadius(energyWeightedRadialDistanceSigma * 2.f);
        electronTreeVariables.m_postShowerStartEstimatedMoliereRadius = estimatedMoliereRadius;
        std::cout << "electronTreeVariables.m_postShowerStartEstimatedMoliereRadius: " << electronTreeVariables.m_postShowerStartEstimatedMoliereRadius << std::endl;

        // now look for any gaps (do this from the middle shower start)
        //const CartesianVector &minFitPosition(isDownstream ? fullSlidingFitResult3D.GetGlobalMinLayerPosition() : fullSlidingFitResult3D.GetGlobalMaxLayerPosition());
        const CartesianVector &fitShowerDirection3D(fullSlidingFitResult3D.GetGlobalMinLayerDirection() * (isDownstream ? 1.f : -1.f));

        FloatVector longitudinalProjections;

        for (const CaloHit * const pCaloHit : postShowerHitList3D)
        {
            const CartesianVector displacement(pCaloHit->GetPositionVector() - middleShowerStart3D);
            longitudinalProjections.push_back(fitShowerDirection3D.GetDotProduct(displacement));
        }

        std::sort(longitudinalProjections.begin(), longitudinalProjections.end());

        const float initialGapSize(longitudinalProjections[0]);
        electronTreeVariables.m_postShowerStartInitialGapSize = initialGapSize;
        std::cout << "electronTreeVariables.m_postShowerStartInitialGapSize: " << electronTreeVariables.m_postShowerStartInitialGapSize << std::endl;

        float maxGapSize(-std::numeric_limits<float>::max());

        for (unsigned int i = 1; i < longitudinalProjections.size(); ++i)
        {
            if ((i != 1) && (longitudinalProjections[i] > 10.f))
                continue;

            const float gapSize(longitudinalProjections[i] - longitudinalProjections[i - 1]);

            if (gapSize > maxGapSize)
                maxGapSize = gapSize;
        }

        electronTreeVariables.m_postShowerStartMaxGapSize = maxGapSize;
        std::cout << "electronTreeVariables.m_postShowerStartMaxGapSize: " << electronTreeVariables.m_postShowerStartMaxGapSize << std::endl;



    }
    catch (...)
    {
        std::cout << "one of the fits didnt work... :(" << std::endl;
        return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FillConnectionPathwayVariables(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU, 
    const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const CartesianVector &nuVertexPosition, CartesianVector &minShowerStart3D,
    CartesianVector &middleShowerStart3D, CartesianVector &maxShowerStart3D, LArConnectionPathwayHelper::ElectronTreeVariables &electronTreeVariables)
{
    electronTreeVariables.m_pathwayLength = (nuVertexPosition - middleShowerStart3D).GetMagnitude();

    electronTreeVariables.m_pathwayShowerStartDelta = (maxShowerStart3D - minShowerStart3D).GetMagnitude();

    electronTreeVariables.m_pathwayMaxScatteringAngle = LArConnectionPathwayHelper::GetLargest3DKink(pAlgorithm, protoShowerU, protoShowerV, protoShowerW, 
        nuVertexPosition, maxShowerStart3D);

    LArConnectionPathwayHelper::FillConnectionPathwayEnergyVariables(pAlgorithm, protoShowerU, protoShowerV, protoShowerW, middleShowerStart3D, electronTreeVariables);

    //////////////////////////////
    std::cout << "electronTreeVariables.m_pathwayLength: " << electronTreeVariables.m_pathwayLength << std::endl;
    std::cout << "electronTreeVariables.m_pathwayShowerStartDelta: " << electronTreeVariables.m_pathwayShowerStartDelta << std::endl;
    std::cout << "electronTreeVariables.m_pathwayMaxScatteringAngle: " << electronTreeVariables.m_pathwayMaxScatteringAngle << std::endl;
    std::cout << "electronTreeVariables.m_pathwayEnergyMeanU: " << electronTreeVariables.m_pathwayEnergyMeanU << std::endl;
    std::cout << "electronTreeVariables.m_pathwayEnergyMeanV: " << electronTreeVariables.m_pathwayEnergyMeanV << std::endl;
    std::cout << "electronTreeVariables.m_pathwayEnergyMeanW: " << electronTreeVariables.m_pathwayEnergyMeanW << std::endl;
    std::cout << "electronTreeVariables.m_pathwayEnergySigmaU: " << electronTreeVariables.m_pathwayEnergySigmaU << std::endl;
    std::cout << "electronTreeVariables.m_pathwayEnergySigmaV: " << electronTreeVariables.m_pathwayEnergySigmaV << std::endl;
    std::cout << "electronTreeVariables.m_pathwayEnergySigmaW: " << electronTreeVariables.m_pathwayEnergySigmaW << std::endl;
    std::cout << "electronTreeVariables.m_pathwayMinEnergyMean: " << electronTreeVariables.m_pathwayMinEnergyMean << std::endl;
    std::cout << "electronTreeVariables.m_pathwayMiddleEnergyMean: " << electronTreeVariables.m_pathwayMiddleEnergyMean << std::endl;
    std::cout << "electronTreeVariables.m_pathwayMaxEnergyMean: " << electronTreeVariables.m_pathwayMaxEnergyMean << std::endl;
    std::cout << "electronTreeVariables.m_pathwayMinEnergySigma: " << electronTreeVariables.m_pathwayMinEnergySigma << std::endl;
    std::cout << "electronTreeVariables.m_pathwayMiddleEnergySigma: " << electronTreeVariables.m_pathwayMiddleEnergySigma << std::endl;
    std::cout << "electronTreeVariables.m_pathwayMaxEnergySigma: " << electronTreeVariables.m_pathwayMaxEnergySigma << std::endl;
    /////////////////////////////

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArConnectionPathwayHelper::GetLargest3DKink(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU, 
    const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const CartesianVector &nuVertexPosition, CartesianVector &maxShowerStart3D)
{
    float maxOpeningAngle(-std::numeric_limits<float>::max());

    CartesianPointVector spinePositionsU, spinePositionsV, spinePositionsW;

    for (const CaloHit *const pCaloHit : protoShowerU.m_spineHitList)
        spinePositionsU.push_back(pCaloHit->GetPositionVector());

    for (const CaloHit *const pCaloHit : protoShowerV.m_spineHitList)
        spinePositionsV.push_back(pCaloHit->GetPositionVector());

    for (const CaloHit *const pCaloHit : protoShowerW.m_spineHitList)
        spinePositionsW.push_back(pCaloHit->GetPositionVector());

    const TwoDSlidingFitResult spineFitU(&spinePositionsU, 10, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));
    const TwoDSlidingFitResult spineFitV(&spinePositionsV, 10, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));
    const TwoDSlidingFitResult spineFitW(&spinePositionsW, 10, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));

    maxOpeningAngle = std::max(maxOpeningAngle, LArConnectionPathwayHelper::GetLargest3DKinkFromView(pAlgorithm, spineFitU, spineFitV, spineFitW, TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, maxShowerStart3D));
    maxOpeningAngle = std::max(maxOpeningAngle, LArConnectionPathwayHelper::GetLargest3DKinkFromView(pAlgorithm, spineFitV, spineFitW, spineFitU, TPC_VIEW_V, TPC_VIEW_W, TPC_VIEW_U, maxShowerStart3D));
    maxOpeningAngle = std::max(maxOpeningAngle, LArConnectionPathwayHelper::GetLargest3DKinkFromView(pAlgorithm, spineFitW, spineFitU, spineFitV, TPC_VIEW_W, TPC_VIEW_U, TPC_VIEW_V, maxShowerStart3D));

    return maxOpeningAngle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArConnectionPathwayHelper::GetLargest3DKinkFromView(pandora::Algorithm *const pAlgorithm, const TwoDSlidingFitResult &spineFit, 
    const TwoDSlidingFitResult &spineFit1, const TwoDSlidingFitResult &spineFit2, const HitType hitType, const HitType hitType1, 
    const HitType hitType2, const CartesianVector &maxShowerStart3D)
{
    const LayerFitResultMap &layerFitResultMap(spineFit.GetLayerFitResultMap());
    const int minLayer(layerFitResultMap.begin()->first), maxLayer(layerFitResultMap.rbegin()->first);

    const int nLayersHalfWindow(spineFit.GetLayerFitHalfWindow());
    const int nLayersSpanned(1 + maxLayer - minLayer);

    if (nLayersSpanned <= 2 * nLayersHalfWindow)
    {
        std::cout << "failed this window cut" << std::endl;
        return -10.f;
    }

    const CartesianVector projectedShowerStart(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), maxShowerStart3D, hitType));

    float showerStartL(0.f), showerStartT(0.f);
    spineFit.GetLocalPosition(projectedShowerStart, showerStartL, showerStartT);
    float maxCentralLayer(spineFit.GetLayer(showerStartL));

    float highestOpeningAngle(-10.f);
    CartesianVector kinkPosition(0.f, 0.f, 0.f);

    for (LayerFitResultMap::const_iterator iter = layerFitResultMap.begin(), iterEnd = layerFitResultMap.end(); iter != iterEnd; ++iter)
    {
        const int iLayer(iter->first);

        if (iLayer > maxCentralLayer)
            continue;

        const float rL(spineFit.GetL(iLayer));
        const float rL1(spineFit.GetL(iLayer - nLayersHalfWindow));
        const float rL2(spineFit.GetL(iLayer + nLayersHalfWindow));

        CartesianVector firstPosition(0.f, 0.f, 0.f), centralPosition(0.f, 0.f, 0.f), secondPosition(0.f, 0.f, 0.f);

        if ((STATUS_CODE_SUCCESS != spineFit.GetGlobalFitPosition(rL1, firstPosition)) ||
            (STATUS_CODE_SUCCESS != spineFit.GetGlobalFitPosition(rL, centralPosition)) ||
            (STATUS_CODE_SUCCESS != spineFit.GetGlobalFitPosition(rL2, secondPosition)))
       {
           continue;
       }

        float firstPositionX(firstPosition.GetX()), centralPositionX(centralPosition.GetX()), secondPositionX(secondPosition.GetX());

        CartesianVector firstPosition1(0.f, 0.f, 0.f), centralPosition1(0.f, 0.f, 0.f), secondPosition1(0.f, 0.f, 0.f);
        CartesianVector firstPosition2(0.f, 0.f, 0.f), centralPosition2(0.f, 0.f, 0.f), secondPosition2(0.f, 0.f, 0.f);

        if ((spineFit1.GetGlobalFitPositionAtX(firstPositionX, firstPosition1) != STATUS_CODE_SUCCESS) ||
            (spineFit1.GetGlobalFitPositionAtX(centralPositionX, centralPosition1) != STATUS_CODE_SUCCESS) ||
            (spineFit1.GetGlobalFitPositionAtX(secondPositionX, secondPosition1) != STATUS_CODE_SUCCESS) ||
            (spineFit2.GetGlobalFitPositionAtX(firstPositionX, firstPosition2) != STATUS_CODE_SUCCESS) ||
            (spineFit2.GetGlobalFitPositionAtX(centralPositionX, centralPosition2) != STATUS_CODE_SUCCESS) ||
            (spineFit2.GetGlobalFitPositionAtX(secondPositionX, secondPosition2) != STATUS_CODE_SUCCESS))
        {
            continue;
        }

        float metric(0.f);

        const CartesianVector &firstPositionU(hitType == TPC_VIEW_U ? firstPosition : hitType1 == TPC_VIEW_U ? firstPosition1 : firstPosition2);
        const CartesianVector &centralPositionU(hitType == TPC_VIEW_U ? centralPosition : hitType1 == TPC_VIEW_U ? centralPosition1 : centralPosition2);
        const CartesianVector &secondPositionU(hitType == TPC_VIEW_U ? secondPosition : hitType1 == TPC_VIEW_U ? secondPosition1 : secondPosition2);

        const CartesianVector &firstPositionV(hitType == TPC_VIEW_V ? firstPosition : hitType1 == TPC_VIEW_V ? firstPosition1 : firstPosition2);
        const CartesianVector &centralPositionV(hitType == TPC_VIEW_V ? centralPosition : hitType1 == TPC_VIEW_V ? centralPosition1 : centralPosition2);
        const CartesianVector &secondPositionV(hitType == TPC_VIEW_V ? secondPosition : hitType1 == TPC_VIEW_V ? secondPosition1 : secondPosition2);

        const CartesianVector &firstPositionW(hitType == TPC_VIEW_W ? firstPosition : hitType1 == TPC_VIEW_W ? firstPosition1 : firstPosition2);
        const CartesianVector &centralPositionW(hitType == TPC_VIEW_W ? centralPosition : hitType1 == TPC_VIEW_W ? centralPosition1 : centralPosition2);
        const CartesianVector &secondPositionW(hitType == TPC_VIEW_W ? secondPosition : hitType1 == TPC_VIEW_W ? secondPosition1 : secondPosition2);

        if (!LArConnectionPathwayHelper::AreShowerStartsConsistent(pAlgorithm, firstPositionU, firstPositionV, firstPositionW, 5.f, 2.f, metric))
            continue;

        if (!LArConnectionPathwayHelper::AreShowerStartsConsistent(pAlgorithm, centralPositionU, centralPositionV, centralPositionW, 5.f, 2.f, metric))
            continue;

        if (!LArConnectionPathwayHelper::AreShowerStartsConsistent(pAlgorithm, secondPositionU, secondPositionV, secondPositionW, 5.f, 2.f, metric))
            continue;

        float chi2(0.0);
        CartesianVector firstPosition3D(0.f, 0.f, 0.f), centralPosition3D(0.f, 0.f, 0.f), secondPosition3D(0.f, 0.f, 0.f);

        LArGeometryHelper::MergeThreePositions3D(pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, firstPositionU, firstPositionV, firstPositionW, firstPosition3D, chi2);
        LArGeometryHelper::MergeThreePositions3D(pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, centralPositionU, centralPositionV, centralPositionW, centralPosition3D, chi2);
        LArGeometryHelper::MergeThreePositions3D(pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, secondPositionU, secondPositionV, secondPositionW, secondPosition3D, chi2);

        const CartesianVector firstDirection3D((centralPosition3D - firstPosition3D).GetUnitVector());
        const CartesianVector secondDirection3D((secondPosition3D - centralPosition3D).GetUnitVector());

        const float openingAngle3D(secondDirection3D.GetOpeningAngle(firstDirection3D) * 180.f / M_PI);

        if (openingAngle3D > highestOpeningAngle)
        {
            highestOpeningAngle = openingAngle3D;
            kinkPosition = centralPosition3D;
        }
    }

    //PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &kinkPosition, "kinkPosition", RED, 2);
    //std::cout << "highest opening angle in view: " << highestOpeningAngle << std::endl;
    //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());

    return highestOpeningAngle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FillConnectionPathwayEnergyVariables(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU, 
    const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, CartesianVector &middleShowerStart3D, LArConnectionPathwayHelper::ElectronTreeVariables &electronTreeVariables)
{
    const CaloHitList connectionPathwayU(LArConnectionPathwayHelper::GetConnectionPathwayProjection(pAlgorithm, protoShowerU, middleShowerStart3D, TPC_VIEW_U));
    const CaloHitList connectionPathwayV(LArConnectionPathwayHelper::GetConnectionPathwayProjection(pAlgorithm, protoShowerV, middleShowerStart3D, TPC_VIEW_V));
    const CaloHitList connectionPathwayW(LArConnectionPathwayHelper::GetConnectionPathwayProjection(pAlgorithm, protoShowerW, middleShowerStart3D, TPC_VIEW_W));

    electronTreeVariables.m_pathwayEnergyMeanU = LArConnectionPathwayHelper::CharacteriseEnergyMean(connectionPathwayU);
    electronTreeVariables.m_pathwayEnergyMeanV = LArConnectionPathwayHelper::CharacteriseEnergyMean(connectionPathwayV);
    electronTreeVariables.m_pathwayEnergyMeanW = LArConnectionPathwayHelper::CharacteriseEnergyMean(connectionPathwayW);

    electronTreeVariables.m_pathwayEnergySigmaU = LArConnectionPathwayHelper::CharacteriseEnergySigma(connectionPathwayU, electronTreeVariables.m_pathwayEnergyMeanU);
    electronTreeVariables.m_pathwayEnergySigmaV = LArConnectionPathwayHelper::CharacteriseEnergySigma(connectionPathwayV, electronTreeVariables.m_pathwayEnergyMeanV);
    electronTreeVariables.m_pathwayEnergySigmaW = LArConnectionPathwayHelper::CharacteriseEnergySigma(connectionPathwayW, electronTreeVariables.m_pathwayEnergyMeanW);

    electronTreeVariables.m_pathwayMinEnergyMean = std::min(std::min(electronTreeVariables.m_pathwayEnergyMeanU, electronTreeVariables.m_pathwayEnergyMeanV), 
        electronTreeVariables.m_pathwayEnergyMeanW);

    electronTreeVariables.m_pathwayMaxEnergyMean = std::max(std::max(electronTreeVariables.m_pathwayEnergyMeanU, electronTreeVariables.m_pathwayEnergyMeanV), 
        electronTreeVariables.m_pathwayEnergyMeanW);

    for (const float mean : {electronTreeVariables.m_pathwayEnergyMeanU, electronTreeVariables.m_pathwayEnergyMeanV, electronTreeVariables.m_pathwayEnergyMeanW})
    {
        if ((std::fabs(mean - electronTreeVariables.m_pathwayMinEnergyMean) > std::numeric_limits<float>::epsilon()) &&
            (std::fabs(mean - electronTreeVariables.m_pathwayMaxEnergyMean) > std::numeric_limits<float>::epsilon()))
        {
            electronTreeVariables.m_pathwayMiddleEnergyMean = mean;
        }
    }

    electronTreeVariables.m_pathwayMinEnergySigma = std::min(std::min(electronTreeVariables.m_pathwayEnergySigmaU, electronTreeVariables.m_pathwayEnergySigmaV), 
        electronTreeVariables.m_pathwayEnergySigmaW);

    electronTreeVariables.m_pathwayMaxEnergySigma = std::max(std::max(electronTreeVariables.m_pathwayEnergySigmaU, electronTreeVariables.m_pathwayEnergySigmaV), 
        electronTreeVariables.m_pathwayEnergySigmaW);

    for (const float sigma : {electronTreeVariables.m_pathwayEnergySigmaU, electronTreeVariables.m_pathwayEnergySigmaV, electronTreeVariables.m_pathwayEnergySigmaW})
    {
        if ((std::fabs(sigma - electronTreeVariables.m_pathwayMinEnergySigma) > std::numeric_limits<float>::epsilon()) &&
            (std::fabs(sigma - electronTreeVariables.m_pathwayMaxEnergySigma) > std::numeric_limits<float>::epsilon()))
        {
            electronTreeVariables.m_pathwayMiddleEnergySigma = sigma;
        }
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CaloHitList LArConnectionPathwayHelper::GetConnectionPathwayProjection(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShower, 
    const CartesianVector &showerStart3D, const HitType hitType)
{
    const CartesianVector projectedShowerStart(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), showerStart3D, hitType));

    CartesianPointVector spinePositions;

    for (const CaloHit *const pCaloHit : protoShower.m_spineHitList)
        spinePositions.push_back(pCaloHit->GetPositionVector());

    const TwoDSlidingFitResult spineFit(&spinePositions, 20, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));

    CaloHitList projectedPathwayHitList;
    float lShowerStart(0.f), tShowerStart(0.f);
    spineFit.GetLocalPosition(projectedShowerStart, lShowerStart, tShowerStart);

    const bool isDownstream(protoShower.m_connectionPathway.m_startDirection.GetZ() > 0.f);

    for (const CaloHit *const pCaloHit : protoShower.m_spineHitList)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());

        float thisL(0.f), thisT(0.f);
        spineFit.GetLocalPosition(hitPosition, thisL, thisT);

        if (isDownstream && (thisL < lShowerStart))
            projectedPathwayHitList.push_back(pCaloHit);

        if (!isDownstream && (thisL > lShowerStart))
            projectedPathwayHitList.push_back(pCaloHit);
    }

    return projectedPathwayHitList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArConnectionPathwayHelper::CharacteriseEnergyMean(const CaloHitList &caloHitList)
{
    float totalEnergy(0.f);

    for (const CaloHit *const pCaloHit : caloHitList)
        totalEnergy += (pCaloHit->GetElectromagneticEnergy() * 1000.f);

    return (caloHitList.size() == 0 ? -10.f : totalEnergy / static_cast<float>(caloHitList.size()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArConnectionPathwayHelper::CharacteriseEnergySigma(const CaloHitList &caloHitList, const float energyMean)
{ 
    float energySigma(0.f);

    for (const CaloHit *const pCaloHit : caloHitList)
        energySigma += std::pow((pCaloHit->GetElectromagneticEnergy() * 1000.f) - energyMean, 2);

    return (caloHitList.size() == 0 ? -10.f : std::sqrt(energySigma / static_cast<float>(caloHitList.size())));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FindShowerStarts3D(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU, 
    const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const LArConnectionPathwayHelper::Consistency &consistency,
    const CartesianVector &nuVertexPosition, CartesianVector &minShowerStart3D, CartesianVector &middleShowerStart3D, CartesianVector &maxShowerStart3D)
{
    // ISOBEL YOU ARE SO LAZY
    float maxSeparation(5.f);

    std::cout << "consistency: " << consistency << std::endl;

    bool uFound(false), vFound(false), wFound(false);
    CartesianVector uShowerStart3D(0.f, 0.f, 0.f), vShowerStart3D(0.f, 0.f, 0.f), wShowerStart3D(0.f, 0.f, 0.f);

    if (consistency == LArConnectionPathwayHelper::Consistency::POSITION)
    {
        LArConnectionPathwayHelper::FindShowerVertexFromPosition(pAlgorithm, protoShowerU, protoShowerV, protoShowerW, uShowerStart3D);
        vShowerStart3D = uShowerStart3D;
        wShowerStart3D = uShowerStart3D;

        std::cout << "FOUND FROM POSITION" << std::endl;

        uFound = true; vFound = true; wFound = true;
    }
    else if (consistency == LArConnectionPathwayHelper::Consistency::DIRECTION)
    {
        std::cout << "IN DIRECTION VERTEX CHECK" << std::endl;

        if (LArConnectionPathwayHelper::FindShowerVertexFromDirection(pAlgorithm, protoShowerU, protoShowerV, protoShowerW, uShowerStart3D, vShowerStart3D, wShowerStart3D))
        {
            std::cout << "FOUND FROM DIRECTION" << std::endl;

            uFound = true; vFound = true; wFound = true;
        }
    }

    if (!uFound || (consistency == LArConnectionPathwayHelper::Consistency::X_PROJECTION))
    {
        if (LArConnectionPathwayHelper::FindShowerVertexFromXProjection(pAlgorithm, protoShowerU, protoShowerV, protoShowerW, maxSeparation, uShowerStart3D))
        {
            std::cout << "FOUND U FROM X PROJECTION" << std::endl;
            uFound = true;
        }

        if (LArConnectionPathwayHelper::FindShowerVertexFromXProjection(pAlgorithm, protoShowerV, protoShowerW, protoShowerU, maxSeparation, vShowerStart3D))
        {
            std::cout << "FOUND V FROM X PROJECTION" << std::endl;
            vFound = true;
        }

        if (LArConnectionPathwayHelper::FindShowerVertexFromXProjection(pAlgorithm, protoShowerW, protoShowerU, protoShowerV, maxSeparation, wShowerStart3D))
        {
            std::cout << "FOUND W FROM X PROJECTION" << std::endl;
            wFound = true;
        }
    }

    CartesianPointVector showerStarts3D;

    if (uFound)
        showerStarts3D.push_back(uShowerStart3D);

    if (vFound)
        showerStarts3D.push_back(vShowerStart3D);

    if (wFound)
        showerStarts3D.push_back(wShowerStart3D);

    std::sort(showerStarts3D.begin(), showerStarts3D.end(), LArConnectionPathwayHelper::SortByDistanceToPoint(nuVertexPosition));

    if (showerStarts3D.empty())
        return false;

    minShowerStart3D = showerStarts3D.front();
    middleShowerStart3D = (showerStarts3D.size() == 3) ? showerStarts3D[1] : showerStarts3D[0];
    maxShowerStart3D = showerStarts3D.back();

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------
// Determining Whether ProtoShowers Are Consistent
//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::AreShowerStartsConsistent(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU, 
    const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const float maxXSeparation, const float maxSeparation)
{
    float metric(0.0);

    const CartesianVector &showerStartU(protoShowerU.m_showerCore.m_startPosition);
    const CartesianVector &showerStartV(protoShowerV.m_showerCore.m_startPosition);
    const CartesianVector &showerStartW(protoShowerW.m_showerCore.m_startPosition);

    return LArConnectionPathwayHelper::AreShowerStartsConsistent(pAlgorithm, showerStartU, showerStartV, showerStartW, maxXSeparation, 
        maxSeparation, metric);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::AreShowerStartsConsistent(pandora::Algorithm *const pAlgorithm, const CartesianVector &showerStartU, 
    const CartesianVector &showerStartV, const CartesianVector &showerStartW, const float maxXSeparation, const float maxSeparation, float &metric)
{
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

    CartesianVector directionU(protoShowerU.m_connectionPathway.m_startDirection);
    CartesianVector directionV(protoShowerV.m_connectionPathway.m_startDirection);
    CartesianVector directionW(protoShowerW.m_connectionPathway.m_startDirection);

    return LArConnectionPathwayHelper::AreDirectionsConsistent(pAlgorithm, directionU, directionV, directionW, maxOpeningAngle, metric);
}

//------------------------------------------------------------------------------------------------------------------------------------------

// not by reference since we might change them...
bool LArConnectionPathwayHelper::AreDirectionsConsistent(pandora::Algorithm *const pAlgorithm, CartesianVector directionU, 
    CartesianVector directionV, CartesianVector directionW, const float maxOpeningAngle, float &metric)
{
    std::cout << "maxOpeningAngle: " << maxOpeningAngle << std::endl;

    const CartesianVector wireAxis(0.f, 0.f, 1.f);

    float wireDeviationU(directionU.GetOpeningAngle(wireAxis));
    wireDeviationU = std::min(wireDeviationU, static_cast<float>(M_PI - wireDeviationU));

    float wireDeviationV(directionV.GetOpeningAngle(wireAxis));
    wireDeviationV = std::min(wireDeviationV, static_cast<float>(M_PI - wireDeviationV));

    float wireDeviationW(directionW.GetOpeningAngle(wireAxis));
    wireDeviationW = std::min(wireDeviationW, static_cast<float>(M_PI - wireDeviationW));

    float radians((2.f * M_PI) / 180.f);
    bool isIsochronous((wireDeviationU < radians) && (wireDeviationV < radians) && (wireDeviationW < radians));

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
    const CartesianVector endU(nuVertexU + (directionU * 100.f));
    const CartesianVector endV(nuVertexV + (directionV * 100.f));
    const CartesianVector endW(nuVertexW + (directionW * 100.f));
    const CartesianVector projectionEndU(nuVertexU + (projectionU * 100.f));
    const CartesianVector projectionEndV(nuVertexV + (projectionV * 100.f));
    const CartesianVector projectionEndW(nuVertexW + (projectionW * 100.f));
    PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &nuVertexU, &projectionEndU, "PROJECTION U", RED, 2, 2);
    PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &nuVertexV, &projectionEndV, "PROJECTION V", RED, 2, 2);
    PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &nuVertexW, &projectionEndW, "PROJECTION W", RED, 2, 2);
    PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &nuVertexU, &endU, "DIRECTION U", BLACK, 2, 2);
    PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &nuVertexV, &endV, "DIRECTION V", BLACK, 2, 2);
    PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &nuVertexW, &endW, "DIRECTION W", BLACK, 2, 2);
    */
    std::cout << "angularDeviationU: " << openingAngleU << std::endl;
    std::cout << "angularDeviationV: " << openingAngleV << std::endl;
    std::cout << "angularDeviationW: " << openingAngleW << std::endl;
    /*
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    */
    /////////////////////////////////

    if ((openingAngleU > maxOpeningAngle) || (openingAngleV > maxOpeningAngle) || (openingAngleW > maxOpeningAngle))
    {
        std::cout << "OPENING ANGLE IS TOO BIG" << std::endl;
        return false;
    }

    metric = (openingAngleU + openingAngleV + openingAngleW) / 3.f;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------
// Finding a 3D Direction
//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FindDirection3D(pandora::Algorithm *const pAlgorithm, const bool isDownstream, const CartesianVector &startPositionU, 
    const CartesianVector &startPositionV, const CartesianVector &startPositionW, const CartesianVector &startPosition3D, const CartesianVector &directionU, 
    const CartesianVector &directionV, const CartesianVector &directionW, CartesianVector &direction3D)
{
    float metric(0.f);
    if (!LArConnectionPathwayHelper::AreDirectionsConsistent(pAlgorithm, directionU, directionV, directionW, 5.f, metric))
    {
        std::cout << "POST SHOWER START DIRECTIONS ARE NOT CONSISTENT" << std::endl;
        return false;
    }

    const CartesianVector xAxis(1.f, 0.f, 0.f);
    const float uDisplacement(std::fabs(1.f / directionU.GetCosOpeningAngle(xAxis)));
    const float vDisplacement(std::fabs(1.f / directionV.GetCosOpeningAngle(xAxis)));
    const float wDisplacement(std::fabs(1.f / directionW.GetCosOpeningAngle(xAxis)));

    const CartesianVector directionSeedU(startPositionU + (directionU * uDisplacement));
    const CartesianVector directionSeedV(startPositionV + (directionV * vDisplacement));
    const CartesianVector directionSeedW(startPositionW + (directionW * wDisplacement));

    if (!LArConnectionPathwayHelper::AreShowerStartsConsistent(pAlgorithm, directionSeedU, directionSeedV, directionSeedW, 5.f, 2.f, metric))
    {
        std::cout << "POST SHOWER DIRECTION SEEDS ARE NOT CONSISTENT" << std::endl;
        return false;
    }

    CartesianVector directionSeed3D(0.f, 0.f, 0.f);

    float chi2(0.0);
    LArGeometryHelper::MergeThreePositions3D(pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, directionSeedU, directionSeedV, directionSeedW, directionSeed3D, chi2);

    direction3D = (directionSeed3D - startPosition3D);

    if (isDownstream && (direction3D.GetZ() < 0.f))
        direction3D = direction3D * (-1.f);

    if (!isDownstream && (direction3D.GetZ() > 0.f))
        direction3D = direction3D * (-1.f);

    return true;
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
    CartesianVector uShowerStart3D(0.f, 0.f, 0.f), vShowerStart3D(0.f, 0.f, 0.f), wShowerStart3D(0.f, 0.f, 0.f);

    if(!LArConnectionPathwayHelper::FindShowerVertexFromDirection(pAlgorithm, protoShowerU, protoShowerV, protoShowerW, uShowerStart3D, 
        vShowerStart3D, wShowerStart3D))
    {
        return false;
    }

    float minSeparation(std::numeric_limits<float>::max());

    for (const CartesianVector showerStart : {uShowerStart3D, vShowerStart3D, wShowerStart3D})
    {
        const float separation((showerStart - nuVertexPosition).GetMagnitude());

        if (separation < minSeparation)
        {
            minSeparation = separation;
            showerStart3D = showerStart;
        }
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FindShowerVertexFromDirection(pandora::Algorithm *const pAlgorithm,
    const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, CartesianVector &uShowerStart3D, 
    CartesianVector &vShowerStart3D, CartesianVector &wShowerStart3D)
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

    std::cout << "transverseU: " << transverseU << std::endl;
    std::cout << "transverseV: " << transverseV << std::endl;
    std::cout << "transverseW: " << transverseW << std::endl;

    if ((transverseU > 1.f) || (transverseV > 1.f) || (transverseW > 1.f))
        return false;

    const CartesianVector xAxis(1.f, 0.f, 0.f);
    const float cosThetaU(std::fabs(peakDirectionU.GetCosOpeningAngle(xAxis)));
    const float cosThetaV(std::fabs(peakDirectionV.GetCosOpeningAngle(xAxis)));
    const float cosThetaW(std::fabs(peakDirectionW.GetCosOpeningAngle(xAxis)));

    const float x1(showerStartU.GetX());
    const CartesianVector showerStartV1(projectedNuVertexV + (peakDirectionV * (std::fabs(x1 - projectedNuVertexV.GetX()) / cosThetaV)));
    const CartesianVector showerStartW1(projectedNuVertexW + (peakDirectionW * (std::fabs(x1 - projectedNuVertexW.GetX()) / cosThetaW)));

    const float x2(showerStartV.GetX());
    const CartesianVector showerStartU2(projectedNuVertexU + (peakDirectionU * (std::fabs(x2 - projectedNuVertexU.GetX()) / cosThetaU)));
    const CartesianVector showerStartW2(projectedNuVertexW + (peakDirectionW * (std::fabs(x2 - projectedNuVertexW.GetX()) / cosThetaW)));

    const float x3(showerStartW.GetX());
    const CartesianVector showerStartU3(projectedNuVertexU + (peakDirectionU * (std::fabs(x3 - projectedNuVertexU.GetX()) / cosThetaU)));
    const CartesianVector showerStartV3(projectedNuVertexV + (peakDirectionV * (std::fabs(x3 - projectedNuVertexV.GetX()) / cosThetaV)));

    float chi2(0.0);

    LArGeometryHelper::MergeThreePositions3D(pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, showerStartU, showerStartV1, showerStartW1, uShowerStart3D, chi2);
    LArGeometryHelper::MergeThreePositions3D(pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, showerStartU2, showerStartV, showerStartW2, vShowerStart3D, chi2);
    LArGeometryHelper::MergeThreePositions3D(pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, showerStartU3, showerStartV3, showerStartW, wShowerStart3D, chi2);

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FindShowerVertexFromXProjection(pandora::Algorithm *const pAlgorithm, const CartesianVector nuVertexPosition, 
    const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const float maxSeparation, CartesianVector &showerStart3D)
{
    CartesianVector showerStartU(0.f, 0.f, 0.f), showerStartV(0.f, 0.f, 0.f), showerStartW(0.f, 0.f, 0.f);

    const bool uFound(LArConnectionPathwayHelper::FindShowerVertexFromXProjection(pAlgorithm, protoShowerU, protoShowerV, protoShowerW, maxSeparation, showerStartU));
    const bool vFound(LArConnectionPathwayHelper::FindShowerVertexFromXProjection(pAlgorithm, protoShowerV, protoShowerW, protoShowerU, maxSeparation, showerStartV));
    const bool wFound(LArConnectionPathwayHelper::FindShowerVertexFromXProjection(pAlgorithm, protoShowerW, protoShowerU, protoShowerV, maxSeparation, showerStartW));

    bool found(false);
    float minSeparation(std::numeric_limits<float>::max());

    if (uFound)
    {
        const float separation((showerStartU - nuVertexPosition).GetMagnitude());

        if (separation < minSeparation)
        {
            found = true;
            minSeparation = separation;
            showerStart3D = showerStartU;
        }
    }

    if (vFound)
    {
        const float separation((showerStartV - nuVertexPosition).GetMagnitude());

        if (separation < minSeparation)
        {
            found = true;
            minSeparation = separation;
            showerStart3D = showerStartV;
        }
    }

    if (wFound)
    {
        const float separation((showerStartW - nuVertexPosition).GetMagnitude());

        if (separation < minSeparation)
        {
            found = true;
            minSeparation = separation;
            showerStart3D = showerStartW;
        }
    }

    return found;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::FindShowerVertexFromXProjection(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShower, 
    const ProtoShower &protoShower1, const ProtoShower &protoShower2, const float maxSeparation, CartesianVector &showerStart3D)
{
    const CartesianVector showerStart(protoShower.m_showerCore.m_startPosition);
    CartesianVector showerStart1(0.f, 0.f, 0.f), showerStart2(0.f, 0.f, 0.f);

    const HitType hitType(protoShower.m_connectionPathway.m_pathwayHitList.front()->GetHitType());
    const HitType hitType1(protoShower1.m_connectionPathway.m_pathwayHitList.front()->GetHitType());
    const HitType hitType2(protoShower2.m_connectionPathway.m_pathwayHitList.front()->GetHitType());

    bool found1(false);
    float lowestL1(std::numeric_limits<float>::max());

    for (const CaloHit *const pCaloHit : protoShower1.m_spineHitList)
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
        return false;

    bool found2(false);
    float lowestL2(std::numeric_limits<float>::max());

    for (const CaloHit *const pCaloHit : protoShower2.m_spineHitList)
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
        return false;

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
        return false;
    }

    LArGeometryHelper::MergeThreePositions3D(pAlgorithm->GetPandora(), hitType, hitType1, hitType2, showerStart, showerStart1, showerStart2, showerStart3D, chi2);

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------
// Sorting Functions
//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::SortByDistanceToPoint::operator()(const CartesianVector &lhs, const CartesianVector &rhs)
{
    return (m_referencePoint.GetDistanceSquared(lhs) < m_referencePoint.GetDistanceSquared(rhs));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArConnectionPathwayHelper::SortByDistanceToPoint::operator()(const CaloHit *const lhs, const CaloHit *const rhs)
{
    return (m_referencePoint.GetDistanceSquared(lhs->GetPositionVector()) < m_referencePoint.GetDistanceSquared(rhs->GetPositionVector()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content
