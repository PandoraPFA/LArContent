/**
 *  @file   larpandoracontent/LArShowerRefinement/ElectronStartRefinementTool.cc
 *
 *  @brief  Implementation of the gamma start refinement tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/AlgorithmTool.h"

#include "PandoraMonitoringApi.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingShowerFitResult.h"

#include "larpandoracontent/LArShowerRefinement/ElectronStartRefinementTool.h"

using namespace pandora;

namespace lar_content
{

ElectronStartRefinementTool::ElectronStartRefinementTool() : 
    m_showerSlidingFitWindow(20),
    m_maxCoincideneTransverseSeparation(2.f),
    m_minSpinePurity(0.7f),
    m_maxAngularDeviation(5.f),
    m_maxXSeparation(5.f),
    m_extendMode(true),
    m_moveVertexMode(true)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

ElectronStartRefinementTool::~ElectronStartRefinementTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ElectronStartRefinementTool::Run(ShowerStartRefinementAlgorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo, 
    const CartesianVector &nuVertexPosition)
{
    // Only apply electron refinement algorithm to showers
    if (!LArPfoHelper::IsShower(pShowerPfo))
        return false;

    // Only consider significant showers
    CaloHitList caloHits3D;
    LArPfoHelper::GetCaloHits(pShowerPfo, TPC_3D, caloHits3D);

    if (caloHits3D.size() < 50)
        return false;

    if (!pAlgorithm->IsElectron(pShowerPfo))
        return false;

    CaloHitList usedHitListU, usedHitListV, usedHitListW;
    ElectronProtoShowerVector protoShowerVectorU, protoShowerVectorV, protoShowerVectorW;
    this->BuildProtoShowers(pAlgorithm, pShowerPfo, nuVertexPosition, TPC_VIEW_U, protoShowerVectorU, usedHitListU);
    this->BuildProtoShowers(pAlgorithm, pShowerPfo, nuVertexPosition, TPC_VIEW_V, protoShowerVectorV, usedHitListV);
    this->BuildProtoShowers(pAlgorithm, pShowerPfo, nuVertexPosition, TPC_VIEW_W, protoShowerVectorW, usedHitListW);

    // here find 3D shower start? if we could then... go straight to extend shower?
    if (!this->IsElectronPathway(pAlgorithm, protoShowerVectorU, protoShowerVectorV, protoShowerVectorW))
    {
        std::cout << "NOT ELECTRON PATHWAY" << std::endl;
        return false;
    }

    // Can maybe do a 'does 3D projection match shower direction?'
    if (!this->ArePathwaysConsistent(pAlgorithm, protoShowerVectorU.front(), protoShowerVectorV.front(), protoShowerVectorW.front()))
    {
        std::cout << "PATHWAYS ARE NOT CONSISTENT" << std::endl;
        return false;
    }

    if (m_extendMode)
    {
        if (this->IsShowerExtendable(pAlgorithm, protoShowerVectorU, protoShowerVectorV, protoShowerVectorW))
        {
            CartesianPointVector peakDirectionsU, peakDirectionsV, peakDirectionsW;
            this->BuildHelperProtoShowers(pAlgorithm, pShowerPfo, nuVertexPosition, TPC_VIEW_U, peakDirectionsU, usedHitListU);
            this->BuildHelperProtoShowers(pAlgorithm, pShowerPfo, nuVertexPosition, TPC_VIEW_V, peakDirectionsV, usedHitListV);
            this->BuildHelperProtoShowers(pAlgorithm, pShowerPfo, nuVertexPosition, TPC_VIEW_W, peakDirectionsW, usedHitListW);

            this->RefineHitsToAdd(pAlgorithm, protoShowerVectorU.front(), nuVertexPosition, TPC_VIEW_U, peakDirectionsU);
            this->RefineHitsToAdd(pAlgorithm, protoShowerVectorV.front(), nuVertexPosition, TPC_VIEW_V, peakDirectionsV);
            this->RefineHitsToAdd(pAlgorithm, protoShowerVectorW.front(), nuVertexPosition, TPC_VIEW_W, peakDirectionsW);

            this->ExtendShower(pAlgorithm, pShowerPfo, protoShowerVectorU, protoShowerVectorV, protoShowerVectorW);
        }
    }

    // set vertex
    if (m_moveVertexMode)
        pAlgorithm->SetElectronMetadata(nuVertexPosition, pShowerPfo);

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronStartRefinementTool::BuildProtoShowers(ShowerStartRefinementAlgorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo,
    const CartesianVector &nuVertexPosition, const HitType hitType, ElectronProtoShowerVector &protoShowerVector, CaloHitList &usedHitList)
{
    /////////////////////////////////
    //std::cout << "Investigating " << (hitType == TPC_VIEW_U ? "U" : hitType == TPC_VIEW_V ? "V" : "W") << " view..." << std::endl;
    /////////////////////////////////

    CaloHitList viewShowerHitList;
    const CaloHitList viewAllHitList(pAlgorithm->GetAllHitsOfType(hitType));
    LArPfoHelper::GetCaloHits(pShowerPfo, hitType, viewShowerHitList);

    const CartesianVector projectedNuVertexPosition(LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertexPosition, hitType));

    CaloHitList viewROIHits;
    this->CollectHitsWithinROI(pAlgorithm, viewShowerHitList, viewAllHitList, projectedNuVertexPosition, viewROIHits);

    if (viewROIHits.empty())
        return;

    // Fill angular decomposition map
    AngularDecompositionMap angularDecompositionMap;
    this->FillAngularDecompositionMap(viewROIHits, projectedNuVertexPosition, angularDecompositionMap);

    if (angularDecompositionMap.empty())
        return;

    // Smooth angular decomposition map
    this->SmoothAngularDecompositionMap(angularDecompositionMap);

    // Find peak directions
    IntVector angularPeakVector;
    this->ObtainPeakVector(angularDecompositionMap, angularPeakVector);

    // Now, for future use, find the shower's vertex <- if layer positions span the vertex then put vertex as nu vertex?? <- might mean we are able to shift the vertex?
    CartesianVector showerVertexPosition(0.f, 0.f, 0.f);
    try
    {
        showerVertexPosition = this->GetShowerVertex(pShowerPfo, hitType, projectedNuVertexPosition);
    }
    catch(...)
    {
        return;
    }

    // Investigate peaks
    IntVector investigatedPeaks;

    for (unsigned int i = 0; i < angularPeakVector.size(); ++i)
    {
        int bestTheta0XZBin(0);

        if (!this->FindBestAngularPeak(angularDecompositionMap, angularPeakVector, investigatedPeaks, bestTheta0XZBin))
            break;

        investigatedPeaks.push_back(bestTheta0XZBin);

        const float theta0XZ(bestTheta0XZBin * m_theta0XZBinSize);
        const CartesianVector peakDirection(std::cos(theta0XZ), 0.f, std::sin(theta0XZ));
        const bool isEndDownstream(peakDirection.GetZ() > 0.f);

        ///////////////////////////////////////////
        //const CartesianVector end(projectedNuVertexPosition + (peakDirection * 10.f));
        //PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &projectedNuVertexPosition, &end, "direction", BLACK, 2,2);
        ///////////////////////////////////////////

        // Find the shower spine
        CaloHitList showerSpineHitList;
        try
        {
            this->FindShowerSpine(pAlgorithm, viewAllHitList, projectedNuVertexPosition, peakDirection, usedHitList, showerSpineHitList);
        }
        catch (...)
        {
            ///////////////////////////////////////////  
            //std::cout << "SHOWER SPINE FIT FAILED" << std::endl;
            //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
            ///////////////////////////////////////////  
            continue;
        }

        // Demand that spine is significant, be lenient here as some have small stubs and a gap
        if (showerSpineHitList.size() < 7) 
        {
            ///////////////////////////////////////////  
            //std::cout << "SHOWER NOT SIGNIFICANT" << std::endl;
            //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
            ///////////////////////////////////////////  
            continue;
        }

        // Does spine coincide with shower pfo
        if (!this->IsSpineCoincident(pAlgorithm, projectedNuVertexPosition, peakDirection, viewShowerHitList, showerVertexPosition, showerSpineHitList))
        {
            ///////////////////////////////////////////  
            //std::cout << "not coincident" << std::endl;
            //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
            ///////////////////////////////////////////  
            continue;
        }

        // Now look for the shower start position, first obtain the longitudinal position of spine hits
        LongitudinalPositionMap longitudinalPositionMap;
        try
        {
            this->ObtainLongitudinalDecomposition(pAlgorithm, showerSpineHitList, longitudinalPositionMap);
        }
        catch(...)
        {
            /////////////////////////////////
            //std::cout << "the super fine fit failed probably a gap in the initial hits??" << std::endl;
            //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
            /////////////////////////////////
            continue;
        }

        // Obtain spine energy profile
        EnergySpectrumMap energySpectrumMap;
        this->GetEnergyDistribution(pAlgorithm, showerSpineHitList, longitudinalPositionMap, energySpectrumMap);

        // Postulate combining the two
        CaloHitList mergedHitList(viewShowerHitList);

        for (const CaloHit *const pCaloHit : showerSpineHitList)
        {
            if (std::find(viewShowerHitList.begin(), viewShowerHitList.end(), pCaloHit) == viewShowerHitList.end())
                mergedHitList.push_back(pCaloHit);
        }

        // I don't like the ProtoShowers being created here (do this so that i can push the vertex forwards if i am sure we have an electron!)
        ProtoShowerVector protoShowerVectorTemp;
        CartesianVector showerStartPosition(0.f, 0.f, 0.f);

        this->FindShowerStart(pAlgorithm, projectedNuVertexPosition, peakDirection, longitudinalPositionMap, energySpectrumMap, showerSpineHitList, 
             showerStartPosition, mergedHitList, isEndDownstream, protoShowerVectorTemp, false); 

        // add in the shower core bit (TO DO -  so then i can check??)
        CaloHitList hitsToAdd;
        const float showerVertexL(std::max(peakDirection.GetDotProduct(showerStartPosition - projectedNuVertexPosition), 
            peakDirection.GetDotProduct(showerVertexPosition - projectedNuVertexPosition)));

        for (const CaloHit *const pCaloHit : showerSpineHitList)
        {
            if (std::find(viewShowerHitList.begin(), viewShowerHitList.end(), pCaloHit) != viewShowerHitList.end())
                continue;

            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            const float showerL(peakDirection.GetDotProduct(hitPosition - projectedNuVertexPosition));

            if ((showerL > 0.f) && (showerL < showerVertexL))
                hitsToAdd.push_back(pCaloHit);
        }

        usedHitList.insert(usedHitList.begin(), showerSpineHitList.begin(), showerSpineHitList.end());

        protoShowerVector.push_back(ElectronProtoShower(protoShowerVectorTemp.front().m_showerCore, protoShowerVectorTemp.front().m_connectionPathway, false, hitsToAdd));

        /////////////////////////////////
        /*  
        PfoList visualize;
        visualize.push_back(pShowerPfo);
        PandoraMonitoringApi::VisualizeParticleFlowObjects(pAlgorithm->GetPandora(), &visualize, "ShowerPfo", BLUE);
        PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &showerSpineHitList, "Shower Spine", BLACK);
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerStartPosition, "Shower start position", BLACK, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &projectedNuVertexPosition, "Nu vertex position", GREEN, 2);

        for (const CaloHit *const pHitToAdd : hitsToAdd)
        {
            const CartesianVector &hitPosition(pHitToAdd->GetPositionVector());
            PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "HIT TO ADD", RED, 2);
        }
        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
        */
        /////////////////////////////////

        break;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronStartRefinementTool::BuildHelperProtoShowers(ShowerStartRefinementAlgorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo,
    const CartesianVector &nuVertexPosition, const HitType tpcView, CartesianPointVector &significantPeakDirections, CaloHitList &usedHitList)
{
    CaloHitList showerHitList;
    LArPfoHelper::GetCaloHits(pShowerPfo, tpcView, showerHitList);

    CaloHitList viewHitList;
    const CaloHitList allHitList(pAlgorithm->GetAllHitsOfType(tpcView));

    for (const CaloHit *const pCaloHit : allHitList)
    {
        if (std::find(showerHitList.begin(), showerHitList.end(), pCaloHit) != showerHitList.end())
            continue;

        viewHitList.push_back(pCaloHit);
    }

    const CartesianVector projectedNuVertexPosition(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), nuVertexPosition, tpcView));

    AngularDecompositionMap angularDecompositionMap;
    this->FillAngularDecompositionMap(viewHitList, projectedNuVertexPosition, angularDecompositionMap);

    if (angularDecompositionMap.empty())
        return;

    this->SmoothAngularDecompositionMap(angularDecompositionMap);

    // Obtain distribution peaks
    IntVector angularPeakVector;
    this->ObtainPeakVector(angularDecompositionMap, angularPeakVector);

    // Investigate peaks, searching for pathways (peaks indexed by theta0XZ bin lower edge)
    IntVector investigatedPeaks;

    for (unsigned int i = 0; i < angularPeakVector.size(); ++i)
    {
        int bestTheta0XZBin(0);

        if (!this->FindBestAngularPeak(angularDecompositionMap, angularPeakVector, investigatedPeaks, bestTheta0XZBin))
            break;

        investigatedPeaks.push_back(bestTheta0XZBin);

        const float theta0XZ(bestTheta0XZBin * m_theta0XZBinSize);
        const CartesianVector peakDirection(std::cos(theta0XZ), 0.f, std::sin(theta0XZ));

        CaloHitList showerSpineHitList;
        try
        {
            this->FindShowerSpine(pAlgorithm, viewHitList, projectedNuVertexPosition, peakDirection, usedHitList, showerSpineHitList);
        }
        catch (...)
        {
            continue;
        }

        // Demand that spine is significant
        if (showerSpineHitList.size() < 7)
            continue;

        usedHitList.insert(usedHitList.begin(), showerSpineHitList.begin(), showerSpineHitList.end());

        significantPeakDirections.push_back(peakDirection);

        ///////////////////////////////////
        /*
        const CartesianVector end(projectedNuVertexPosition + (peakDirection * 10.f));
        PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &projectedNuVertexPosition, &end, "peak direction", VIOLET, 2,2);

        for (const CaloHit *const pCaloHit : showerSpineHitList)
        {
            const CartesianVector hitPosition(pCaloHit->GetPositionVector());
            PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "shower spine", BLACK, 2);
        }
        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
        */
        ///////////////////////////////////
    }

    ///////////////////////////////////
    //std::cout << "Peak directions..." << std::endl;
    //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    ///////////////////////////////////
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronStartRefinementTool::RefineHitsToAdd(ShowerStartRefinementAlgorithm *const pAlgorithm, ElectronProtoShower &protoShower,
    const CartesianVector &nuVertexPosition, const HitType hitType, const CartesianPointVector &significantPeakDirections)
{
    const CartesianVector projectedNuVertexPosition(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), nuVertexPosition, hitType));
    const CartesianVector &peakDirection(protoShower.m_connectionPathway.m_startDirection);

    CaloHitList refinedHitList;

    for (const CaloHit *const pHitToAdd : protoShower.m_hitsToAdd)
    {
        bool found(false);

        const CartesianVector &hitPosition(pHitToAdd->GetPositionVector());
        const CartesianVector displacement(hitPosition - projectedNuVertexPosition);
        const float thisT((peakDirection.GetCrossProduct(displacement)).GetMagnitudeSquared());

        for (const CartesianVector &significantPeakDirection: significantPeakDirections)
        {
            const float otherT((significantPeakDirection.GetCrossProduct(displacement)).GetMagnitudeSquared());

            if ((otherT < thisT) || (otherT < 0.5f))
            {
                found = true;
                break;
            }
        }

        if (!found)
            refinedHitList.push_back(pHitToAdd);
    }

    //////////////////////////////
    /*
    std::cout << "BEFORE HIT REFINEMENT" << std::endl;
    for (const CaloHit *const pCaloHit :  protoShower.m_hitsToAdd)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "before hit refinement", GREEN, 2);
    }
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    */
    //////////////////////////////

    protoShower.m_hitsToAdd = this->FindContinuousPath(refinedHitList, projectedNuVertexPosition);

    //////////////////////////////
    /*
    std::cout << "AFTER HIT REFINEMENT" << std::endl;
    for (const CaloHit *const pCaloHit :  protoShower.m_hitsToAdd)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "after hit refinement", RED, 2);
    }
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    */
    //////////////////////////////
}

//------------------------------------------------------------------------------------------------------------------------------------------

CaloHitList ElectronStartRefinementTool::FindContinuousPath(const CaloHitList &refinedHitList, const CartesianVector &projectedNuVertexPosition)
{
    CaloHitList continuousHitList;

    CaloHitVector refinedHitVector(refinedHitList.begin(), refinedHitList.end());
    std::sort(refinedHitVector.begin(), refinedHitVector.end(), SortByDistanceToPoint(projectedNuVertexPosition));

    unsigned int startIndex(refinedHitVector.size());

    for (unsigned int i = 0; i < refinedHitVector.size(); ++i)
    {
        CaloHitList connectedHitList;
        connectedHitList.push_back(refinedHitVector[i]);

        bool found(true);

        while(found)
        {
            found = false;

            for (unsigned int j = (i + 1); j < refinedHitVector.size(); ++j)
            {
                const CaloHit *const pCaloHit(refinedHitVector[j]);

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

        if ((connectedHitList.size() >= 2) || (connectedHitList.size() == refinedHitVector.size()))
        {
            startIndex = i;
            break;
        }
    }

    for (unsigned int i = startIndex; i < refinedHitVector.size(); ++i)
        continuousHitList.push_back(refinedHitVector[i]);

    return continuousHitList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronStartRefinementTool::CollectHitsWithinROI(ShowerStartRefinementAlgorithm *const pAlgorithm, const CaloHitList &showerHitList, 
    const CaloHitList &allHitList, const CartesianVector &projectedNuVertexPosition, CaloHitList &collectedHits)
{
    float lowestTheta(std::numeric_limits<float>::max()), highestTheta((-1.f) * std::numeric_limits<float>::max());

    this->GetAngularExtrema(pAlgorithm, showerHitList, projectedNuVertexPosition, lowestTheta, highestTheta);
    this->CollectHitsWithinExtrema(pAlgorithm, allHitList, projectedNuVertexPosition, lowestTheta, highestTheta, collectedHits);

    ///////////////////////////////
    /*
    for (const CaloHit *const pCaloHit : collectedHits)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "ROI HIT", SPRING, 2);
    }
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    */
    ///////////////////////////////
}

//------------------------------------------------------------------------------------------------------------------------------------------
// i dont know if this is better off as a simple x thing? i.e. between nu vertex and shower edge?
void ElectronStartRefinementTool::GetAngularExtrema(ShowerStartRefinementAlgorithm *const pAlgorithm, const CaloHitList &viewHitList, 
    const CartesianVector &projectedNuVertexPosition, float &lowestTheta, float &highestTheta)
{
    lowestTheta = std::numeric_limits<float>::max();
    highestTheta = (-1.f) * std::numeric_limits<float>::max();

    const CartesianVector xAxis(1.f, 0.f, 0.f);

    for (const CaloHit *const pCaloHit : viewHitList)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        const CartesianVector displacementVector(hitPosition - projectedNuVertexPosition);

        float theta0XZ(xAxis.GetOpeningAngle(displacementVector));

        if (displacementVector.GetZ() < 0.f)
            theta0XZ += (3.14 - theta0XZ);

        if (theta0XZ < lowestTheta)
            lowestTheta = theta0XZ;

        if (theta0XZ > highestTheta)
            highestTheta = theta0XZ;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronStartRefinementTool::CollectHitsWithinExtrema(ShowerStartRefinementAlgorithm *const pAlgorithm, const CaloHitList &viewHitList, 
    const CartesianVector &projectedNuVertexPosition, const float lowestTheta, const float highestTheta, CaloHitList &collectedHits)
{
    const CartesianVector xAxis(1.f, 0.f, 0.f);

    for (const CaloHit *const pCaloHit : viewHitList)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        const CartesianVector displacementVector(hitPosition - projectedNuVertexPosition);

        float theta0XZ(xAxis.GetOpeningAngle(displacementVector));

        if (displacementVector.GetZ() < 0.f)
            theta0XZ += (3.14 - theta0XZ);

        if ((theta0XZ < lowestTheta) || (theta0XZ > highestTheta))
            continue;

        collectedHits.push_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector ElectronStartRefinementTool::GetShowerVertex(const ParticleFlowObject *const pShowerPfo, const HitType hitType, 
    const CartesianVector &projectedNuVertexPosition)
{
    const TwoDSlidingShowerFitResult twoDShowerSlidingFit(this->FitShower(pShowerPfo, hitType));

    const CartesianVector &minLayerPosition(twoDShowerSlidingFit.GetShowerFitResult().GetGlobalMinLayerPosition());
    const CartesianVector &maxLayerPosition(twoDShowerSlidingFit.GetShowerFitResult().GetGlobalMaxLayerPosition());

    if ((projectedNuVertexPosition.GetZ() > minLayerPosition.GetZ()) && (projectedNuVertexPosition.GetZ() < maxLayerPosition.GetZ()))
    {
        //std::cout << "pfo spans nu vertex" << std::endl;
        return projectedNuVertexPosition;
    }

    const float minSeparation((projectedNuVertexPosition - minLayerPosition).GetMagnitudeSquared());
    const float maxSeparation((projectedNuVertexPosition - maxLayerPosition).GetMagnitudeSquared());
    const CartesianVector &showerVertexPosition(minSeparation < maxSeparation ? minLayerPosition : maxLayerPosition);

    return showerVertexPosition;
}

//------------------------------------------------------------------------------------------------------------------------------------------

TwoDSlidingShowerFitResult ElectronStartRefinementTool::FitShower(const ParticleFlowObject *const pShowerPfo, const HitType hitType)
{
    ClusterList viewCusterList;
    LArPfoHelper::GetClusters(pShowerPfo, hitType, viewCusterList);

    const TwoDSlidingShowerFitResult twoDShowerSlidingFit(viewCusterList.front(), m_showerSlidingFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));

    return twoDShowerSlidingFit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ElectronStartRefinementTool::IsSpineCoincident(ShowerStartRefinementAlgorithm *const pAlgorithm, const CartesianVector &projectedNuVertexPosition, 
    const CartesianVector &peakDirection, const CaloHitList &viewShowerHitList, const CartesianVector &showerVertexPosition, const CaloHitList &showerSpineHitList)
{
    CaloHitList postShowerVertexSpineHits;
    const float showerDistanceFromNuSquared((showerVertexPosition - projectedNuVertexPosition).GetMagnitudeSquared());

    for (const CaloHit * const pSpineHit : showerSpineHitList)
    {
        const CartesianVector &hitPosition(pSpineHit->GetPositionVector());
        const float separationSquared((hitPosition - projectedNuVertexPosition).GetMagnitudeSquared());

        if (separationSquared > showerDistanceFromNuSquared)
            postShowerVertexSpineHits.push_back(pSpineHit);
    }

    /////////////////////////////////
    /*
    CartesianVector end(projectedNuVertexPosition + (peakDirection * 100.0));
    PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &projectedNuVertexPosition, &end, "DIRECTION", PINK, 2, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerVertexPosition, "SPINE VERTEX POSITION", GREEN, 2);
    */
    /////////////////////////////////

    // Check whether shower is an extension of found pathway
    const CartesianVector displacement(showerVertexPosition - projectedNuVertexPosition);

    const float longitudinalSeparation(peakDirection.GetDotProduct(displacement));

    if (longitudinalSeparation < (-1.f))
        return false;

    const float transverseSeparation((peakDirection.GetCrossProduct(displacement)).GetMagnitude());

    if (transverseSeparation > m_maxCoincideneTransverseSeparation)
    {
        /////////////////////////////////
        //std::cout << "transverse separation issue" << std::endl;
        //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
        /////////////////////////////////
        return false;
    }

    // Check whether shower spine is pure
    const CaloHitList sharedHitList(LArMCParticleHelper::GetSharedHits(postShowerVertexSpineHits, viewShowerHitList));
    const float spinePurity(postShowerVertexSpineHits.size() == 0 ? 1.f : static_cast<float>(sharedHitList.size()) / static_cast<float>(postShowerVertexSpineHits.size()));

    if (spinePurity < m_minSpinePurity)
    {
        /////////////////////////////////
        //std::cout <<  "spine purity issue" << std::endl;
        //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
        /////////////////////////////////
        return false;
    }

    /////////////////////////////////
    //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    /////////////////////////////////

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------
// I can add something here to make sure the shower actually has a shower bit once i sort out the shower core stuff?

bool ElectronStartRefinementTool::IsElectronPathway(ShowerStartRefinementAlgorithm *const pAlgorithm, const ElectronProtoShowerVector &protoShowerVectorU, 
    const ElectronProtoShowerVector &protoShowerVectorV, const ElectronProtoShowerVector &protoShowerVectorW)
{
    // Return if have not found a connection pathway in each view
    if ((protoShowerVectorU.size() != 1) || (protoShowerVectorV.size() != 1) || (protoShowerVectorW.size() != 1))
        return false;

    /*
    const ElectronProtoShower &protoShowerU(protoShowerVectorU.front());
    const ElectronProtoShower &protoShowerV(protoShowerVectorV.front());
    const ElectronProtoShower &protoShowerW(protoShowerVectorW.front());

    const bool isElectronPathwayU(pAlgorithm->IsElectronPathway(protoShowerU.m_connectionPathway.m_pathwayHitList));
    const bool isElectronPathwayV(pAlgorithm->IsElectronPathway(protoShowerV.m_connectionPathway.m_pathwayHitList));
    const bool isElectronPathwayW(pAlgorithm->IsElectronPathway(protoShowerW.m_connectionPathway.m_pathwayHitList));

    if ((isElectronPathwayU && isElectronPathwayV) || (isElectronPathwayU && isElectronPathwayW) || (isElectronPathwayV && isElectronPathwayW))
        return true;
    */
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ElectronStartRefinementTool::IsShowerExtendable(ShowerStartRefinementAlgorithm *const pAlgorithm, const ElectronProtoShowerVector &protoShowerVectorU, 
    const ElectronProtoShowerVector &protoShowerVectorV, const ElectronProtoShowerVector &protoShowerVectorW)
{
    // Return if have not found a connection pathway in each view
    if ((protoShowerVectorU.size() != 1) || (protoShowerVectorV.size() != 1) || (protoShowerVectorW.size() != 1))
        return false;

    const ElectronProtoShower &protoShowerU(protoShowerVectorU.front());
    const ElectronProtoShower &protoShowerV(protoShowerVectorV.front());
    const ElectronProtoShower &protoShowerW(protoShowerVectorW.front());

    const CaloHitList &hitsToAddU(protoShowerU.m_hitsToAdd);
    const CaloHitList &hitsToAddV(protoShowerV.m_hitsToAdd);
    const CaloHitList &hitsToAddW(protoShowerW.m_hitsToAdd);

    // Return if no hits to add...
    if (hitsToAddU.empty() && hitsToAddV.empty() && hitsToAddW.empty())
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//DUDE what if i change this to be like, does the 3D direction match the 3D shower??? since the direction comparison is pretty tempremental? 
// could then even pick up two view issues? if that turns out to be a problem?
bool ElectronStartRefinementTool::ArePathwaysConsistent(ShowerStartRefinementAlgorithm *const pAlgorithm, const ProtoShower &protoShowerU, 
    const ProtoShower &protoShowerV, const ProtoShower &protoShowerW)
{
    if (this->AreShowerStartsConsistent(pAlgorithm, protoShowerU, protoShowerV, protoShowerW, 2.f))
        return true;
    else if (this->AreDirectionsConsistent(pAlgorithm, protoShowerU, protoShowerV, protoShowerW, 5.f))
        return true;
    else
        return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ElectronStartRefinementTool::AreShowerStartsConsistent(ShowerStartRefinementAlgorithm *const pAlgorithm, const ProtoShower &protoShowerU, 
    const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const float allowance)
{
    const CartesianVector &showerStartU(protoShowerU.m_showerCore.m_startPosition);
    const CartesianVector &showerStartV(protoShowerV.m_showerCore.m_startPosition);
    const CartesianVector &showerStartW(protoShowerW.m_showerCore.m_startPosition);

    const float xSeparationUV(std::fabs(showerStartU.GetX() - showerStartV.GetX()));
    const float xSeparationUW(std::fabs(showerStartU.GetX() - showerStartW.GetX()));
    const float xSeparationVW(std::fabs(showerStartV.GetX() - showerStartW.GetX()));

    if ((xSeparationUV > m_maxXSeparation) || (xSeparationUW > m_maxXSeparation) || (xSeparationVW > m_maxXSeparation))
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

    if (((separationU + separationV + separationW) / 3.f) > allowance)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ElectronStartRefinementTool::AreDirectionsConsistent(ShowerStartRefinementAlgorithm *const pAlgorithm, const ProtoShower &protoShowerU, 
    const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const float allowance)
{
    const CartesianVector &directionU(protoShowerU.m_connectionPathway.m_startDirection);
    const CartesianVector &directionV(protoShowerV.m_connectionPathway.m_startDirection);
    const CartesianVector &directionW(protoShowerW.m_connectionPathway.m_startDirection);

    if (directionU.GetX() * directionV.GetX() < 0.f)
        return false;

    if (directionU.GetX() * directionW.GetX() < 0.f)
        return false;

    if (directionV.GetX() * directionW.GetX() < 0.f)
        return false;

    const CartesianVector projectionU(LArGeometryHelper::MergeTwoDirections(pAlgorithm->GetPandora(), TPC_VIEW_V, TPC_VIEW_W, directionV, directionW));
    const CartesianVector projectionV(LArGeometryHelper::MergeTwoDirections(pAlgorithm->GetPandora(), TPC_VIEW_W, TPC_VIEW_U, directionW, directionU));
    const CartesianVector projectionW(LArGeometryHelper::MergeTwoDirections(pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, directionU, directionV));

    const float openingAngleU(directionU.GetOpeningAngle(projectionU) * 180.f / M_PI);
    const float openingAngleV(directionV.GetOpeningAngle(projectionV) * 180.f / M_PI);
    const float openingAngleW(directionW.GetOpeningAngle(projectionW) * 180.f / M_PI);

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

    if (((openingAngleU + openingAngleV + openingAngleW) / 3.f) > allowance)
        return false;

    /*
    if ((openingAngleU > allowance) || (openingAngleV > allowance) || (openingAngleW > allowance))
        return false;
    */
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronStartRefinementTool::ExtendShower(ShowerStartRefinementAlgorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo,
    const ElectronProtoShowerVector &protoShowerVectorU, const ElectronProtoShowerVector &protoShowerVectorV, const ElectronProtoShowerVector &protoShowerVectorW)
{
    const CaloHitList &hitsToAddU(protoShowerVectorU.front().m_hitsToAdd);
    const CaloHitList &hitsToAddV(protoShowerVectorV.front().m_hitsToAdd);
    const CaloHitList &hitsToAddW(protoShowerVectorW.front().m_hitsToAdd);

    ElectronProtoShowerVector showersToExtend;

    if (!hitsToAddU.empty())
        showersToExtend.push_back(protoShowerVectorU.front());

    if (!hitsToAddV.empty())
        showersToExtend.push_back(protoShowerVectorV.front());

    if (!hitsToAddW.empty())
        showersToExtend.push_back(protoShowerVectorW.front());

    if (showersToExtend.size() == 1)
        this->ExtendShowerOneView(pAlgorithm, pShowerPfo, showersToExtend);
    else if (showersToExtend.size() == 2)
        this->ExtendShowerTwoView(pAlgorithm, pShowerPfo, showersToExtend);
    else if (showersToExtend.size() == 3)
        this->ExtendShowerThreeView(pAlgorithm, pShowerPfo, showersToExtend);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronStartRefinementTool::ExtendShowerOneView(ShowerStartRefinementAlgorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo, 
    const ElectronProtoShowerVector &showersToExtend)
{
    const CaloHitList &hitsToAdd(showersToExtend.front().m_hitsToAdd);

    /////////////////////////////////
    /*
    std::cout << "ONE VIEW HIT ADDITION" << std::endl;
    for (const CaloHit *const pCaloHit : hitsToAdd)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "ADDED HIT", BLUE, 2);
    }
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    */
    /////////////////////////////////

    pAlgorithm->AddElectronPathway(pShowerPfo, hitsToAdd);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronStartRefinementTool::ExtendShowerTwoView(ShowerStartRefinementAlgorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo, 
    const ElectronProtoShowerVector &showersToExtend)
{
    const ElectronProtoShower &protoShower1(showersToExtend[0]), &protoShower2(showersToExtend[1]);

    const CaloHitList &hitsToAdd1(protoShower1.m_hitsToAdd);
    const CaloHitList &hitsToAdd2(protoShower2.m_hitsToAdd);

    /////////////////////////////////
    /*
    std::cout << "TWO VIEW HIT ADDITION" << std::endl;
    for (const CaloHit *const pCaloHit : hitsToAdd1)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "ADDED HIT", BLUE, 2);
    }

    for (const CaloHit *const pCaloHit : hitsToAdd2)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "ADDED HIT", BLUE, 2);
    }

    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    */
    /////////////////////////////////

    pAlgorithm->AddElectronPathway(pShowerPfo, hitsToAdd1);
    pAlgorithm->AddElectronPathway(pShowerPfo, hitsToAdd2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronStartRefinementTool::ExtendShowerThreeView(ShowerStartRefinementAlgorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo, 
    const ElectronProtoShowerVector &showersToExtend)
{
    const ElectronProtoShower &protoShower1(showersToExtend[0]), &protoShower2(showersToExtend[1]), &protoShower3(showersToExtend[2]);

    const CaloHitList &hitsToAdd1(protoShower1.m_hitsToAdd);
    const CaloHitList &hitsToAdd2(protoShower2.m_hitsToAdd);
    const CaloHitList &hitsToAdd3(protoShower3.m_hitsToAdd);

    /////////////////////////////////
    /*
    std::cout << "THREE VIEW HIT ADDITION" << std::endl;
    for (const CaloHit *const pCaloHit : hitsToAdd1)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "ADDED HIT", BLUE, 2);
    }

    for (const CaloHit *const pCaloHit : hitsToAdd2)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "ADDED HIT", BLUE, 2);
    }

    for (const CaloHit *const pCaloHit : hitsToAdd3)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "ADDED HIT", BLUE, 2);
    }

    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    */
    /////////////////////////////////

    pAlgorithm->AddElectronPathway(pShowerPfo, hitsToAdd1);
    pAlgorithm->AddElectronPathway(pShowerPfo, hitsToAdd2);
    pAlgorithm->AddElectronPathway(pShowerPfo, hitsToAdd3);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ElectronStartRefinementTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxCoincideneTransverseSeparation", m_maxCoincideneTransverseSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinSpinePurity", m_minSpinePurity));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxAngularDeviation", m_maxAngularDeviation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxXSeparation", m_maxXSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ExtendMode", m_extendMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MoveVertexMode", m_moveVertexMode));

    ShowerStartRefinementBaseTool::ReadSettings(xmlHandle);

    return STATUS_CODE_SUCCESS;
}
}

//------------------------------------------------------------------------------------------------------------------------------------------


            /* // be a bit nicer
            CartesianVector showerExtremalPoint(0.f, 0.f, 0.f);
            showerExtremalPoint = isEndDownstream ? twoDShowerSlidingFit.GetShowerFitResult().GetGlobalMaxLayerPosition() : 
                twoDShowerSlidingFit.GetShowerFitResult().GetGlobalMinLayerPosition();

            CartesianVector displacement(showerStartPosition - showerExtremalPoint);

            const float l(peakDirection.GetDotProduct(displacement));

            if (l > -1.f)
            {
                PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerStartPosition, "Shower Start", BLACK, 2);
                PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerExtremalPoint, "Shower Extremal", GREEN, 2);
                std::cout << "l: " << l << std::endl;
                /////////////////////////////////
                std::cout << "NO SHOWER START FOUND" << std::endl;
                PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
                /////////////////////////////////
                continue;
            }
            */
