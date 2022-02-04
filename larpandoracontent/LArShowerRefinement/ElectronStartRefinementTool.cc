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
    m_minSpinePurity(0.7f)
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

    ProtoShowerVector protoShowerVectorU, protoShowerVectorV, protoShowerVectorW;

    this->BuildProtoShowers(pAlgorithm, pShowerPfo, nuVertexPosition, TPC_VIEW_U, protoShowerVectorU);
    this->BuildProtoShowers(pAlgorithm, pShowerPfo, nuVertexPosition, TPC_VIEW_V, protoShowerVectorV);
    this->BuildProtoShowers(pAlgorithm, pShowerPfo, nuVertexPosition, TPC_VIEW_W, protoShowerVectorW);

    // can the shower be extended?
    if (!this->IsShowerExtendable(pAlgorithm, pShowerPfo, protoShowerVectorU, protoShowerVectorV, protoShowerVectorW))
        return false;

    this->ExtendShower(pAlgorithm, pShowerPfo, protoShowerVectorU, protoShowerVectorV, protoShowerVectorW);

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronStartRefinementTool::BuildProtoShowers(ShowerStartRefinementAlgorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo,
    const CartesianVector &nuVertexPosition, const HitType hitType, ProtoShowerVector &protoShowerVector)
{
    /////////////////////////////////
    std::cout << "Investigating " << (hitType == TPC_VIEW_U ? "U" : hitType == TPC_VIEW_V ? "V" : "W") << " view..." << std::endl;
    /////////////////////////////////

    CaloHitList viewShowerHitList;
    const CaloHitList viewAllHitList(this->GetAllHitsOfType(pAlgorithm, hitType));
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

    // Now, for future use, find the shower's vertex
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
        const CartesianVector end(projectedNuVertexPosition + (peakDirection * 10.f));
        PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &projectedNuVertexPosition, &end, "direction", BLACK, 2,2);
        ///////////////////////////////////////////

        // Find the shower spine
        CaloHitList showerSpineHitList;
        try
        {
            CaloHitList unavailableHits;
            this->FindShowerSpine(pAlgorithm, viewAllHitList, projectedNuVertexPosition, peakDirection, unavailableHits, showerSpineHitList);
        }
        catch (...)
        {
            ///////////////////////////////////////////  
            std::cout << "SHOWER SPINE FIT FAILED" << std::endl;
            PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
            ///////////////////////////////////////////  
            continue;
        }

        // Demand that spine is significant, be lenient here as some have small stubs and a gap
        if (showerSpineHitList.size() < 7)
        {
            ///////////////////////////////////////////  
            std::cout << "SHOWER NOT SIGNIFICANT" << std::endl;
            PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
            ///////////////////////////////////////////  
            continue;
        }

        // Does spine coincide with shower pfo
        if (!this->IsSpineCoincident(pAlgorithm, projectedNuVertexPosition, peakDirection, viewShowerHitList, showerVertexPosition, showerSpineHitList))
        {
            ///////////////////////////////////////////  
            std::cout << "not coincident" << std::endl;
            PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
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
            std::cout << "the super fine fit failed probably a gap in the initial hits??" << std::endl;
            PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
            /////////////////////////////////
            continue;
        }

        // Obtain spine energy profile
        EnergySpectrumMap energySpectrumMap;
        this->GetEnergyDistribution(pAlgorithm, showerSpineHitList, longitudinalPositionMap, energySpectrumMap);

        // Postulate combining the two
        CaloHitList mergedHitList(viewShowerHitList);
        mergedHitList.insert(mergedHitList.begin(), viewShowerHitList.begin(), viewShowerHitList.end());

        // I don't like the ProtoShowers being created here
        CartesianVector showerStartPosition(0.f, 0.f, 0.f);
        this->FindShowerStart(pAlgorithm, projectedNuVertexPosition, peakDirection, longitudinalPositionMap, energySpectrumMap, showerSpineHitList, 
            showerStartPosition, mergedHitList, isEndDownstream, protoShowerVector); 

        /////////////////////////////////
        PfoList visualize;
        visualize.push_back(pShowerPfo);
        PandoraMonitoringApi::VisualizeParticleFlowObjects(pAlgorithm->GetPandora(), &visualize, "ShowerPfo", BLUE);
        PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &showerSpineHitList, "Shower Spine", BLACK);
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &projectedNuVertexPosition, "Nu vertex position", GREEN, 2);

        CaloHitList hitsToAdd;
        this->GetHitsToAdd(pAlgorithm, protoShowerVector.front(), pShowerPfo, hitsToAdd);

        for (const CaloHit *const pHitToAdd : hitsToAdd)
        {
            const CartesianVector &hitPosition(pHitToAdd->GetPositionVector());
            PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "HIT TO ADD", RED, 2);
        }
        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
        /////////////////////////////////

        break;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

CaloHitList ElectronStartRefinementTool::GetAllHitsOfType(ShowerStartRefinementAlgorithm *const pAlgorithm, const HitType hitType)
{
    CaloHitList viewHitList;

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*pAlgorithm, "CaloHitList2D", pCaloHitList));

    if (!pCaloHitList || pCaloHitList->empty())
    {
        if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
            std::cout << "ElectronStartRefinementTool: unable to find calo hit list " << "CaloHitList2D" << std::endl;

        return viewHitList;
    }

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        if (pCaloHit->GetHitType() == hitType)
            viewHitList.push_back(pCaloHit);
    }

    return viewHitList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronStartRefinementTool::CollectHitsWithinROI(ShowerStartRefinementAlgorithm *const pAlgorithm, const CaloHitList &showerHitList, 
    const CaloHitList &allHitList, const CartesianVector &projectedNuVertexPosition, CaloHitList &collectedHits)
{
    float lowestTheta(std::numeric_limits<float>::max()), highestTheta((-1.f) * std::numeric_limits<float>::max());

    this->GetAngularExtrema(pAlgorithm, showerHitList, projectedNuVertexPosition, lowestTheta, highestTheta);
    this->CollectHitsWithinExtrema(pAlgorithm, allHitList, projectedNuVertexPosition, lowestTheta, highestTheta, collectedHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

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
            theta0XZ += 3.14;

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
            theta0XZ += 3.14;

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
    CartesianVector end(projectedNuVertexPosition + (peakDirection * 100.0));
    PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &projectedNuVertexPosition, &end, "DIRECTION", PINK, 2, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerVertexPosition, "SPINE VERTEX POSITION", GREEN, 2);
    /////////////////////////////////

    // Check whether shower is an extension of found pathway
    const CartesianVector displacement(showerVertexPosition - projectedNuVertexPosition);
    const float transverseSeparation((peakDirection.GetCrossProduct(displacement)).GetMagnitude());

    if (transverseSeparation > m_maxCoincideneTransverseSeparation)
    {
        /////////////////////////////////
        std::cout << "transverse separation issue" << std::endl;
        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
        /////////////////////////////////
        return false;
    }

    // Check whether shower spine is pure
    const CaloHitList sharedHitList(LArMCParticleHelper::GetSharedHits(postShowerVertexSpineHits, viewShowerHitList));
    const float spinePurity(postShowerVertexSpineHits.size() == 0 ? 1.f : static_cast<float>(sharedHitList.size()) / static_cast<float>(postShowerVertexSpineHits.size()));

    if (spinePurity < m_minSpinePurity)
    {
        /////////////////////////////////
        std::cout <<  "spine purity issue" << std::endl;
        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
        /////////////////////////////////
        return false;
    }

    /////////////////////////////////
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    /////////////////////////////////

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ElectronStartRefinementTool::IsShowerExtendable(ShowerStartRefinementAlgorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo,
    const ProtoShowerVector &protoShowerVectorU, const ProtoShowerVector &protoShowerVectorV, const ProtoShowerVector &protoShowerVectorW)
{
    // throw if have a size more than one?
    if (protoShowerVectorU.empty() || protoShowerVectorV.empty() || protoShowerVectorW.empty())
        return false;

    CaloHitList hitsToAddU, hitsToAddV, hitsToAddW;

    this->GetHitsToAdd(pAlgorithm, protoShowerVectorU.front(), pShowerPfo, hitsToAddU);
    this->GetHitsToAdd(pAlgorithm, protoShowerVectorV.front(), pShowerPfo, hitsToAddV);
    this->GetHitsToAdd(pAlgorithm, protoShowerVectorW.front(), pShowerPfo, hitsToAddW);

    // Return if no hits to add...
    if (hitsToAddU.empty() && hitsToAddV.empty() && hitsToAddW.empty())
        return false;

    // Check whether angular peak direction agrees?
    const CartesianVector &peakDirectionU(protoShowerVectorU.front().m_connectionPathway.m_startDirection);
    const CartesianVector &peakDirectionV(protoShowerVectorV.front().m_connectionPathway.m_startDirection);
    const CartesianVector &peakDirectionW(protoShowerVectorW.front().m_connectionPathway.m_startDirection);

    if (!this->AreDirectionsConsistent(pAlgorithm, peakDirectionU, peakDirectionV, peakDirectionW, 3.f))
        return false;

    // Shower start is unlikely to be in exacty the same place, so check that the general direction agrees
    const CartesianVector &pathwayDirectionU(protoShowerVectorU.front().m_showerCore.m_startPosition - protoShowerVectorU.front().m_connectionPathway.m_startPosition);
    const CartesianVector &pathwayDirectionV(protoShowerVectorV.front().m_showerCore.m_startPosition - protoShowerVectorV.front().m_connectionPathway.m_startPosition);
    const CartesianVector &pathwayDirectionW(protoShowerVectorW.front().m_showerCore.m_startPosition - protoShowerVectorW.front().m_connectionPathway.m_startPosition);

    if (!this->AreDirectionsConsistent(pAlgorithm, pathwayDirectionU, pathwayDirectionV, pathwayDirectionW, 5.f))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronStartRefinementTool::GetHitsToAdd(ShowerStartRefinementAlgorithm *const pAlgorithm, const ProtoShower &protoShower, 
    const ParticleFlowObject *const pShowerPfo, CaloHitList &hitsToAdd)
{
    const HitType hitType(protoShower.m_connectionPathway.m_pathwayHitList.front()->GetHitType());

    CaloHitList showerHitList;
    LArPfoHelper::GetCaloHits(pShowerPfo, hitType, showerHitList);

    for (const CaloHit *const pPathwayHit : protoShower.m_connectionPathway.m_pathwayHitList)
    {
        if (std::find(showerHitList.begin(), showerHitList.end(), pPathwayHit) == showerHitList.end())
            hitsToAdd.push_back(pPathwayHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ElectronStartRefinementTool::AreDirectionsConsistent(ShowerStartRefinementAlgorithm *const pAlgorithm, const CartesianVector &directionU, 
    const CartesianVector &directionV, const CartesianVector &directionW, const float allowance)
{
    const CartesianVector projectionU(LArGeometryHelper::MergeTwoDirections(pAlgorithm->GetPandora(), TPC_VIEW_V, TPC_VIEW_W, directionV, directionW));
    const CartesianVector projectionV(LArGeometryHelper::MergeTwoDirections(pAlgorithm->GetPandora(), TPC_VIEW_W, TPC_VIEW_U, directionW, directionU));
    const CartesianVector projectionW(LArGeometryHelper::MergeTwoDirections(pAlgorithm->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, directionU, directionV));

    const float openingAngleU(directionU.GetOpeningAngle(projectionU) * 180 / M_PI);
    const float openingAngleV(directionV.GetOpeningAngle(projectionV) * 180 / M_PI);
    const float openingAngleW(directionW.GetOpeningAngle(projectionW) * 180 / M_PI);

    if ((openingAngleU > allowance) || (openingAngleV > allowance) || (openingAngleW > allowance))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronStartRefinementTool::ExtendShower(ShowerStartRefinementAlgorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo,
    const ProtoShowerVector &protoShowerVectorU, const ProtoShowerVector &protoShowerVectorV, const ProtoShowerVector &protoShowerVectorW)
{
    CaloHitList hitsToAddU, hitsToAddV, hitsToAddW;

    this->GetHitsToAdd(pAlgorithm, protoShowerVectorU.front(), pShowerPfo, hitsToAddU);
    this->GetHitsToAdd(pAlgorithm, protoShowerVectorV.front(), pShowerPfo, hitsToAddV);
    this->GetHitsToAdd(pAlgorithm, protoShowerVectorW.front(), pShowerPfo, hitsToAddW);

    ProtoShowerVector showersToExtend;

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
    const ProtoShowerVector &showersToExtend)
{
    // don't check whether it is an electron pathway? just go ahead? <- maybe have a cautious mode?

    CaloHitList hitsToAdd;
    this->GetHitsToAdd(pAlgorithm, showersToExtend.front(), pShowerPfo, hitsToAdd);

    pAlgorithm->AddElectronPathway(pShowerPfo, hitsToAdd);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronStartRefinementTool::ExtendShowerTwoView(ShowerStartRefinementAlgorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo, 
    const ProtoShowerVector &showersToExtend)
{
    const ProtoShower &protoShower1(showersToExtend[0]), &protoShower2(showersToExtend[1]);

    CaloHitList hitsToAdd1, hitsToAdd2;
    this->GetHitsToAdd(pAlgorithm, protoShower1, pShowerPfo, hitsToAdd1);
    this->GetHitsToAdd(pAlgorithm, protoShower2, pShowerPfo, hitsToAdd2);

    const bool isElectronPathway1(pAlgorithm->IsElectronPathway(hitsToAdd1));
    const bool isElectronPathway2(pAlgorithm->IsElectronPathway(hitsToAdd2));

    if (isElectronPathway1 && isElectronPathway2)
    {
        pAlgorithm->AddElectronPathway(pShowerPfo, hitsToAdd1);
        pAlgorithm->AddElectronPathway(pShowerPfo, hitsToAdd2);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronStartRefinementTool::ExtendShowerThreeView(ShowerStartRefinementAlgorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo, 
    const ProtoShowerVector &showersToExtend)
{
    const ProtoShower &protoShower1(showersToExtend[0]), &protoShower2(showersToExtend[1]), &protoShower3(showersToExtend[2]);

    CaloHitList hitsToAdd1, hitsToAdd2, hitsToAdd3;
    this->GetHitsToAdd(pAlgorithm, protoShower1, pShowerPfo, hitsToAdd1);
    this->GetHitsToAdd(pAlgorithm, protoShower2, pShowerPfo, hitsToAdd2);
    this->GetHitsToAdd(pAlgorithm, protoShower3, pShowerPfo, hitsToAdd3);

    const bool isElectronPathway1(pAlgorithm->IsElectronPathway(hitsToAdd1));
    const bool isElectronPathway1(pAlgorithm->IsElectronPathway(hitsToAdd2));
    const bool isElectronPathway1(pAlgorithm->IsElectronPathway(hitsToAdd3));

    if ((isElectronPathway1 && isElectronPathway2) || (isElectronPathway1 && isElectronPathway3) || (isElectronPathway2 && isElectronPathway3))
    {
        pAlgorithm->AddElectronPathway(pShowerPfo, hitsToAdd1);
        pAlgorithm->AddElectronPathway(pShowerPfo, hitsToAdd2);
        pAlgorithm->AddElectronPathway(pShowerPfo, hitsToAdd3);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ElectronStartRefinementTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxCoincideneTransverseSeparation", m_maxCoincideneTransverseSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinSpinePurity", m_minSpinePurity));

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
