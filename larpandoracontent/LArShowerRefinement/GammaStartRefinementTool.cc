/**
 *  @file   larpandoracontent/LArShowerRefinement/GammaStartRefinementTool.cc
 *
 *  @brief  Implementation of the gamma start refinement tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/AlgorithmTool.h"

#include "PandoraMonitoringApi.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingShowerFitResult.h"

#include "larpandoracontent/LArShowerRefinement/GammaStartRefinementTool.h"

using namespace pandora;

namespace lar_content
{

GammaStartRefinementTool::GammaStartRefinementTool() : 
    m_counter(0), 
    m_theta0XZBinSize(0.005f),
    m_pathwaySearchRegion(14.f),
    m_smoothingWindow(3),
    m_showerCounter(0), 
    m_microSlidingFitWindow(20),
    m_minSigmaDeviation(5.f),
    m_trackSearchWindow(5.f),
    m_nInitialEnergyBins(5),
    m_minTrackBlipMean(3.f),
    m_showerSlidingFitWindow(20),
    m_molliereRadius(9.f),
    m_minShowerOpeningAngle(3.f)
{
}

GammaStartRefinementTool::~GammaStartRefinementTool()
{
}
    //// ISOBEL NEED TO DECIDE TO DO NOTHING IF PATHWAY IS LIKE 80% OF TOTAL PARTICLE FLOW OBJECT
//------------------------------------------------------------------------------------------------------------------------------------------

bool GammaStartRefinementTool::Run(ShowerStartRefinementAlgorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo, const CartesianVector &nuVertexPosition)
{
    // only apply gamma refinement algorithm to showers
    if (!LArPfoHelper::IsShower(pShowerPfo))
        return false;

    if (!this->HasPathToNuVertex(pShowerPfo, nuVertexPosition))
        return false;

    ////////////////////////////////
    // temporary  - actually works really well so keep this
    CaloHitList caloHits3D;
    LArPfoHelper::GetCaloHits(pShowerPfo, TPC_3D, caloHits3D);

    if (caloHits3D.size() < 100)
        return false;
    ////////////////////////////////

    ProtoShowerVector protoShowerVectorU, protoShowerVectorV, protoShowerVectorW;

    this->BuildProtoShowers(pAlgorithm, pShowerPfo, nuVertexPosition, TPC_VIEW_U, protoShowerVectorU);

    std::cout << "CONNECTION PATHWAYS U" << std::endl;
    for (const ProtoShower protoShower : protoShowerVectorU)
    {
        const CaloHitList &connectionPathway(protoShower.m_connectionPathway.m_pathwayHitList);
        const CartesianVector &showerStartPosition(protoShower.m_showerCore.m_startPosition);

        PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &connectionPathway, "ConnectionPathway", BLACK);
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerStartPosition, "ShowerStartPosition", BLACK, 2);
    }
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());

    this->BuildProtoShowers(pAlgorithm, pShowerPfo, nuVertexPosition, TPC_VIEW_V, protoShowerVectorV);

    std::cout << "CONNECTION PATHWAYS V" << std::endl;
    for (const ProtoShower protoShower : protoShowerVectorV)
    {
        const CaloHitList &connectionPathway(protoShower.m_connectionPathway.m_pathwayHitList);
        const CartesianVector &showerStartPosition(protoShower.m_showerCore.m_startPosition);

        PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &connectionPathway, "ConnectionPathway", BLACK);
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerStartPosition, "ShowerStartPosition", BLACK, 2);
    }

    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());

    this->BuildProtoShowers(pAlgorithm, pShowerPfo, nuVertexPosition, TPC_VIEW_W, protoShowerVectorW);

    std::cout << "CONNECTION PATHWAYS W" << std::endl;
    for (const ProtoShower protoShower : protoShowerVectorW)
    {
        const CaloHitList &connectionPathway(protoShower.m_connectionPathway.m_pathwayHitList);
        const CartesianVector &showerStartPosition(protoShower.m_showerCore.m_startPosition);

        PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &connectionPathway, "ConnectionPathway", BLACK);
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerStartPosition, "ShowerStartPosition", BLACK, 2);
    }
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
  
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::BuildProtoShowers(ShowerStartRefinementAlgorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo,
    const CartesianVector &nuVertexPosition, const HitType tpcView, ProtoShowerVector &protoShowerVector)
{
    std::cout << "Investigating " << (tpcView == TPC_VIEW_U ? "U" : tpcView == TPC_VIEW_V ? "V" : "W") << " view..." << std::endl;

    CaloHitList viewShowerHitList;
    LArPfoHelper::GetCaloHits(pShowerPfo, tpcView, viewShowerHitList);

    const CartesianVector projectedNuVertexPosition(LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertexPosition, tpcView));

    // Fill angular decomposition map
    AngularDecompositionMap angularDecompositionMap;
    this->FillAngularDecompositionMap(pAlgorithm, viewShowerHitList, projectedNuVertexPosition, angularDecompositionMap);

    if (angularDecompositionMap.empty())
        return;

    this->SmoothAngularDecompositionMap(angularDecompositionMap);

    // Obtain distribution peaks
    IntVector angularPeakVector;
    this->ObtainPeakVector(angularDecompositionMap, angularPeakVector);

    // Investigate peaks, searching for pathways (peaks indexed by theta0XZ bin lower edge)
    IntVector investigatedPeaks;

    // Keep track of 'taken' hits
    CaloHitList unavailableHits;

    for (unsigned int i = 0; i < angularPeakVector.size(); ++i)
    {
        int bestTheta0XZBin(0);

        if (!this->FindBestAngularPeak(angularDecompositionMap, angularPeakVector, investigatedPeaks, bestTheta0XZBin))
            break;

        investigatedPeaks.push_back(bestTheta0XZBin);

        const float theta0XZ(bestTheta0XZBin * m_theta0XZBinSize);
        const CartesianVector peakDirection(std::cos(theta0XZ), 0.f, std::sin(theta0XZ));
        const bool isEndDownstream(peakDirection.GetZ() > 0.f);

        //////////////////////////////////////
        //std::cout << "bestTheta0XZBin: " << bestTheta0XZBin << std::endl;
        //std::cout << "peakDirection: " << peakDirection << std::endl;
        //const CartesianVector projection(projectedNuVertexPosition + (peakDirection * 14.f));
        //PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &projectedNuVertexPosition, &projection, "direction", BLACK, 2, 1);
        //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
        //////////////////////////////////////

        CaloHitList showerSpineHitList;
        this->FindShowerSpine(pAlgorithm, viewShowerHitList, projectedNuVertexPosition, peakDirection, unavailableHits, showerSpineHitList);

        // Demand that spine is significant
        if (showerSpineHitList.size() < 20)
        {
            std::cout << "Found shower spine is insignificant" << std::endl;
            continue;
        }

        // Obtail longitudinal position of spine hits
        LongitudinalPositionMap longitudinalPositionMap;
        this->ObtainLongitudinalDecomposition(pAlgorithm, showerSpineHitList, longitudinalPositionMap);

        // Obtain spine energy profile
        EnergySpectrumMap energySpectrumMap;
        this->GetEnergyDistribution(pAlgorithm, showerSpineHitList, longitudinalPositionMap, energySpectrumMap);

        // Find shower start position - and pathway! - and possible post shower hits? (ISOBEL TODO)
        CartesianVector showerStartPosition(0.f, 0.f, 0.f);
        this->FindShowerStart(pAlgorithm, longitudinalPositionMap, energySpectrumMap, showerSpineHitList, showerStartPosition, viewShowerHitList, isEndDownstream, protoShowerVector);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::FillAngularDecompositionMap(ShowerStartRefinementAlgorithm *const pAlgorithm, const CaloHitList &viewShowerHitList, 
    const CartesianVector &projectedNuVertexPosition, AngularDecompositionMap &angularDecompositionMap)
{
    const CartesianVector xAxis(1.f, 0.f, 0.f);

    for (const CaloHit *const pCaloHit : viewShowerHitList)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        const CartesianVector displacementVector(hitPosition - projectedNuVertexPosition);

        if (displacementVector.GetMagnitudeSquared() > (m_pathwaySearchRegion * m_pathwaySearchRegion))
            continue;

        float theta0XZ(xAxis.GetOpeningAngle(displacementVector));

        if (displacementVector.GetZ() < 0.f)
            theta0XZ *= -1.f;

        const int theta0XZFactor(std::floor(theta0XZ / m_theta0XZBinSize));

        if (angularDecompositionMap.find(theta0XZFactor) == angularDecompositionMap.end())
            angularDecompositionMap[theta0XZFactor] = 1;
        else
            angularDecompositionMap[theta0XZFactor] += 1;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::SmoothAngularDecompositionMap(AngularDecompositionMap &angularDecompositionMap)
{
    //for (auto &entry : angularDecompositionMap)
    //std::cout << "binNo: " << entry.first << ", calo hit list size: " << entry.second << std::endl;

    const int loopMin = (-1) * (m_smoothingWindow - 1) / 2;
    const int loopMax = (m_smoothingWindow - 1) / 2;

    const AngularDecompositionMap angularDecompositionMapTemp(angularDecompositionMap);

    angularDecompositionMap.clear();

    for (const auto &entry : angularDecompositionMapTemp)
    {
        const int currentBin(entry.first);
        float total(0.f);
        int binCount(0);

        for (int binOffset = loopMin; binOffset <= loopMax; ++binOffset)
        {
            ++binCount;

            const int contributingBin = currentBin + binOffset;
            total += (angularDecompositionMapTemp.find(contributingBin) == angularDecompositionMapTemp.end()) ? 0.f : angularDecompositionMapTemp.at(contributingBin);
        }

        angularDecompositionMap[currentBin] = total / static_cast<float>(binCount);
    }

    //std::cout << "AFTER SMOOTHING" << std::endl;
    //for (auto &entry : angularDecompositionMap)
    //std::cout << "binNo: " << entry.first << ", calo hit list size: " << entry.second << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::ObtainPeakVector(AngularDecompositionMap &angularDecompositionMap, IntVector &viewPeakVector)
{
    for (const auto &entry : angularDecompositionMap)
    {
        const float bin(entry.first);
        const float binWeight(entry.second);

        int precedingBin(bin - 1);
        bool foundPreceeding(false);

        while (!foundPreceeding)
        {
            if ((angularDecompositionMap.find(precedingBin) == angularDecompositionMap.end()) ||
                (std::fabs(angularDecompositionMap.at(precedingBin) - angularDecompositionMap.at(bin)) > std::numeric_limits<float>::epsilon()))
            {
                foundPreceeding = true;
                break;
            }

            --precedingBin;
        }

        // ISOBEL SHOULD ALWAYS FIND ONE!?
        if (!foundPreceeding)
        {
            std::cout << "ISOBEL DID NOT FIND PRECEEDING BIN" << std::endl;
            throw;
        }

        int followingBin(bin + 1);
        bool foundFollowing(false);

        while (!foundFollowing)
        {
            if ((angularDecompositionMap.find(followingBin) == angularDecompositionMap.end()) ||
                (std::fabs(angularDecompositionMap.at(followingBin) - angularDecompositionMap.at(bin)) > std::numeric_limits<float>::epsilon()))
            {
                foundFollowing = true;
                break;
            }

            ++followingBin;
        }

        // ISOBEL SHOULD ALWAYS FIND ONE!?
        if (!foundFollowing)
        {
            std::cout << "SIOBEL DID NOT FIND FOLLOWING BIN" << std::endl;
            throw;
        }

        const float precedingBinWeight(angularDecompositionMap.find(precedingBin) == angularDecompositionMap.end() ? 0.f : angularDecompositionMap.at(precedingBin));
        const float followingBinWeight(angularDecompositionMap.find(followingBin) == angularDecompositionMap.end() ? 0.f : angularDecompositionMap.at(followingBin));

        if ((binWeight < precedingBinWeight) || (binWeight < followingBinWeight))
            continue;

        viewPeakVector.push_back(bin);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool GammaStartRefinementTool::FindBestAngularPeak(AngularDecompositionMap &angularDecompositionMap, IntVector &viewPeakVector, 
    IntVector &investigatedPeaks, int &bestTheta0XZBin)
{
    bool found(false);
    float bestWeight(0.f);

    for (const int theta0XZBin : viewPeakVector)
    {
        if (std::find(investigatedPeaks.begin(), investigatedPeaks.end(), theta0XZBin) != investigatedPeaks.end())
            continue;

        if (angularDecompositionMap.find(theta0XZBin) == angularDecompositionMap.end())
        {
            std::cout << "ISOBEL SHOULD NEVER HAPPEN" << std::endl;
            throw;
        }

        const float binWeight(angularDecompositionMap.at(theta0XZBin));

        if (binWeight > bestWeight)
        {
            found = true;
            bestWeight = binWeight;
            bestTheta0XZBin = theta0XZBin;
        }
    }

    return found;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::GetEnergyDistribution(ShowerStartRefinementAlgorithm *const pAlgorithm, const CaloHitList &showerSpineHitList, 
    const LongitudinalPositionMap &longitudinalPositionMap, EnergySpectrumMap &energySpectrumMap)
{
    float e0(0.f);

    for (const CaloHit *pCaloHit : showerSpineHitList)
        e0 += pCaloHit->GetElectromagneticEnergy();



    for (const CaloHit *pCaloHit : showerSpineHitList)
    {
        const float fractionalEnergy(pCaloHit->GetElectromagneticEnergy() / e0);
        const float projection(longitudinalPositionMap.at(pCaloHit));

        const int longitudinalIndex = std::floor(projection / m_longitudinalCoordinateBinSize);

        if (energySpectrumMap.find(longitudinalIndex) == energySpectrumMap.end())
            energySpectrumMap[longitudinalIndex] = fractionalEnergy;
        else
            energySpectrumMap[longitudinalIndex] += fractionalEnergy;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::ObtainLongitudinalDecomposition(ShowerStartRefinementAlgorithm *const pAlgorithm, const CaloHitList &showerSpineHitList, 
    LongitudinalPositionMap &longitudinalPositionMap)
{
    CartesianPointVector hitPositions;

    for (const CaloHit *const pCaloHit : showerSpineHitList)
        hitPositions.push_back(pCaloHit->GetPositionVector());

    const TwoDSlidingFitResult twoDSlidingFit(&hitPositions, m_microSlidingFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const LayerFitResultMap &layerFitResultMap(twoDSlidingFit.GetLayerFitResultMap());

    // Find hits in each layer
    LayerToHitMap layerToHitMap;

    for (const CaloHit *const pCaloHit : showerSpineHitList)
    {
        float hitL(0.f), hitT(0.f);
        const CartesianVector hitPosition(pCaloHit->GetPositionVector());

        twoDSlidingFit.GetLocalPosition(hitPosition, hitL, hitT);
        layerToHitMap[twoDSlidingFit.GetLayer(hitL)].push_back(pCaloHit);
    }

    // Find the longitudinal distance of each hit along shower spine fit
    float runningDistance(0);

    for (auto iter = layerToHitMap.begin(); iter != layerToHitMap.end(); ++iter)
    {
        const int layer(iter->first);
        std::cout << "layer: " << layer << std::endl;

        //CartesianVector blob(0.f, 0.f, 0.f);
        //twoDSlidingFit.GetGlobalPosition(layerFitResultMap.at(layer).GetL(), layerFitResultMap.at(layer).GetFitT(), blob);
        //PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &blob, std::to_string(layer), BLUE, 2);

        const float layerL(layerFitResultMap.at(layer).GetL());

        const int higherLayer(std::next(iter) == layerToHitMap.end() ? layer : std::next(iter)->first);
        const int middleLayer(layer);
        const int lowerLayer(iter == layerToHitMap.begin() ? layer : std::prev(iter)->first);

        CartesianVector lowerLayerPosition(0.f, 0.f, 0.f), middleLayerPosition(0.f, 0.f, 0.f), higherLayerPosition(0.f, 0.f, 0.f);

        twoDSlidingFit.GetGlobalPosition(layerFitResultMap.at(lowerLayer).GetL(), layerFitResultMap.at(lowerLayer).GetFitT(), lowerLayerPosition);
        twoDSlidingFit.GetGlobalPosition(layerFitResultMap.at(middleLayer).GetL(), layerFitResultMap.at(middleLayer).GetFitT(), middleLayerPosition);
        twoDSlidingFit.GetGlobalPosition(layerFitResultMap.at(higherLayer).GetL(), layerFitResultMap.at(higherLayer).GetFitT(), higherLayerPosition);

        const float layerLength = std::next(iter) == layerToHitMap.end() ? 0.f : iter == layerToHitMap.begin() ? 0.f : (middleLayerPosition - lowerLayerPosition).GetMagnitude();

        for (const CaloHit *const pCaloHit : iter->second)
        {
            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());

            float hitL(0.f), hitT(0.f);
            twoDSlidingFit.GetLocalPosition(hitPosition, hitL, hitT);

            CartesianVector localLowerLayerPosition(0.f, 0.f, 0.f), localHigherLayerPosition(0.f, 0.f, 0.f);

            if (hitL < layerL)
            {
                localLowerLayerPosition = lowerLayerPosition;
                localHigherLayerPosition = iter == layerToHitMap.begin() ? higherLayerPosition : middleLayerPosition;
            }
            else
            {
                localLowerLayerPosition = std::next(iter) == layerToHitMap.end() ? lowerLayerPosition : middleLayerPosition;
                localHigherLayerPosition = higherLayerPosition;
            }

            const CartesianVector displacement((higherLayerPosition - lowerLayerPosition).GetUnitVector());

            float longitudinalDisplament = (hitL > layerL ? layerLength : 0.f);
            longitudinalDisplament += (displacement.GetDotProduct(hitPosition - lowerLayerPosition) + runningDistance);

            longitudinalPositionMap[pCaloHit] = longitudinalDisplament;
        }

        runningDistance += layerLength;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool GammaStartRefinementTool::FindShowerStart(ShowerStartRefinementAlgorithm *const pAlgorithm, const LongitudinalPositionMap &longitudinalPositionMap, 
    const EnergySpectrumMap &energySpectrumMap, const CaloHitList &showerSpineHitList, CartesianVector &showerStartPosition, const CaloHitList &showerPfoHitList, 
    const bool isEndDownstream, ProtoShowerVector &protoShowerVector)
{
    if (energySpectrumMap.size() < (m_nInitialEnergyBins + 1))
        return false;

    float meanEnergy(0.f), energySigma(0.f);
    this->CharacteriseInitialEnergy(energySpectrumMap, meanEnergy, energySigma);

    auto energySpectrumIter(energySpectrumMap.begin());

    for (unsigned int i = 0; i < energySpectrumMap.size(); ++i)
    {
        // Bypass bins used to asses the initial energy
        if (i < m_nInitialEnergyBins)
        {
            ++energySpectrumIter;
            continue;
        }

        const float longitudinalCoordinate(energySpectrumIter->first * m_longitudinalCoordinateBinSize);
        const float energyDeviation((energySpectrumIter->second - meanEnergy) / energySigma);

        std::cout << "energyDeviation (sigma): " << energyDeviation << std::endl;

        // Use energy and local topology to assess whether we are at the shower start
        if ((energyDeviation > m_minSigmaDeviation) && this->IsShowerTopology(pAlgorithm, longitudinalCoordinate,
            showerPfoHitList, showerSpineHitList, isEndDownstream))
        {
            break;
        }

        ++energySpectrumIter;
    }

    const int longitudinalStartBin(energySpectrumIter->first);
    const float longitudinalStartCoordinate(longitudinalStartBin * m_longitudinalCoordinateBinSize);

    CartesianVector showerStartDirection(0.f, 0.f, 0.f);
    this->ConvertLongitudinalProjectionToGlobalPosition(pAlgorithm, showerSpineHitList, longitudinalStartCoordinate, showerStartPosition, showerStartDirection);

    std::cout << "longitudinalStartCoordinate: " << longitudinalStartCoordinate << std::endl;

    CaloHitList pathwayHitList;
    for (const CaloHit *const pCaloHit : showerSpineHitList)
    {
        if (longitudinalPositionMap.find(pCaloHit) == longitudinalPositionMap.end())
            continue;

        if (longitudinalPositionMap.at(pCaloHit) > longitudinalStartCoordinate)
            continue;

        std::cout << "longitudinalPositionMap.at(pCaloHit): " << longitudinalPositionMap.at(pCaloHit) << std::endl;

        pathwayHitList.push_back(pCaloHit);
    }

    /*
    for (const CaloHit *const pCaloHit : pathwayHitList)
    {
        const CartesianVector &hit(pCaloHit->GetPositionVector());
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hit, "path", BLUE, 2);
    }
    */

    protoShowerVector.push_back(ProtoShower(ShowerCore(showerStartPosition, showerStartDirection, CaloHitList()), 
        ConnectionPathway(CartesianVector(0.f, 0.f, 0.f), CartesianVector(0.f, 0.f, 0.f), pathwayHitList)));

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::CharacteriseInitialEnergy(const EnergySpectrumMap &energySpectrumMap, float &meanEnergy, float &energySigma)
{
    auto energySpectrumIter(energySpectrumMap.begin());

    for (unsigned int i = 0 ; i < m_nInitialEnergyBins; ++i)
    {
        meanEnergy += energySpectrumIter->second;
        ++energySpectrumIter;
    }

    meanEnergy /= static_cast<float>(m_nInitialEnergyBins);

    energySpectrumIter = energySpectrumMap.begin();

    for (unsigned int i = 0; i < m_nInitialEnergyBins; ++i)
    {
        energySigma += std::pow(energySpectrumIter->second - meanEnergy, 2);
        ++energySpectrumIter;
    }

    energySigma = std::sqrt(energySigma / m_nInitialEnergyBins);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool GammaStartRefinementTool::IsShowerTopology(ShowerStartRefinementAlgorithm *const pAlgorithm, const float longitudinalDistance, 
    const CaloHitList &showerPfoHits, const CaloHitList &showerSpineHitList, const bool isEndDownstream)
{
    CartesianVector showerStartPosition(0.f, 0.f, 0.f), showerStartDirection(0.f, 0.f, 0.f);
    this->ConvertLongitudinalProjectionToGlobalPosition(pAlgorithm, showerSpineHitList,  longitudinalDistance, showerStartPosition, showerStartDirection);

    // Characterise the shower
    bool isBetween(false), doesStraddle(false);
    CartesianVector positiveEdgeStart(0.f, 0.f, 0.f), positiveEdgeEnd(0.f, 0.f, 0.f), positiveEdgeDirection(0.f, 0.f, 0.f);
    CartesianVector negativeEdgeStart(0.f, 0.f, 0.f), negativeEdgeEnd(0.f, 0.f, 0.f), negativeEdgeDirection(0.f, 0.f, 0.f);

    if (this->CharacteriseShower(pAlgorithm, showerPfoHits, showerSpineHitList, showerStartPosition, showerStartDirection, isEndDownstream,
        positiveEdgeStart, positiveEdgeEnd, negativeEdgeStart, negativeEdgeEnd, isBetween, doesStraddle) != STATUS_CODE_SUCCESS)
    {
        std::cout << "failed to characterise the shower" << std::endl;
        return false;
    }

    if (!isBetween)
    {
        std::cout << "shower start is not inbetween start edge points :( " << std::endl;
        return false;
    }

    if (!doesStraddle)
    {
        std::cout << "shower does not straddle core " << std::endl;
        return false;
    }

    positiveEdgeStart = showerStartPosition;
    negativeEdgeStart = showerStartPosition;
    positiveEdgeDirection = positiveEdgeEnd - positiveEdgeStart;
    negativeEdgeDirection = negativeEdgeEnd - negativeEdgeStart;

    const float showerOpeningAngle(positiveEdgeDirection.GetOpeningAngle(negativeEdgeDirection) * 180.f / 3.14);
    const float positiveEdgeDeviation(positiveEdgeDirection.GetOpeningAngle(showerStartDirection) * 180.f / 3.14);
    const float negativeEdgeDeviation(negativeEdgeDirection.GetOpeningAngle(showerStartDirection) * 180.f / 3.14);

    std::cout << "edge opening angle: " << (positiveEdgeDirection.GetOpeningAngle(negativeEdgeDirection) * 180.f / 3.14) << std::endl;
    std::cout << "positive opening angle: " << (positiveEdgeDirection.GetOpeningAngle(showerStartDirection) * 180.f / 3.14) << std::endl;
    std::cout << "negative opening angle: " << (negativeEdgeDirection.GetOpeningAngle(showerStartDirection) * 180.f / 3.14) << std::endl;

    if (showerOpeningAngle < m_minShowerOpeningAngle)
    {
        std::cout << "shower opening angle is too small" << std::endl;
        return false;
    }

    if (positiveEdgeDeviation > 45.f)
    {
        std::cout << "positive edge is not in the direction of the shower" << std::endl;
    }

    if (negativeEdgeDeviation > 45.f)
    {
        std::cout << "negative edge is not in the direction of the shower" << std::endl;
    }

    PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &positiveEdgeStart, &positiveEdgeEnd, "positive direction", BLACK, 2, 1);
    PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &negativeEdgeStart, &negativeEdgeEnd, "negative direction", BLACK, 2, 1);
    //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

    void GammaStartRefinementTool::ConvertLongitudinalProjectionToGlobalPosition(ShowerStartRefinementAlgorithm *const pAlgorithm, const CaloHitList &showerSpineHitList, const float longitudinalDistance, 
    CartesianVector &globalPosition, CartesianVector &globalDirection)
{
    CartesianPointVector hitPositions;

    for (const CaloHit *const pCaloHit : showerSpineHitList)
        hitPositions.push_back(pCaloHit->GetPositionVector());

    const TwoDSlidingFitResult twoDSlidingFit(&hitPositions, m_microSlidingFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const LayerFitResultMap &layerFitResultMap(twoDSlidingFit.GetLayerFitResultMap());

    float runningDistance(0.f);

    CartesianVector previousLayerPosition(0.f, 0.f, 0.f);
    twoDSlidingFit.GetGlobalPosition(layerFitResultMap.begin()->second.GetL(), layerFitResultMap.begin()->second.GetFitT(), previousLayerPosition);

    const CartesianVector jam(twoDSlidingFit.GetGlobalMinLayerPosition());
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &jam, "MIN LAYER POSITION", RED, 2);

    int showerStartLayer(0);

    for (auto iter = std::next(layerFitResultMap.begin()); iter != layerFitResultMap.end(); ++iter)
    {
        const int layer(iter->first);
        showerStartLayer = layer;

        CartesianVector layerPosition(0.f, 0.f, 0.f);
        twoDSlidingFit.GetGlobalPosition(layerFitResultMap.at(layer).GetL(), layerFitResultMap.at(layer).GetFitT(), layerPosition);

        runningDistance += (layerPosition - previousLayerPosition).GetMagnitude();

        if (runningDistance > longitudinalDistance)
            break;

        previousLayerPosition = layerPosition;
    }

    const float lCoordinate(layerFitResultMap.at(showerStartLayer).GetL()), tCoordinate(layerFitResultMap.at(showerStartLayer).GetFitT());
    const float localGradient(layerFitResultMap.at(showerStartLayer).GetGradient());

    twoDSlidingFit.GetGlobalPosition(lCoordinate, tCoordinate, globalPosition);
    twoDSlidingFit.GetGlobalDirection(localGradient, globalDirection);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode GammaStartRefinementTool::FillHaloHitPositionVector(const CaloHitList &viewShowerHitList, const CaloHitList &showerSpineHitList, const CartesianVector &showerStartPosition, 
    const CartesianVector &showerStartDirection, const bool isEndDownstream, CartesianPointVector &haloHitPositionVector)
{
    // Find halo hits
    CaloHitList haloHitList;

    for (const CaloHit *const pCaloHit : viewShowerHitList)
    {
        if (std::find(showerSpineHitList.begin(), showerSpineHitList.end(), pCaloHit) != showerSpineHitList.end())
            continue;

        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        const float l(showerStartDirection.GetDotProduct(hitPosition - showerStartPosition));

        if ((isEndDownstream && (l < 0.f)) || (!isEndDownstream && (l > 0.f)))
            continue;

        if (showerStartDirection.GetCrossProduct(hitPosition - showerStartPosition).GetMagnitude() > m_molliereRadius)
            continue;

        haloHitPositionVector.push_back(pCaloHit->GetPositionVector());
        haloHitList.push_back(pCaloHit);
    }

    // Searching for a continuous shower, so identify any gaps
    CartesianPointVector coordinateListP, coordinateListN;

    try 
    {
        // Perform shower sliding fit
        const TwoDSlidingShowerFitResult twoDShowerSlidingFit(&haloHitPositionVector, m_showerSlidingFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const LayerFitResultMap &layerFitResultMapS(twoDShowerSlidingFit.GetShowerFitResult().GetLayerFitResultMap());
        const LayerFitResultMap &layerFitResultMapP(twoDShowerSlidingFit.GetPositiveEdgeFitResult().GetLayerFitResultMap());
        const LayerFitResultMap &layerFitResultMapN(twoDShowerSlidingFit.GetNegativeEdgeFitResult().GetLayerFitResultMap());

        float showerStartL(0.f), showerStartT(0.f);
        twoDShowerSlidingFit.GetShowerFitResult().GetLocalPosition(showerStartPosition, showerStartL, showerStartT);

        const int startLayer(twoDShowerSlidingFit.GetShowerFitResult().GetLayer(showerStartL));

        for (LayerFitResultMap::const_iterator iterS = layerFitResultMapS.begin(); iterS != layerFitResultMapS.end(); ++iterS)
        {
            const int layer(iterS->first);

            if (layer < startLayer)
                continue;

            LayerFitResultMap::const_iterator iterP = layerFitResultMapP.find(layer);
            LayerFitResultMap::const_iterator iterN = layerFitResultMapN.find(layer);

            if (layerFitResultMapP.end() != iterP)
            {
                CartesianVector positiveEdgePosition(0.f, 0.f, 0.f);
                twoDShowerSlidingFit.GetShowerFitResult().GetGlobalPosition(iterP->second.GetL(), iterP->second.GetFitT(), positiveEdgePosition);
                coordinateListP.push_back(positiveEdgePosition);
            }

            if (layerFitResultMapN.end() != iterN)
            {
                CartesianVector negativeEdgePosition(0.f, 0.f, 0.f);
                twoDShowerSlidingFit.GetShowerFitResult().GetGlobalPosition(iterN->second.GetL(), iterN->second.GetFitT(), negativeEdgePosition);
                coordinateListN.push_back(negativeEdgePosition);
            }
        }

        std::sort(coordinateListP.begin(), coordinateListP.end(), SortByDistanceToPoint(showerStartPosition));
        std::sort(coordinateListN.begin(), coordinateListN.end(), SortByDistanceToPoint(showerStartPosition));

        float pMaximumL(std::numeric_limits<float>::max());
        CartesianVector previousCoordinateP(coordinateListP.front());

        for (auto iterP = std::next(coordinateListP.begin()); iterP != coordinateListP.end(); ++iterP)
        {
            const CartesianVector &coordinateP(*iterP);
            const float separationSquared((coordinateP - previousCoordinateP).GetMagnitudeSquared());

            if (separationSquared > (5.f * 5.f))
            {
                std::cout << "gap in positive hits" << std::endl;
                break;
            }

            float thisT(0.f);
            twoDShowerSlidingFit.GetShowerFitResult().GetLocalPosition(coordinateP, pMaximumL, thisT);

            previousCoordinateP = coordinateP;
        }

        float nMaximumL(std::numeric_limits<float>::max());
        CartesianVector previousCoordinateN(coordinateListN.front());

        for (auto iterN = std::next(coordinateListN.begin()); iterN != coordinateListN.end(); ++iterN)
        {
            const CartesianVector &coordinateN(*iterN);
            const float separationSquared((coordinateN - previousCoordinateN).GetMagnitudeSquared());

            if (separationSquared > (5.f * 5.f))
            {
                std::cout << "gap in negative hits" << std::endl;
                break;
            }

            float thisT(0.f);
            twoDShowerSlidingFit.GetShowerFitResult().GetLocalPosition(coordinateN, nMaximumL, thisT);

            previousCoordinateN = coordinateN;
        }

        // Now refind the halo hit position vector
        haloHitPositionVector.clear();

        for (const CaloHit *const pCaloHit : haloHitList)
        {
            float thisL(0.f), thisT(0.f);
            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());

            twoDShowerSlidingFit.GetShowerFitResult().GetLocalPosition(hitPosition, thisL, thisT);

            if (thisL > std::max(pMaximumL, nMaximumL))
                continue;

            haloHitPositionVector.push_back(pCaloHit->GetPositionVector());
        }
    }
    catch (const StatusCodeException &)
    {
        std::cout << "couldn't perform first fit" << std::endl;
        return STATUS_CODE_FAILURE;
    }

    return STATUS_CODE_SUCCESS;
}



//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode GammaStartRefinementTool::CharacteriseShower(ShowerStartRefinementAlgorithm *const pAlgorithm, const CaloHitList &showerPfoHits, 
    const CaloHitList &showerSpineHits, const CartesianVector &showerStartPosition, const CartesianVector &showerStartDirection, const bool isEndDownstream, 
    CartesianVector &positiveEdgeStart, CartesianVector &positiveEdgeEnd, CartesianVector &negativeEdgeStart, CartesianVector &negativeEdgeEnd, bool &isBetween, bool &doesStraddle)
{
    CartesianPointVector haloHitPositionVector;

    if (this->FillHaloHitPositionVector(showerPfoHits, showerSpineHits, showerStartPosition, showerStartDirection, isEndDownstream, haloHitPositionVector) != STATUS_CODE_SUCCESS)
        return STATUS_CODE_FAILURE;

    //////////////////////////////////////
    /*
    CartesianVector end(showerStartPosition + (showerStartDirection * 10.f));
    PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &showerStartPosition, &end, "direction", VIOLET, 2, 1);

    for (auto &entry : haloHitPositionVector)
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &entry, "halo", RED, 2);

    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    */
    //////////////////////////////////////

    // now for reals...
    try
    {
        // Perform shower sliding fit
        const TwoDSlidingShowerFitResult twoDShowerSlidingFit(&haloHitPositionVector, m_showerSlidingFitWindow,  LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const LayerFitResultMap &layerFitResultMapS(twoDShowerSlidingFit.GetShowerFitResult().GetLayerFitResultMap());
        const LayerFitResultMap &layerFitResultMapP(twoDShowerSlidingFit.GetPositiveEdgeFitResult().GetLayerFitResultMap());
        const LayerFitResultMap &layerFitResultMapN(twoDShowerSlidingFit.GetNegativeEdgeFitResult().GetLayerFitResultMap());

        float showerStartL(0.f), showerStartT(0.f);
        twoDShowerSlidingFit.GetShowerFitResult().GetLocalPosition(showerStartPosition, showerStartL, showerStartT);

        const int startLayer(twoDShowerSlidingFit.GetShowerFitResult().GetLayer(showerStartL));
        const int endLayer(twoDShowerSlidingFit.GetShowerFitResult().GetMaxLayer());

        //////////////////////////////////////
        CartesianPointVector coordinateListP, coordinateListN;

        int layerCount(0);
        bool isFirstBetween(false), isLastBetween(false);

        for (LayerFitResultMap::const_iterator iterS = layerFitResultMapS.begin(); iterS != layerFitResultMapS.end(); ++iterS)
        {
            int layer(iterS->first);

            if (layer < startLayer)
                continue;

            LayerFitResultMap::const_iterator iterP = layerFitResultMapP.find(layer);
            LayerFitResultMap::const_iterator iterN = layerFitResultMapN.find(layer);

            CartesianVector positiveEdgePosition(0.f, 0.f, 0.f), negativeEdgePosition(0.f, 0.f, 0.f);
            if (layerFitResultMapP.end() != iterP)
            {
                twoDShowerSlidingFit.GetShowerFitResult().GetGlobalPosition(iterP->second.GetL(), iterP->second.GetFitT(), positiveEdgePosition);
                coordinateListP.push_back(positiveEdgePosition);
                //PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &positiveEdgePosition, "positive edge", GREEN, 2);
            }

            if (layerFitResultMapN.end() != iterN)
            {
                twoDShowerSlidingFit.GetShowerFitResult().GetGlobalPosition(iterN->second.GetL(), iterN->second.GetFitT(), negativeEdgePosition);
                coordinateListN.push_back(negativeEdgePosition);
                //PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &negativeEdgePosition, "negative edge", RED, 2);
            }

            if ((layerFitResultMapP.end() != iterP) && (layerFitResultMapN.end() != iterN))
            {
                const CartesianVector positiveDisplacement((positiveEdgePosition - showerStartPosition).GetUnitVector());
                const float positiveOpeningAngleFromCore(showerStartDirection.GetOpeningAngle(positiveDisplacement));

                const CartesianVector positiveClockwiseRotation(positiveDisplacement.GetZ() * std::sin(positiveOpeningAngleFromCore) + 
                    positiveDisplacement.GetX() * std::cos(positiveOpeningAngleFromCore), 0.f, positiveDisplacement.GetZ() * std::cos(positiveOpeningAngleFromCore) - 
                    positiveDisplacement.GetX() * std::sin(positiveOpeningAngleFromCore));

                const CartesianVector positiveAnticlockwiseRotation(positiveDisplacement.GetZ() * std::sin(-1.f * positiveOpeningAngleFromCore) + 
                    positiveDisplacement.GetX() * std::cos(-1.f * positiveOpeningAngleFromCore), 0.f, positiveDisplacement.GetZ() * std::cos(-1.f * positiveOpeningAngleFromCore) - 
                    positiveDisplacement.GetX() * std::sin(-1.f * positiveOpeningAngleFromCore));

                const float positiveClockwiseT((positiveClockwiseRotation - showerStartDirection).GetMagnitude());
                const float positiveAnticlockwiseT((positiveAnticlockwiseRotation - showerStartDirection).GetMagnitude());
                const bool isPositiveClockwise(positiveClockwiseT < positiveAnticlockwiseT);

                const CartesianVector negativeDisplacement((negativeEdgePosition - showerStartPosition).GetUnitVector());
                const float negativeOpeningAngleFromCore(showerStartDirection.GetOpeningAngle(negativeDisplacement));

                const CartesianVector negativeClockwiseRotation(negativeDisplacement.GetZ() * std::sin(negativeOpeningAngleFromCore) + 
                    negativeDisplacement.GetX() * std::cos(negativeOpeningAngleFromCore), 0.f, negativeDisplacement.GetZ() * std::cos(negativeOpeningAngleFromCore) - 
                    negativeDisplacement.GetX() * std::sin(negativeOpeningAngleFromCore));
                const CartesianVector negativeAnticlockwiseRotation(negativeDisplacement.GetZ() * std::sin(-1.f * negativeOpeningAngleFromCore) + 
                    negativeDisplacement.GetX() * std::cos(-1.f * negativeOpeningAngleFromCore), 0.f, negativeDisplacement.GetZ() * std::cos(-1.f * negativeOpeningAngleFromCore) - 
                    negativeDisplacement.GetX() * std::sin(-1.f * negativeOpeningAngleFromCore));

                const float negativeClockwiseT((negativeClockwiseRotation - showerStartDirection).GetMagnitude());
                const float negativeAnticlockwiseT((negativeAnticlockwiseRotation - showerStartDirection).GetMagnitude());
                const bool isNegativeClockwise(negativeClockwiseT < negativeAnticlockwiseT);

                ++layerCount;

                if (layerCount == 1)
                    isFirstBetween = (isPositiveClockwise != isNegativeClockwise);

                isLastBetween = (isPositiveClockwise != isNegativeClockwise);

                if (!doesStraddle)
                    doesStraddle = (isPositiveClockwise != isNegativeClockwise);
            }
        }

        isBetween = (isFirstBetween || isLastBetween);

        //////////////////////////
        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
        //////////////////////////

        // Find extremal coordinates

        // positive
        int positiveStartLayer(10000);
        int positiveEndLayer(10000);

        for (auto &entry : layerFitResultMapP)
        {
            const int bestStartSeparation(std::abs(startLayer - positiveStartLayer));
            const int thisStartSeparation(std::abs(startLayer - entry.first));

            if (((bestStartSeparation - thisStartSeparation) > 0) && (thisStartSeparation < bestStartSeparation))
                positiveStartLayer = entry.first;

            const int bestEndSeparation(std::abs(endLayer - positiveEndLayer));
            const int thisEndSeparation(std::abs(endLayer - entry.first));

            if (((bestEndSeparation - thisEndSeparation) > 0) && (thisEndSeparation < bestEndSeparation))
                positiveEndLayer = entry.first;
        }
        
        if (std::abs(startLayer - positiveStartLayer) > 5)
            return STATUS_CODE_FAILURE;
        
        // negative
        int negativeStartLayer(10000);
        int negativeEndLayer(10000);

        for (auto &entry : layerFitResultMapN)
        {
            const int bestStartSeparation(std::abs(startLayer - negativeStartLayer));
            const int thisStartSeparation(std::abs(startLayer - entry.first));

            if (((bestStartSeparation - thisStartSeparation) > 0) && (thisStartSeparation < bestStartSeparation))
                negativeStartLayer = entry.first;

            const int bestEndSeparation(std::abs(endLayer - negativeEndLayer));
            const int thisEndSeparation(std::abs(endLayer - entry.first));

            if (((bestEndSeparation - thisEndSeparation) > 0) && (thisEndSeparation < bestEndSeparation))
                negativeEndLayer = entry.first;
        }

        if (std::abs(startLayer - negativeStartLayer) > 5)
            return STATUS_CODE_FAILURE;

        const float showerStartPositiveLocalL(layerFitResultMapP.at(positiveStartLayer).GetL()), showerStartPositiveLocalT(layerFitResultMapP.at(positiveStartLayer).GetFitT());
        const float showerEndPositiveLocalL(layerFitResultMapP.at(positiveEndLayer).GetL()), showerEndPositiveLocalT(layerFitResultMapP.at(positiveEndLayer).GetFitT());
        const float showerStartNegativeLocalL(layerFitResultMapN.at(negativeStartLayer).GetL()), showerStartNegativeLocalT(layerFitResultMapN.at(negativeStartLayer).GetFitT());
        const float showerEndNegativeLocalL(layerFitResultMapN.at(negativeEndLayer).GetL()), showerEndNegativeLocalT(layerFitResultMapN.at(negativeEndLayer).GetFitT());

        twoDShowerSlidingFit.GetShowerFitResult().GetGlobalPosition(showerStartPositiveLocalL, showerStartPositiveLocalT, positiveEdgeStart);
        twoDShowerSlidingFit.GetShowerFitResult().GetGlobalPosition(showerEndPositiveLocalL, showerEndPositiveLocalT, positiveEdgeEnd);
        twoDShowerSlidingFit.GetShowerFitResult().GetGlobalPosition(showerStartNegativeLocalL, showerStartNegativeLocalT, negativeEdgeStart);
        twoDShowerSlidingFit.GetShowerFitResult().GetGlobalPosition(showerEndNegativeLocalL, showerEndNegativeLocalT, negativeEdgeEnd);

        /////////////////////////////////////////
        //PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &positiveEdgeStart, "Positive Edge", BLACK, 2);
        //PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &negativeEdgeStart, "Negative Edge", BLACK, 2);
        //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
        /////////////////////////////////////////

        // Swap start and end if needs be
        if ((showerStartPosition - positiveEdgeEnd).GetMagnitudeSquared() < (showerStartPosition - positiveEdgeStart).GetMagnitude())
        {
            const CartesianVector startTemp(positiveEdgeStart), endTemp(positiveEdgeEnd);

            positiveEdgeStart = endTemp;
            positiveEdgeEnd = startTemp;
        }

        if ((showerStartPosition - negativeEdgeEnd).GetMagnitudeSquared() < (showerStartPosition - negativeEdgeStart).GetMagnitude())
        {
            const CartesianVector startTemp(negativeEdgeStart), endTemp(negativeEdgeEnd);

            negativeEdgeStart = endTemp;
            negativeEdgeEnd = startTemp;
        }
    }
    catch (const StatusCodeException &)
    {
        std::cout << "couldn't perform second shower fit" << std::endl;
        return STATUS_CODE_FAILURE;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::RemoveConnectionPathway(const ProtoShower &protoShower)
{
    std::cout << protoShower.m_showerCore.m_startPosition.GetX() << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode GammaStartRefinementTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "Theta0XZBinSize", m_theta0XZBinSize));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "SmoothingWindow", m_smoothingWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "PathwaySearchRegion", m_pathwaySearchRegion));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MicroSlidingFitWindow", m_microSlidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinSigmaDeviation", m_minSigmaDeviation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "TrackSearchWindow", m_trackSearchWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "NInitialEnergyBins", m_nInitialEnergyBins));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinTrackBlipMean", m_minTrackBlipMean));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ShowerSlidingFitWindow", m_showerSlidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MolliereRadius", m_molliereRadius));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinShowerOpeningAngle", m_minShowerOpeningAngle));

    ShowerStartRefinementBaseTool::ReadSettings(xmlHandle);

    return STATUS_CODE_SUCCESS;
}


//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::FillTree(ShowerStartRefinementAlgorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo, const CartesianVector &nuVertexPosition)
{
    CaloHitList caloHits2D;
    LArPfoHelper::GetCaloHits(pShowerPfo, TPC_VIEW_U, caloHits2D);
    LArPfoHelper::GetCaloHits(pShowerPfo, TPC_VIEW_V, caloHits2D);
    LArPfoHelper::GetCaloHits(pShowerPfo, TPC_VIEW_W, caloHits2D);

    const CartesianVector nuVertexPositionU(LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertexPosition, TPC_VIEW_U));
    const CartesianVector nuVertexPositionV(LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertexPosition, TPC_VIEW_V));
    const CartesianVector nuVertexPositionW(LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertexPosition, TPC_VIEW_W));

    ++m_counter;

    FloatVector deviationAngleU, deviationAngleV, deviationAngleW;
    FloatVector xCoordinate, yCoordinate, zCoordinate;
    FloatVector lCoordinateU, tCoordinateU;
    FloatVector lCoordinateV, tCoordinateV;
    FloatVector lCoordinateW, tCoordinateW;
    FloatVector hitEnergyU, hitEnergyV, hitEnergyW;

    float lNuU(-999.f), tNuU(-999.f);
    float lNuV(-999.f), tNuV(-999.f);
    float lNuW(-999.f), tNuW(-999.f);

    CartesianPointVector initialHitPositionsU, initialHitPositionsV, initialHitPositionsW;

    for (const CaloHit *const pCaloHit : caloHits2D)
    {
        const HitType hitType(pCaloHit->GetHitType());
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        const CartesianVector nuVertexPositionProjection(hitType == TPC_VIEW_U ? nuVertexPositionU : hitType == TPC_VIEW_V ? nuVertexPositionV : nuVertexPositionW);
        CartesianPointVector &initialHitPositions(hitType == TPC_VIEW_U ? initialHitPositionsU : hitType == TPC_VIEW_V ? initialHitPositionsV : initialHitPositionsW);

        if ((hitPosition - nuVertexPositionProjection).GetMagnitude() > 14.f)
            continue;

        initialHitPositions.push_back(pCaloHit->GetPositionVector());
    }

    const TwoDSlidingFitResult twoDSlidingFitU(&initialHitPositionsU, 10000, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const TwoDSlidingFitResult twoDSlidingFitV(&initialHitPositionsV, 10000, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const TwoDSlidingFitResult twoDSlidingFitW(&initialHitPositionsW, 10000, LArGeometryHelper::GetWireZPitch(this->GetPandora()));

    twoDSlidingFitU.GetLocalPosition(nuVertexPosition, lNuU, tNuU);
    twoDSlidingFitV.GetLocalPosition(nuVertexPosition, lNuV, tNuV);
    twoDSlidingFitW.GetLocalPosition(nuVertexPosition, lNuW, tNuW);

    const CartesianVector xAxis(1.f, 0.f, 0.f);

    for (const CaloHit *const pCaloHit : caloHits2D)
    {
        const HitType hitType(pCaloHit->GetHitType());
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        const CartesianVector nuVertexPositionProjection(hitType == TPC_VIEW_U ? nuVertexPositionU : hitType == TPC_VIEW_V ? nuVertexPositionV : nuVertexPositionW);

        if ((hitPosition - nuVertexPositionProjection).GetMagnitude() > 14.f)
            continue;

        FloatVector &hitEnergy(hitType == TPC_VIEW_U ? hitEnergyU : hitType == TPC_VIEW_V ? hitEnergyV : hitEnergyW);

        hitEnergy.push_back(pCaloHit->GetElectromagneticEnergy());
        xCoordinate.push_back(hitPosition.GetX());
        yCoordinate.push_back(hitPosition.GetY());
        zCoordinate.push_back(hitPosition.GetZ());

        const TwoDSlidingFitResult &twoDSlidingFit(hitType == TPC_VIEW_U ? twoDSlidingFitU : hitType == TPC_VIEW_V ? twoDSlidingFitV : twoDSlidingFitW);
        FloatVector &lCoordinate(hitType == TPC_VIEW_U ? lCoordinateU : hitType == TPC_VIEW_V ? lCoordinateV : lCoordinateW);
        FloatVector &tCoordinate(hitType == TPC_VIEW_U ? tCoordinateU : hitType == TPC_VIEW_V ? tCoordinateV : tCoordinateW);

        float l, t;
        twoDSlidingFit.GetLocalPosition(hitPosition, l, t);

        lCoordinate.push_back(l);
        tCoordinate.push_back(t);

        FloatVector &deviationAngle(hitType == TPC_VIEW_U ? deviationAngleU : hitType == TPC_VIEW_V ? deviationAngleV : deviationAngleW);

        float deviationAngleTemp = xAxis.GetOpeningAngle(hitPosition - nuVertexPositionProjection);

        if (hitPosition.GetZ() < 0.f)
            deviationAngleTemp *= -1.f;

        deviationAngle.push_back(deviationAngleTemp);
    }

    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "ShowerCounter", m_counter));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "DeviationAngleU", &deviationAngleU));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "DeviationAngleV", &deviationAngleV));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "DeviationAngleW", &deviationAngleW));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "XCoordinate", &xCoordinate));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "YCoordinate", &yCoordinate));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "ZCoordinate", &zCoordinate));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "LCoordinateU", &lCoordinateU));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "LCoordinateV", &lCoordinateV));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "LCoordinateW", &lCoordinateW));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "TCoordinateU", &tCoordinateU));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "TCoordinateV", &tCoordinateV));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "TCoordinateW", &tCoordinateW));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "HitEnergyU", &hitEnergyU));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "HitEnergyV", &hitEnergyV));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "HitEnergyW", &hitEnergyW));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "LCoordinateNuU", lNuU));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "LCoordinateNuV", lNuV));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "LCoordinateNuW", lNuW));
    PANDORA_MONITORING_API(FillTree(pAlgorithm->GetPandora(), "ShowerDistribution"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::AngularDistributionTree(ShowerStartRefinementAlgorithm *const pAlgorithm, const AngularDecompositionMap &angularDecompositionMapU, 
    const AngularDecompositionMap &angularDecompositionMapV, const AngularDecompositionMap &angularDecompositionMapW)
{
    IntVector deviationBinU, deviationBinV, deviationBinW;
    FloatVector deviationWeightU, deviationWeightV, deviationWeightW;

    for (const auto &entry : angularDecompositionMapU)
    {
        deviationBinU.push_back(entry.first);
        deviationWeightU.push_back(entry.second);
    }

    for (const auto &entry : angularDecompositionMapV)
    {
        deviationBinV.push_back(entry.first);
        deviationWeightV.push_back(entry.second);
    }

    for (const auto &entry : angularDecompositionMapW)
    {
        deviationBinW.push_back(entry.first);
        deviationWeightW.push_back(entry.second);
    }

    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "DeviationAngle", "DeviationBinU", &deviationBinU));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "DeviationAngle", "DeviationBinV", &deviationBinV));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "DeviationAngle", "DeviationBinW", &deviationBinW));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "DeviationAngle", "DeviationWeightU", &deviationWeightU));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "DeviationAngle", "DeviationWeightV", &deviationWeightV));
    PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "DeviationAngle", "DeviationWeightW", &deviationWeightW));
    PANDORA_MONITORING_API(FillTree(pAlgorithm->GetPandora(), "DeviationAngle"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content
