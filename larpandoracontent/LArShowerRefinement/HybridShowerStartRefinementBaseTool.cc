/**
 *  @file   larpandoracontent/LArShowerRefinement/HybridShowerStartRefinementBaseTool.cc
 *
 *  @brief  Implementation of the shower start refinement base algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArHelpers/LArHitWidthHelper.h"
#include "larpandoracontent/LArHelpers/LArConnectionPathwayHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingShowerFitResult.h"

#include "larpandoracontent/LArShowerRefinement/LArProtoShower.h"
#include "larpandoracontent/LArShowerRefinement/HybridShowerStartRefinementBaseTool.h"

using namespace pandora;

namespace lar_content
{

HybridShowerStartRefinementBaseTool::HybridShowerStartRefinementBaseTool() : 
    m_maxDistanceForConnection(3.f),
    m_pathwaySearchRegion(14.f),
    m_theta0XZBinSize(0.005f),
    m_smoothingWindow(3),
    m_growingFitInitialLength(10.f),
    m_macroSlidingFitWindow(1000),
    m_growingFitSegmentLength(5.0f),
    m_distanceToLine(1.f),
    m_initialFitDistanceToLine(0.5f),
    m_maxFittingHits(20),
    m_longitudinalCoordinateBinSize(1.f),
    m_hitConnectionDistance(1.f),
    m_minInitialHitsFound(10),
    m_microSlidingFitWindow(20),
    m_nInitialEnergyBins(5),
    m_minSigmaDeviation(5.f),
    m_molliereRadius(9.f),
    m_showerSlidingFitWindow(20),
    m_minShowerOpeningAngle(3.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool HybridShowerStartRefinementBaseTool::HasPathToNuVertex(const ParticleFlowObject *const pShowerPfo, const CartesianVector &neutrinoVertex) const
{
    if (this->HasPathToNuVertex(pShowerPfo, neutrinoVertex, TPC_VIEW_U))
        return true;

    if (this->HasPathToNuVertex(pShowerPfo, neutrinoVertex, TPC_VIEW_V))
        return true;

    if (this->HasPathToNuVertex(pShowerPfo, neutrinoVertex, TPC_VIEW_W))
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool HybridShowerStartRefinementBaseTool::HasPathToNuVertex(const ParticleFlowObject *const pShowerPfo, const CartesianVector &neutrinoVertex, const HitType &hitType) const
{
    // this algorithm has to run after we have 3D hits but not after shower vertex creation
    // therefore, find the closest space point to the neutrino vertex

    ClusterList clusterList;
    LArPfoHelper::GetClusters(pShowerPfo, hitType, clusterList);

    const CartesianVector projectedNuVertexPosition(LArGeometryHelper::ProjectPosition(this->GetPandora(), neutrinoVertex, hitType));

    for (const Cluster *const pCluster : clusterList)
    {
        CaloHitList caloHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

        for (const CaloHit *const caloHit : caloHitList)
        {
            const CartesianVector &hitPosition(caloHit->GetPositionVector());
            const float separationSquared((projectedNuVertexPosition - hitPosition).GetMagnitudeSquared());

            if (separationSquared < (m_maxDistanceForConnection * m_maxDistanceForConnection))
                return true;
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

/*
bool HybridShowerStartRefinementBaseTool::HasPathToNuVertex(const ParticleFlowObject *const pShowerPfo, const CartesianVector &neutrinoVertex) const
{
    // this algorithm has to run after we have 3D hits but not after shower vertex creation
    // therefore, find the closest space point to the neutrino vertex

    ClusterList clusterList3D;
    LArPfoHelper::GetClusters(pShowerPfo, TPC_3D, clusterList3D);

    for (const Cluster *const pCluster3D : clusterList3D)
    {
        CaloHitList caloHitList3D;
        pCluster3D->GetOrderedCaloHitList().FillCaloHitList(caloHitList3D);

        for (const CaloHit *const caloHit3D : caloHitList3D)
        {
            const CartesianVector &hitPosition(caloHit3D->GetPositionVector());
            const float separationSquared((neutrinoVertex - hitPosition).GetMagnitudeSquared());

            if (separationSquared < (m_maxDistanceForConnection * m_maxDistanceForConnection))
                return true;
        }
    }

    return false;
}
*/

//------------------------------------------------------------------------------------------------------------------------------------------

void HybridShowerStartRefinementBaseTool::FillAngularDecompositionMap(const CaloHitList &viewShowerHitList, const CartesianVector &projectedNuVertexPosition, 
    AngularDecompositionMap &angularDecompositionMap)
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

void HybridShowerStartRefinementBaseTool::SmoothAngularDecompositionMap(AngularDecompositionMap &angularDecompositionMap)
{
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
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HybridShowerStartRefinementBaseTool::ObtainPeakVector(AngularDecompositionMap &angularDecompositionMap, IntVector &viewPeakVector)
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
            throw StatusCodeException(STATUS_CODE_FAILURE);
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
            throw StatusCodeException(STATUS_CODE_FAILURE);
        }

        const float precedingBinWeight(angularDecompositionMap.find(precedingBin) == angularDecompositionMap.end() ? 0.f : angularDecompositionMap.at(precedingBin));
        const float followingBinWeight(angularDecompositionMap.find(followingBin) == angularDecompositionMap.end() ? 0.f : angularDecompositionMap.at(followingBin));

        if ((binWeight < precedingBinWeight) || (binWeight < followingBinWeight))
            continue;

        viewPeakVector.push_back(bin);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool HybridShowerStartRefinementBaseTool::FindBestAngularPeak(AngularDecompositionMap &angularDecompositionMap, IntVector &viewPeakVector, 
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
            throw StatusCodeException(STATUS_CODE_FAILURE);
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

void HybridShowerStartRefinementBaseTool::FindShowerSpine(const HybridShowerStartRefinementAlgorithm *pAlgorithm, const CaloHitList &viewShowerHitList, 
    const CartesianVector &projectedNuVertexPosition, const CartesianVector &initialDirection, CaloHitList &unavailableHitList, CaloHitList &showerSpineHitList)
{
    // Construct initial fit
    float highestL(0.f);
    CartesianPointVector runningFitPositionVector;

    for (const CaloHit *const pCaloHit : viewShowerHitList)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        const CartesianVector &displacementVector(hitPosition - projectedNuVertexPosition);
        const float l(initialDirection.GetDotProduct(displacementVector));
        const float t(initialDirection.GetCrossProduct(displacementVector).GetMagnitude());

        if ((l < m_growingFitInitialLength) && (l > 0.f) && (t < m_initialFitDistanceToLine))
        {
            if (l > highestL)
                highestL = l;

            if (std::find(unavailableHitList.begin(), unavailableHitList.end(), pCaloHit) == unavailableHitList.end())
            {
                showerSpineHitList.push_back(pCaloHit);
            }

            ////////////////////////////
            //PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "INITIAL HIT", VIOLET, 2);
            ////////////////////////////
            runningFitPositionVector.push_back(hitPosition);
        }
    }

    // Require significant number of initial hits
    if (runningFitPositionVector.size() < m_minInitialHitsFound)
    {
        ////////////////////////////////////
        //std::cout << "Not enough initial hits to form fit" << std::endl;
        //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
        ////////////////////////////////////

        showerSpineHitList.clear();
        return;
    }

    ///////////////////
    //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    ////////////////////

    // Perform a running fit to follow pathway
    const bool isEndDownstream(initialDirection.GetZ() > 0.f);
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    CartesianVector extrapolatedStartPosition(projectedNuVertexPosition), extrapolatedDirection(initialDirection);
    CartesianVector extrapolatedEndPosition(extrapolatedStartPosition + (extrapolatedDirection * highestL));

    unsigned int count(0);
    bool hitsCollected(true);

    // make this a boolean
    while (hitsCollected)
    {
        ++count;

        try
        {
            const int excessHitsInFit(runningFitPositionVector.size() - m_maxFittingHits);

            if (excessHitsInFit > 0)
            {
                // Remove furthest away hits
                std::sort(runningFitPositionVector.begin(), runningFitPositionVector.end(), LArConnectionPathwayHelper::SortByDistanceToPoint(extrapolatedEndPosition));

                for (int i = 0; i < excessHitsInFit; ++i)
                    runningFitPositionVector.erase(std::prev(runningFitPositionVector.end()));
            }

            const TwoDSlidingFitResult extrapolatedFit(&runningFitPositionVector, m_macroSlidingFitWindow, slidingFitPitch);            

            extrapolatedStartPosition = count == 1 ? extrapolatedEndPosition : isEndDownstream ? extrapolatedFit.GetGlobalMaxLayerPosition() : extrapolatedFit.GetGlobalMinLayerPosition();
            extrapolatedDirection = isEndDownstream ? extrapolatedFit.GetGlobalMaxLayerDirection() : extrapolatedFit.GetGlobalMinLayerDirection() * (-1.f);
            extrapolatedEndPosition = extrapolatedStartPosition + (extrapolatedDirection * m_growingFitSegmentLength);

            hitsCollected = this->CollectSubsectionHits(pAlgorithm, extrapolatedFit, extrapolatedStartPosition, extrapolatedEndPosition, extrapolatedDirection,
                isEndDownstream, viewShowerHitList, runningFitPositionVector, unavailableHitList, showerSpineHitList);

            // As a final effort, reduce the sliding fit window
            if (!hitsCollected)
            {
                const TwoDSlidingFitResult microExtrapolatedFit(&runningFitPositionVector, 5.f, slidingFitPitch);
            
                extrapolatedStartPosition = count == 1 ? extrapolatedStartPosition : isEndDownstream ? microExtrapolatedFit.GetGlobalMaxLayerPosition() : 
                    microExtrapolatedFit.GetGlobalMinLayerPosition();
                extrapolatedDirection = isEndDownstream ? microExtrapolatedFit.GetGlobalMaxLayerDirection() : microExtrapolatedFit.GetGlobalMinLayerDirection() * (-1.f);
                extrapolatedEndPosition = extrapolatedStartPosition + (extrapolatedDirection * m_growingFitSegmentLength);

                hitsCollected = this->CollectSubsectionHits(pAlgorithm, extrapolatedFit, extrapolatedStartPosition, extrapolatedEndPosition, extrapolatedDirection,
                    isEndDownstream, viewShowerHitList, runningFitPositionVector, unavailableHitList, showerSpineHitList);
            }
        }
        catch (const StatusCodeException &)
        {
            ////////////////////////////////////
            //std::cout << "couldn't fit runningFitPositionVector" << std::endl;
            ////////////////////////////////////

            return;
        }
    }

    ////////////////////////////////////
    //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    ////////////////////////////////////
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool HybridShowerStartRefinementBaseTool::CollectSubsectionHits(const HybridShowerStartRefinementAlgorithm *pAlgorithm, const TwoDSlidingFitResult &extrapolatedFit, 
    const CartesianVector &extrapolatedStartPosition, const CartesianVector &extrapolatedEndPosition, const CartesianVector &extrapolatedDirection, 
    const bool isEndDownstream, const CaloHitList &viewShowerHitList, CartesianPointVector &runningFitPositionVector, CaloHitList &unavailableHitList, 
    CaloHitList &showerSpineHitList)
{ 
    ////////////////////////////
    //PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &extrapolatedStartPosition, "start", GREEN, 2);
    //PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &extrapolatedEndPosition, "end", GREEN, 2);
    ////////////////////////////

    float extrapolatedStartL(0.f), extrapolatedStartT(0.f);
    extrapolatedFit.GetLocalPosition(extrapolatedStartPosition, extrapolatedStartL, extrapolatedStartT);

    float extrapolatedEndL(0.f), extrapolatedEndT(0.f);
    extrapolatedFit.GetLocalPosition(extrapolatedEndPosition, extrapolatedEndL, extrapolatedEndT);

    CaloHitList collectedHits;

    for (const CaloHit *const pCaloHit : viewShowerHitList)
    {
        if (std::find(showerSpineHitList.begin(), showerSpineHitList.end(), pCaloHit) != showerSpineHitList.end())
            continue;

        if (std::find(unavailableHitList.begin(), unavailableHitList.end(), pCaloHit) != unavailableHitList.end())
            continue;

        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());

        float hitL(0.f), hitT(0.f);
        extrapolatedFit.GetLocalPosition(hitPosition, hitL, hitT);

        // Assess whether hit is within section boundaries
        if (isEndDownstream && ((hitL < extrapolatedStartL) || (hitL > extrapolatedEndL)))
            continue;

        if (!isEndDownstream && ((hitL > extrapolatedStartL) || (hitL < extrapolatedEndL)))
            continue;

        // Assess whether hit is close to connecting line - taking account hit width if necessary
        if (this->IsCloseToLine(hitPosition, extrapolatedStartPosition, extrapolatedDirection, m_distanceToLine))
        {
            collectedHits.push_back(pCaloHit);
        }
        else
        {
            const CartesianVector closestPointInHit(LArHitWidthHelper::GetClosestPointToLine2D(extrapolatedStartPosition, extrapolatedDirection, pCaloHit));

            if (this->IsCloseToLine(closestPointInHit, extrapolatedStartPosition, extrapolatedDirection, m_distanceToLine))
                collectedHits.push_back(pCaloHit);
        }
    }

    const int nInitialHits(showerSpineHitList.size());
    this->CollectConnectedHits(pAlgorithm, collectedHits, extrapolatedStartPosition, extrapolatedDirection, runningFitPositionVector, unavailableHitList, showerSpineHitList);
    const int nFinalHits(showerSpineHitList.size());

    ////////////////////////////////////
    //PandoraMonitoringApi::ViewEvent(this->GetPandora());
    ////////////////////////////////////

    return (nFinalHits != nInitialHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HybridShowerStartRefinementBaseTool::CollectConnectedHits(const HybridShowerStartRefinementAlgorithm *pAlgorithm, const CaloHitList &collectedHits, 
    const CartesianVector &extrapolatedStartPosition, const CartesianVector &extrapolatedDirection, CartesianPointVector &runningFitPositionVector, 
    CaloHitList &/*unavailableHitList*/, CaloHitList &showerSpineHitList)
{
    // Now add connected hits - taking into account hit width if necessary
    bool found = true;

    while (found)
    {
        found = false;

        for (const CaloHit *const pCaloHit : collectedHits)
        {
            if (std::find(showerSpineHitList.begin(), showerSpineHitList.end(), pCaloHit) != showerSpineHitList.end())
                continue;

            CartesianVector hitPosition(pCaloHit->GetPositionVector());

            if (this->GetClosestDistance(hitPosition, runningFitPositionVector) > m_hitConnectionDistance)
            {
                hitPosition = LArHitWidthHelper::GetClosestPointToLine2D(extrapolatedStartPosition, extrapolatedDirection, pCaloHit);
                        
                if (LArHitWidthHelper::GetClosestDistance(pCaloHit, showerSpineHitList) > m_hitConnectionDistance)
                    continue;
            }

            found = true;

            runningFitPositionVector.push_back(hitPosition);
            showerSpineHitList.push_back(pCaloHit);
 
            ////////////////////////////////
            //PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "added hit", BLUE, 2);
            ////////////////////////////////
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool HybridShowerStartRefinementBaseTool::IsCloseToLine(const CartesianVector &hitPosition, const CartesianVector &lineStart, 
    const CartesianVector &lineDirection, const float distanceToLine) const
{
    const float transverseDistanceFromLine(lineDirection.GetCrossProduct(hitPosition - lineStart).GetMagnitude());
    
    if (transverseDistanceFromLine > distanceToLine)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float HybridShowerStartRefinementBaseTool::GetClosestDistance(const CartesianVector &position, const CartesianPointVector &testPositions) const
{
    float closestDistanceSqaured(std::numeric_limits<float>::max());

    for (const CartesianVector &testPosition : testPositions)
    {
        const float separationSquared((testPosition - position).GetMagnitudeSquared());

        if (separationSquared < closestDistanceSqaured)
            closestDistanceSqaured = separationSquared;
    }

    return std::sqrt(closestDistanceSqaured);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HybridShowerStartRefinementBaseTool::ObtainLongitudinalDecomposition(HybridShowerStartRefinementAlgorithm *const pAlgorithm, const CaloHitList &showerSpineHitList, 
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

        /////////////////////////////////
        /*
        CartesianVector blob(0.f, 0.f, 0.f);
        twoDSlidingFit.GetGlobalPosition(layerFitResultMap.at(layer).GetL(), layerFitResultMap.at(layer).GetFitT(), blob);
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &blob, std::to_string(layer), BLUE, 2);
        */
        /////////////////////////////////

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

    /////////////////////////////////
    //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    /////////////////////////////////
}


//------------------------------------------------------------------------------------------------------------------------------------------

void HybridShowerStartRefinementBaseTool::GetEnergyDistribution(HybridShowerStartRefinementAlgorithm *const pAlgorithm, const CaloHitList &showerSpineHitList, 
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

bool HybridShowerStartRefinementBaseTool::FindShowerStart(HybridShowerStartRefinementAlgorithm *const pAlgorithm, const CartesianVector &projectedNuVertexPosition, 
    const CartesianVector &startDirection, const LongitudinalPositionMap &longitudinalPositionMap, const EnergySpectrumMap &energySpectrumMap, 
    const CaloHitList &showerSpineHitList, CartesianVector &showerStartPosition, const CaloHitList &showerPfoHitList, const bool isEndDownstream, 
    ProtoShowerVector &protoShowerVector, const bool isHelper)
{
    float meanEnergy(0.f), energySigma(0.f);

    if (energySpectrumMap.size() > m_nInitialEnergyBins)
        this->CharacteriseInitialEnergy(energySpectrumMap, isEndDownstream, meanEnergy, energySigma);

    auto energySpectrumIter(isEndDownstream ? energySpectrumMap.begin() : std::prev(energySpectrumMap.end()));

    bool showerStartFound(false);
    for (unsigned int i = 0; i < energySpectrumMap.size(); ++i)
    {
        // Bypass bins used to asses the initial energy
        if (i < m_nInitialEnergyBins)
        {
            isEndDownstream ? ++energySpectrumIter : --energySpectrumIter;
            continue;
        }

        const float longitudinalCoordinate(energySpectrumIter->first * m_longitudinalCoordinateBinSize);
        const float energyDeviation((energySpectrumIter->second - meanEnergy) / energySigma);

        /////////////////////////////////
        /*        
        std::cout << "energyDeviation (sigma): " << energyDeviation << std::endl;

        CartesianVector jam(0.f, 0.f, 0.f), jam2(0.f, 0.f, 0.f);
        this->ConvertLongitudinalProjectionToGlobalPosition(pAlgorithm, showerSpineHitList, longitudinalCoordinate, jam, jam2);
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &jam, std::to_string(longitudinalCoordinate), GREEN, 2);
        */
        /////////////////////////////////

        // Use energy and local topology to assess whether we are at the shower start
        if ((energyDeviation > m_minSigmaDeviation) && this->IsShowerTopology(pAlgorithm, longitudinalCoordinate,
            showerPfoHitList, showerSpineHitList, isEndDownstream))
        {
            showerStartFound = true;
            break;
        }
        
        if (i != (energySpectrumMap.size() - 1))
            isEndDownstream ? ++energySpectrumIter : --energySpectrumIter;
    }

    const int longitudinalStartBin(energySpectrumIter->first);
    const float longitudinalStartCoordinate(longitudinalStartBin * m_longitudinalCoordinateBinSize);

    CartesianVector showerStartDirection(0.f, 0.f, 0.f);
    this->ConvertLongitudinalProjectionToGlobalPosition(pAlgorithm, showerSpineHitList, longitudinalStartCoordinate, showerStartPosition, showerStartDirection);

    showerStartDirection = isEndDownstream ? showerStartDirection : showerStartDirection * (-1.f);

    CaloHitList pathwayHitList, showerHitList;
    for (const CaloHit *const pCaloHit : showerPfoHitList)
    {
        if (longitudinalPositionMap.find(pCaloHit) != longitudinalPositionMap.end())
        {
            if ((isEndDownstream && (longitudinalPositionMap.at(pCaloHit) < longitudinalStartCoordinate)) || 
                (!isEndDownstream && (longitudinalPositionMap.at(pCaloHit) > longitudinalStartCoordinate)))
            {
                pathwayHitList.push_back(pCaloHit);
            }
            else
            {
                showerHitList.push_back(pCaloHit);
            }
        }
        else
        {
            const float t(showerStartDirection.GetCrossProduct(showerStartPosition - pCaloHit->GetPositionVector()).GetMagnitude());
            const float l(showerStartDirection.GetDotProduct(pCaloHit->GetPositionVector() - showerStartPosition));

            if ((l > 0.f) && (t < m_molliereRadius))
                showerHitList.push_back(pCaloHit);
        }
    }
    
    showerHitList.clear();


    /////////////////////////////////
    /*
    for (const CaloHit *const pCaloHit : pathwayHitList)
    {
        const CartesianVector &hit(pCaloHit->GetPositionVector());
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hit, "path", BLUE, 2);
    }
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());    
    */
    /////////////////////////////////

    protoShowerVector.push_back(ProtoShower(ShowerCore(showerStartPosition, showerStartDirection, showerHitList), 
        ConnectionPathway(projectedNuVertexPosition, startDirection, pathwayHitList), showerSpineHitList, isHelper, CaloHitList(), CartesianPointVector()));

    return showerStartFound;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HybridShowerStartRefinementBaseTool::CharacteriseInitialEnergy(const EnergySpectrumMap &energySpectrumMap, const bool isEndDownstream, float &meanEnergy, float &energySigma)
{
    auto energySpectrumIter(isEndDownstream ? energySpectrumMap.begin() : std::prev(energySpectrumMap.end()));

    for (unsigned int i = 0 ; i < m_nInitialEnergyBins; ++i)
    {
        meanEnergy += energySpectrumIter->second;
        isEndDownstream ? ++energySpectrumIter : --energySpectrumIter; 
    }

    meanEnergy /= static_cast<float>(m_nInitialEnergyBins);

    energySpectrumIter = isEndDownstream ? energySpectrumMap.begin() : std::prev(energySpectrumMap.end());

    for (unsigned int i = 0; i < m_nInitialEnergyBins; ++i)
    {
        energySigma += std::pow(energySpectrumIter->second - meanEnergy, 2);
        isEndDownstream ? ++energySpectrumIter : --energySpectrumIter; 
    }

    energySigma = std::sqrt(energySigma / m_nInitialEnergyBins);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool HybridShowerStartRefinementBaseTool::IsShowerTopology(HybridShowerStartRefinementAlgorithm *const pAlgorithm, const float longitudinalDistance, 
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
        /////////////////////////////////
        //std::cout << "failed to characterise the shower" << std::endl;
        /////////////////////////////////
        return false;
    }

    if (!isBetween)
    {
        /////////////////////////////////
        //std::cout << "shower start is not inbetween start edge points :( " << std::endl;
        /////////////////////////////////
        return false;
    }

    if (!doesStraddle)
    {
        /////////////////////////////////
        //std::cout << "shower does not straddle core " << std::endl;
        /////////////////////////////////
        return false;
    }

    positiveEdgeStart = showerStartPosition;
    negativeEdgeStart = showerStartPosition;
    positiveEdgeDirection = positiveEdgeEnd - positiveEdgeStart;
    negativeEdgeDirection = negativeEdgeEnd - negativeEdgeStart;

    const float showerOpeningAngle(positiveEdgeDirection.GetOpeningAngle(negativeEdgeDirection) * 180.f / 3.14);
    const float positiveEdgeDeviation(positiveEdgeDirection.GetOpeningAngle(showerStartDirection) * 180.f / 3.14);
    const float negativeEdgeDeviation(negativeEdgeDirection.GetOpeningAngle(showerStartDirection) * 180.f / 3.14);

    /////////////////////////////////
    /*
    std::cout << "edge opening angle: " << (positiveEdgeDirection.GetOpeningAngle(negativeEdgeDirection) * 180.f / 3.14) << std::endl;
    std::cout << "positive opening angle: " << (positiveEdgeDirection.GetOpeningAngle(showerStartDirection) * 180.f / 3.14) << std::endl;
    std::cout << "negative opening angle: " << (negativeEdgeDirection.GetOpeningAngle(showerStartDirection) * 180.f / 3.14) << std::endl;

    PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &positiveEdgeStart, &positiveEdgeEnd, "positive direction", BLACK, 2, 1);
    PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &negativeEdgeStart, &negativeEdgeEnd, "negative direction", BLACK, 2, 1);
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    */
    /////////////////////////////////

    if (showerOpeningAngle < m_minShowerOpeningAngle)
    {
        /////////////////////////////////
        //std::cout << "shower opening angle is too small" << std::endl;
        /////////////////////////////////
        return false;
    }

    if (positiveEdgeDeviation > 45.f)
    {
        /////////////////////////////////
        //std::cout << "positive edge is not in the direction of the shower" << std::endl;
        /////////////////////////////////
    }

    if (negativeEdgeDeviation > 45.f)
    {
        /////////////////////////////////
        //std::cout << "negative edge is not in the direction of the shower" << std::endl;
        /////////////////////////////////
    }
    
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HybridShowerStartRefinementBaseTool::ConvertLongitudinalProjectionToGlobalPosition(HybridShowerStartRefinementAlgorithm *const pAlgorithm, 
    const CaloHitList &showerSpineHitList, const float longitudinalDistance, CartesianVector &globalPosition, CartesianVector &globalDirection)
{
    CartesianPointVector hitPositions;

    for (const CaloHit *const pCaloHit : showerSpineHitList)
        hitPositions.push_back(pCaloHit->GetPositionVector());

    const TwoDSlidingFitResult twoDSlidingFit(&hitPositions, m_microSlidingFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const LayerFitResultMap &layerFitResultMap(twoDSlidingFit.GetLayerFitResultMap());

    float runningDistance(0.f);

    CartesianVector previousLayerPosition(0.f, 0.f, 0.f);
    twoDSlidingFit.GetGlobalPosition(layerFitResultMap.begin()->second.GetL(), layerFitResultMap.begin()->second.GetFitT(), previousLayerPosition);

    /////////////////////////////////
    /*
    const CartesianVector jam(twoDSlidingFit.GetGlobalMinLayerPosition());
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &jam, "MIN LAYER POSITION", RED, 2);
    */
    /////////////////////////////////

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

StatusCode HybridShowerStartRefinementBaseTool::CharacteriseShower(HybridShowerStartRefinementAlgorithm *const pAlgorithm, const CaloHitList &showerPfoHits, 
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
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &entry, "halo", BLACK, 2);
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

        CartesianPointVector coordinateListP, coordinateListN;

        int layerCount(0);
        bool isFirstBetween(false), isLastBetween(false);

        for (LayerFitResultMap::const_iterator iterS = layerFitResultMapS.begin(); iterS != layerFitResultMapS.end(); ++iterS)
        {
            int layer(iterS->first);

            if (isEndDownstream && (layer < startLayer))
                continue;

            if (!isEndDownstream && (layer > startLayer))
                continue;

            LayerFitResultMap::const_iterator iterP = layerFitResultMapP.find(layer);
            LayerFitResultMap::const_iterator iterN = layerFitResultMapN.find(layer);

            CartesianVector positiveEdgePosition(0.f, 0.f, 0.f), negativeEdgePosition(0.f, 0.f, 0.f);
            if (layerFitResultMapP.end() != iterP)
            {
                twoDShowerSlidingFit.GetShowerFitResult().GetGlobalPosition(iterP->second.GetL(), iterP->second.GetFitT(), positiveEdgePosition);
                coordinateListP.push_back(positiveEdgePosition);

                /////////////////////////////////
                //PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &positiveEdgePosition, "positive edge", GREEN, 2);
                /////////////////////////////////
            }

            if (layerFitResultMapN.end() != iterN)
            {
                twoDShowerSlidingFit.GetShowerFitResult().GetGlobalPosition(iterN->second.GetL(), iterN->second.GetFitT(), negativeEdgePosition);
                coordinateListN.push_back(negativeEdgePosition);

                /////////////////////////////////
                //PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &negativeEdgePosition, "negative edge", RED, 2);
                /////////////////////////////////
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
        /*
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &positiveEdgeStart, "Positive Edge", BLACK, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &negativeEdgeStart, "Negative Edge", BLACK, 2);
        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
        */
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
        /////////////////////////////////
        //std::cout << "couldn't perform second shower fit" << std::endl;
        /////////////////////////////////

        return STATUS_CODE_FAILURE;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HybridShowerStartRefinementBaseTool::FillHaloHitPositionVector(const CaloHitList &viewShowerHitList, const CaloHitList &showerSpineHitList, const CartesianVector &showerStartPosition, 
    const CartesianVector &showerStartDirection, const bool isEndDownstream, CartesianPointVector &haloHitPositionVector)
{
    // Find halo hits
    CaloHitList haloHitList;

    for (const CaloHit *const pCaloHit : viewShowerHitList)
    {
        if (std::find(showerSpineHitList.begin(), showerSpineHitList.end(), pCaloHit) != showerSpineHitList.end())
            continue;

        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        float l(showerStartDirection.GetDotProduct(hitPosition - showerStartPosition));

        l *= isEndDownstream ? 1.f : -1.f;

        if (l < 0.f)
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
      
            if (isEndDownstream && (layer < startLayer))
                continue;

            if (!isEndDownstream && (layer > startLayer))
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

        if ((coordinateListP.size() == 0) || (coordinateListN.size() == 0))
            return STATUS_CODE_FAILURE;

        std::sort(coordinateListP.begin(), coordinateListP.end(), LArConnectionPathwayHelper::SortByDistanceToPoint(showerStartPosition));
        std::sort(coordinateListN.begin(), coordinateListN.end(), LArConnectionPathwayHelper::SortByDistanceToPoint(showerStartPosition));

        float pMaximumL(std::numeric_limits<float>::max());
        CartesianVector previousCoordinateP(coordinateListP.front());

        for (auto iterP = std::next(coordinateListP.begin()); iterP != coordinateListP.end(); ++iterP)
        {
            const CartesianVector &coordinateP(*iterP);
            const float separationSquared((coordinateP - previousCoordinateP).GetMagnitudeSquared());

            if (separationSquared > (3.f * 3.f))
            {
                /////////////////////////////////
                //std::cout << "gap in positive hits" << std::endl;
                /////////////////////////////////

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

            if (separationSquared > (3.f * 3.f))
            {
                /////////////////////////////////
                //std::cout << "gap in negative hits" << std::endl;
                /////////////////////////////////
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

            // 0.85f because the end of showers can get really narrow and thus fail the opening angle check
            if (isEndDownstream && (thisL > (0.85 * std::max(pMaximumL, nMaximumL))))
                continue;

            if (!isEndDownstream && (thisL < (0.85 * std::min(pMaximumL, nMaximumL))))
                continue;

            haloHitPositionVector.push_back(pCaloHit->GetPositionVector());
        }

        /////////////////////////////////
        /*
        std::cout << "coordinateListP.size(): " << coordinateListP.size() << std::endl;
        std::cout << "coordinateListN.size(): " << coordinateListP.size() << std::endl;
        std::cout << "haloHitList.size(): " << haloHitList.size() << std::endl;
        */
        /////////////////////////////////
    }
    catch (const StatusCodeException &)
    {
        /////////////////////////////////
        std::cout << "couldn't perform first fit" << std::endl;
        /////////////////////////////////

        return STATUS_CODE_FAILURE;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------------------------------------


StatusCode HybridShowerStartRefinementBaseTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxDistanceForConnection", m_maxDistanceForConnection));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "PathwaySearchRegion", m_pathwaySearchRegion));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "Theta0XZBinSize", m_theta0XZBinSize));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "SmoothingWindow", m_smoothingWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "GrowingFitInitialLength", m_growingFitInitialLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MacroSlidingFitWindow", m_macroSlidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "GrowingFitSegmentLength", m_growingFitSegmentLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "DistanceToLine", m_distanceToLine));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "InitialFitDistanceToLine", m_initialFitDistanceToLine));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxFittingsHits", m_maxFittingHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "LongitudinalCoordinateBinSize", m_longitudinalCoordinateBinSize));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "HitConnectionDistance", m_hitConnectionDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "HitInitialHitsFound", m_minInitialHitsFound));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MicroSlidingFitWindow", m_microSlidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "NInitialEnergyBins", m_nInitialEnergyBins));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinSigmaDeviation", m_minSigmaDeviation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MolliereRadius", m_molliereRadius));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ShowerSlidingFitWindow", m_showerSlidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinShowerOpeningAngle", m_minShowerOpeningAngle));

    return STATUS_CODE_SUCCESS;
}


} // namespace lar_content
