/**
 *  @file   larpandoracontent/LArShowerRefinement/ShowerStartFinderTool.cc
 *
 *  @brief  Implementation of the shower start finder tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingShowerFitResult.h"

#include "larpandoracontent/LArHelpers/LArConnectionPathwayHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArShowerRefinement/ShowerStartFinderTool.h"

using namespace pandora;

namespace lar_content
{

ShowerStartFinderTool::ShowerStartFinderTool() :
    m_spineSlidingFitWindow(10),
    m_longitudinalCoordinateBinSize(1.f),
    m_nInitialEnergyBins(5),
    m_minSigmaDeviation(1.f),
    m_maxEdgeGap(3.f),
    m_longitudinalShowerFraction(0.85f),
    m_minShowerOpeningAngle(2.f),
    m_molliereRadius(4.5f),
    m_showerSlidingFitWindow(1000),
    m_maxLayerSeparation(5)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerStartFinderTool::Run(const ParticleFlowObject *const pShowerPfo, const CartesianVector &peakDirection,
    const HitType hitType, const CaloHitList &showerSpineHitList, CartesianVector &showerStartPosition, CartesianVector &showerStartDirection)
{
    try
    {
        // Fit the shower spine
        CartesianPointVector hitPositions;
        for (const CaloHit *const pCaloHit : showerSpineHitList)
            hitPositions.push_back(pCaloHit->GetPositionVector());

        const TwoDSlidingFitResult spineTwoDSlidingFit(
            &hitPositions, m_spineSlidingFitWindow, LArGeometryHelper::GetWirePitch(this->GetPandora(), hitType));

        // First obtain the longitudinal position of spine hits
        LongitudinalPositionMap longitudinalPositionMap;
        this->ObtainLongitudinalDecomposition(spineTwoDSlidingFit, showerSpineHitList, longitudinalPositionMap);

        // Obtain spine energy profile
        EnergySpectrumMap energySpectrumMap;
        this->GetEnergyDistribution(showerSpineHitList, longitudinalPositionMap, energySpectrumMap);

        // Now find the shower start position/direction
        const bool isEndDownstream(peakDirection.GetZ() > 0.f);

        this->FindShowerStartAndDirection(pShowerPfo, hitType, spineTwoDSlidingFit, energySpectrumMap, showerSpineHitList, isEndDownstream,
            showerStartPosition, showerStartDirection);
    }
    catch (...)
    {
        return STATUS_CODE_FAILURE;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerStartFinderTool::ObtainLongitudinalDecomposition(const TwoDSlidingFitResult &spineTwoDSlidingFit,
    const CaloHitList &showerSpineHitList, LongitudinalPositionMap &longitudinalPositionMap) const
{
    // Find hits in each layer
    LayerToHitMap layerToHitMap;

    for (const CaloHit *const pCaloHit : showerSpineHitList)
    {
        float hitL(0.f), hitT(0.f);
        const CartesianVector hitPosition(pCaloHit->GetPositionVector());

        spineTwoDSlidingFit.GetLocalPosition(hitPosition, hitL, hitT);
        layerToHitMap[spineTwoDSlidingFit.GetLayer(hitL)].push_back(pCaloHit);
    }

    // Walk through layers and store the longitudinal distance of each hit
    // IMPORTANT: layer positions correspond to the middle of the layer
    // IMPORTANT: the longitudinal distance is measured relative to the first layer position (where l = 0)
    // IMPORTANT: at any point, the running distance is that at the lowest neighbour central position
    float runningDistance(0);
    const LayerFitResultMap &layerFitResultMap(spineTwoDSlidingFit.GetLayerFitResultMap());

    for (auto iter = layerToHitMap.begin(); iter != layerToHitMap.end(); ++iter)
    {
        const int layer(iter->first);
        const float layerL(layerFitResultMap.at(layer).GetL());

        // Get global positions of layer and its neighbours
        const int higherLayer(std::next(iter) == layerToHitMap.end() ? layer : std::next(iter)->first);
        const int middleLayer(layer);
        const int lowerLayer(iter == layerToHitMap.begin() ? layer : std::prev(iter)->first);
        CartesianVector lowerLayerPosition(0.f, 0.f, 0.f), middleLayerPosition(0.f, 0.f, 0.f), higherLayerPosition(0.f, 0.f, 0.f);

        spineTwoDSlidingFit.GetGlobalPosition(layerFitResultMap.at(lowerLayer).GetL(), layerFitResultMap.at(lowerLayer).GetFitT(), lowerLayerPosition);
        spineTwoDSlidingFit.GetGlobalPosition(layerFitResultMap.at(middleLayer).GetL(), layerFitResultMap.at(middleLayer).GetFitT(), middleLayerPosition);
        spineTwoDSlidingFit.GetGlobalPosition(layerFitResultMap.at(higherLayer).GetL(), layerFitResultMap.at(higherLayer).GetFitT(), higherLayerPosition);

        const float layerLength = std::next(iter) == layerToHitMap.end()
                                      ? 0.f
                                      : iter == layerToHitMap.begin() ? 0.f : (middleLayerPosition - lowerLayerPosition).GetMagnitude();

        for (const CaloHit *const pCaloHit : iter->second)
        {
            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());

            float hitL(0.f), hitT(0.f);
            spineTwoDSlidingFit.GetLocalPosition(hitPosition, hitL, hitT);

            // Is the hit below or above the middle of the layer?
            // Therefore, determine which layer positions the hit is between
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

            // Determine the axis to project onto
            // Calculate the distance of axis between its two nearest neighbours
            const CartesianVector displacement((higherLayerPosition - lowerLayerPosition).GetUnitVector());
            float longitudinalDisplacement = (displacement.GetDotProduct(hitPosition - lowerLayerPosition) + runningDistance);

            // If we've passed the central position, we need to add on the layer length
            longitudinalDisplacement += (hitL > layerL ? layerLength : 0.f);

            longitudinalPositionMap[pCaloHit] = longitudinalDisplacement;
        }

        runningDistance += layerLength;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerStartFinderTool::GetEnergyDistribution(const CaloHitList &showerSpineHitList,
    const LongitudinalPositionMap &longitudinalPositionMap, EnergySpectrumMap &energySpectrumMap) const
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

void ShowerStartFinderTool::FindShowerStartAndDirection(const ParticleFlowObject *const pShowerPfo, const HitType hitType,
    const TwoDSlidingFitResult &spineTwoDSlidingFit, const EnergySpectrumMap &energySpectrumMap, const CaloHitList &showerSpineHitList,
    const bool isEndDownstream, CartesianVector &showerStartPosition, CartesianVector &showerStartDirection) const
{
    // Search for shower start at significant energy deviations
    int longitudinalStartBin(0);

    if (isEndDownstream)
    {
        longitudinalStartBin = this->FindShowerStartLongitudinalCoordinate(pShowerPfo, hitType, spineTwoDSlidingFit, energySpectrumMap,
            showerSpineHitList, isEndDownstream, energySpectrumMap.begin(), std::prev(energySpectrumMap.end()));
    }
    else
    {
        longitudinalStartBin = this->FindShowerStartLongitudinalCoordinate(pShowerPfo, hitType, spineTwoDSlidingFit, energySpectrumMap,
            showerSpineHitList, isEndDownstream, energySpectrumMap.rbegin(), std::prev(energySpectrumMap.rend()));
    }

    const float longitudinalStartCoordinate(longitudinalStartBin * m_longitudinalCoordinateBinSize);

    this->ConvertLongitudinalProjectionToGlobal(spineTwoDSlidingFit, longitudinalStartCoordinate, showerStartPosition, showerStartDirection);

    showerStartDirection = isEndDownstream ? showerStartDirection : showerStartDirection * (-1.f);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
int ShowerStartFinderTool::FindShowerStartLongitudinalCoordinate(const ParticleFlowObject *const pShowerPfo, const HitType hitType,
    const TwoDSlidingFitResult &spineTwoDSlidingFit, const EnergySpectrumMap &energySpectrumMap, const CaloHitList &showerSpineHitList,
    const bool isEndDownstream, const T startIter, const T endIter) const
{
    // Characterise initial energy
    float meanEnergy(0.f), energySigma(0.f);

    if (energySpectrumMap.size() > m_nInitialEnergyBins)
        this->CharacteriseInitialEnergy(energySpectrumMap, isEndDownstream, meanEnergy, energySigma);

    // Now loop through energySpectrumMap
    T iter(startIter);

    std::advance(iter, (energySpectrumMap.size() > m_nInitialEnergyBins) ? m_nInitialEnergyBins : (energySpectrumMap.size() - 1));

    for (; iter != endIter; iter++)
    {
        const float longitudinalCoordinate(iter->first * m_longitudinalCoordinateBinSize);
        const float energyDeviation((iter->second - meanEnergy) / energySigma);

        // Use energy and local topology to assess whether we are at the shower start
        if ((energyDeviation > m_minSigmaDeviation) &&
            this->IsShowerTopology(pShowerPfo, hitType, spineTwoDSlidingFit, longitudinalCoordinate, showerSpineHitList, isEndDownstream))
        {
            break;
        }
    }

    return iter->first;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerStartFinderTool::CharacteriseInitialEnergy(
    const EnergySpectrumMap &energySpectrumMap, const bool isEndDownstream, float &meanEnergy, float &energySigma) const
{
    auto energySpectrumIter(isEndDownstream ? energySpectrumMap.begin() : std::prev(energySpectrumMap.end()));

    for (unsigned int i = 0; i < m_nInitialEnergyBins; ++i)
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

bool ShowerStartFinderTool::IsShowerTopology(const ParticleFlowObject *const pShowerPfo, const HitType hitType, const TwoDSlidingFitResult &spineTwoDSlidingFit,
    const float longitudinalDistance, const CaloHitList &showerSpineHitList, const bool isEndDownstream) const
{
    CartesianVector showerStartPosition(0.f, 0.f, 0.f), showerStartDirection(0.f, 0.f, 0.f);
    this->ConvertLongitudinalProjectionToGlobal(spineTwoDSlidingFit, longitudinalDistance, showerStartPosition, showerStartDirection);

    // Build shower
    CartesianPointVector showerRegionPositionVector;
    if (this->BuildShowerRegion(pShowerPfo, hitType, showerSpineHitList, showerStartPosition, showerStartDirection, isEndDownstream,
            showerRegionPositionVector) != STATUS_CODE_SUCCESS)
    {
        return false;
    }

    // Characterise the shower
    bool isBetween(false);
    CartesianVector positiveEdgeStart(0.f, 0.f, 0.f), positiveEdgeEnd(0.f, 0.f, 0.f), positiveEdgeDirection(0.f, 0.f, 0.f);
    CartesianVector negativeEdgeStart(0.f, 0.f, 0.f), negativeEdgeEnd(0.f, 0.f, 0.f), negativeEdgeDirection(0.f, 0.f, 0.f);

    if (this->CharacteriseShowerTopology(showerRegionPositionVector, showerStartPosition, hitType, isEndDownstream, showerStartDirection,
            positiveEdgeStart, positiveEdgeEnd, negativeEdgeStart, negativeEdgeEnd, isBetween) != STATUS_CODE_SUCCESS)
    {
        return false;
    }

    if (!isBetween)
        return false;

    positiveEdgeStart = showerStartPosition;
    negativeEdgeStart = showerStartPosition;
    positiveEdgeDirection = positiveEdgeEnd - positiveEdgeStart;
    negativeEdgeDirection = negativeEdgeEnd - negativeEdgeStart;

    const float showerOpeningAngle(positiveEdgeDirection.GetOpeningAngle(negativeEdgeDirection) * 180.f / 3.14);

    if (showerOpeningAngle < m_minShowerOpeningAngle)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerStartFinderTool::ConvertLongitudinalProjectionToGlobal(const TwoDSlidingFitResult &spineTwoDSlidingFit,
    const float longitudinalDistance, CartesianVector &globalPosition, CartesianVector &globalDirection) const
{
    const LayerFitResultMap &layerFitResultMap(spineTwoDSlidingFit.GetLayerFitResultMap());

    float runningDistance(0.f);
    int showerStartLayer(0);
    CartesianVector previousLayerPosition(0.f, 0.f, 0.f);
    spineTwoDSlidingFit.GetGlobalPosition(layerFitResultMap.begin()->second.GetL(), layerFitResultMap.begin()->second.GetFitT(), previousLayerPosition);

    for (auto iter = std::next(layerFitResultMap.begin()); iter != layerFitResultMap.end(); ++iter)
    {
        const int layer(iter->first);
        showerStartLayer = layer;

        CartesianVector layerPosition(0.f, 0.f, 0.f);
        spineTwoDSlidingFit.GetGlobalPosition(layerFitResultMap.at(layer).GetL(), layerFitResultMap.at(layer).GetFitT(), layerPosition);

        runningDistance += (layerPosition - previousLayerPosition).GetMagnitude();

        if (runningDistance > longitudinalDistance)
            break;

        previousLayerPosition = layerPosition;
    }

    const float lCoordinate(layerFitResultMap.at(showerStartLayer).GetL()), tCoordinate(layerFitResultMap.at(showerStartLayer).GetFitT());
    const float localGradient(layerFitResultMap.at(showerStartLayer).GetGradient());

    spineTwoDSlidingFit.GetGlobalPosition(lCoordinate, tCoordinate, globalPosition);
    spineTwoDSlidingFit.GetGlobalDirection(localGradient, globalDirection);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerStartFinderTool::BuildShowerRegion(const ParticleFlowObject *const pShowerPfo, const HitType hitType,
    const CaloHitList &showerSpineHitList, const CartesianVector &showerStartPosition, const CartesianVector &showerStartDirection,
    const bool isEndDownstream, CartesianPointVector &showerRegionPositionVector) const
{
    CaloHitList viewShowerHitList;
    LArPfoHelper::GetCaloHits(pShowerPfo, hitType, viewShowerHitList);

    // Find shower region hits
    CaloHitList showerRegionHitList;

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

        showerRegionPositionVector.push_back(pCaloHit->GetPositionVector());
        showerRegionHitList.push_back(pCaloHit);
    }

    // Searching for a continuous shower, so truncate at first gap
    CartesianPointVector coordinateListP, coordinateListN;
    try
    {
        // Perform shower sliding fit
        const TwoDSlidingShowerFitResult twoDShowerSlidingFit(
            &showerRegionPositionVector, m_showerSlidingFitWindow, LArGeometryHelper::GetWirePitch(this->GetPandora(), hitType));
        const LayerFitResultMap &layerFitResultMapS(twoDShowerSlidingFit.GetShowerFitResult().GetLayerFitResultMap());
        const LayerFitResultMap &layerFitResultMapP(twoDShowerSlidingFit.GetPositiveEdgeFitResult().GetLayerFitResultMap());
        const LayerFitResultMap &layerFitResultMapN(twoDShowerSlidingFit.GetNegativeEdgeFitResult().GetLayerFitResultMap());

        float showerStartL(0.f), showerStartT(0.f);
        twoDShowerSlidingFit.GetShowerFitResult().GetLocalPosition(showerStartPosition, showerStartL, showerStartT);

        const int startLayer(twoDShowerSlidingFit.GetShowerFitResult().GetLayer(showerStartL));

        // Check for gaps in the shower
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

            if (separationSquared > (m_maxEdgeGap * m_maxEdgeGap))
                break;

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

            if (separationSquared > (m_maxEdgeGap * m_maxEdgeGap))
                break;

            float thisT(0.f);
            twoDShowerSlidingFit.GetShowerFitResult().GetLocalPosition(coordinateN, nMaximumL, thisT);

            previousCoordinateN = coordinateN;
        }

        // Now truncate the halo hit position vector
        showerRegionPositionVector.clear();

        for (const CaloHit *const pCaloHit : showerRegionHitList)
        {
            float thisL(0.f), thisT(0.f);
            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());

            twoDShowerSlidingFit.GetShowerFitResult().GetLocalPosition(hitPosition, thisL, thisT);

            // because the end of showers can get really narrow and thus fail the opening angle check
            if (isEndDownstream && (thisL > (m_longitudinalShowerFraction * std::max(pMaximumL, nMaximumL))))
                continue;

            if (!isEndDownstream && (thisL < (m_longitudinalShowerFraction * std::min(pMaximumL, nMaximumL))))
                continue;

            showerRegionPositionVector.push_back(pCaloHit->GetPositionVector());
        }
    }
    catch (const StatusCodeException &)
    {
        return STATUS_CODE_FAILURE;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerStartFinderTool::CharacteriseShowerTopology(const CartesianPointVector &showerRegionPositionVector, const CartesianVector &showerStartPosition,
    const HitType hitType, const bool isEndDownstream, const CartesianVector &showerStartDirection, CartesianVector &positiveEdgeStart,
    CartesianVector &positiveEdgeEnd, CartesianVector &negativeEdgeStart, CartesianVector &negativeEdgeEnd, bool &isBetween) const
{
    try
    {
        // Perform shower sliding fit
        const TwoDSlidingShowerFitResult twoDShowerSlidingFit(
            &showerRegionPositionVector, m_showerSlidingFitWindow, LArGeometryHelper::GetWirePitch(this->GetPandora(), hitType));
        const LayerFitResultMap &layerFitResultMapS(twoDShowerSlidingFit.GetShowerFitResult().GetLayerFitResultMap());
        const LayerFitResultMap &layerFitResultMapP(twoDShowerSlidingFit.GetPositiveEdgeFitResult().GetLayerFitResultMap());
        const LayerFitResultMap &layerFitResultMapN(twoDShowerSlidingFit.GetNegativeEdgeFitResult().GetLayerFitResultMap());

        float showerStartL(0.f), showerStartT(0.f);
        twoDShowerSlidingFit.GetShowerFitResult().GetLocalPosition(showerStartPosition, showerStartL, showerStartT);

        const int startLayer(twoDShowerSlidingFit.GetShowerFitResult().GetLayer(showerStartL));
        const int endLayer(twoDShowerSlidingFit.GetShowerFitResult().GetMaxLayer());

        // Does the shower look to originate from the shower start position
        int layerCount(0);
        bool isFirstBetween(false), isLastBetween(false);
        CartesianPointVector coordinateListP, coordinateListN;

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
            }

            if (layerFitResultMapN.end() != iterN)
            {
                twoDShowerSlidingFit.GetShowerFitResult().GetGlobalPosition(iterN->second.GetL(), iterN->second.GetFitT(), negativeEdgePosition);
                coordinateListN.push_back(negativeEdgePosition);
            }

            if ((layerFitResultMapP.end() != iterP) && (layerFitResultMapN.end() != iterN))
            {
                const CartesianVector positiveDisplacement((positiveEdgePosition - showerStartPosition).GetUnitVector());
                const bool isPositiveClockwise(this->IsClockwiseRotation(showerStartDirection, positiveDisplacement));

                const CartesianVector negativeDisplacement((negativeEdgePosition - showerStartPosition).GetUnitVector());
                const bool isNegativeClockwise(this->IsClockwiseRotation(showerStartDirection, negativeDisplacement));

                ++layerCount;

                if (layerCount == 1)
                    isFirstBetween = (isPositiveClockwise != isNegativeClockwise);

                isLastBetween = (isPositiveClockwise != isNegativeClockwise);
            }
        }

        isBetween = (isFirstBetween || isLastBetween);

        // Find the extremal points of the shower
        if (this->GetBoundaryExtremalPoints(twoDShowerSlidingFit, layerFitResultMapP, showerStartPosition, startLayer, endLayer,
                positiveEdgeStart, positiveEdgeEnd) != STATUS_CODE_SUCCESS)
        {
            return STATUS_CODE_FAILURE;
        }

        if (this->GetBoundaryExtremalPoints(twoDShowerSlidingFit, layerFitResultMapN, showerStartPosition, startLayer, endLayer,
                negativeEdgeStart, negativeEdgeEnd) != STATUS_CODE_SUCCESS)
        {
            return STATUS_CODE_FAILURE;
        }
    }
    catch (const StatusCodeException &)
    {
        return STATUS_CODE_FAILURE;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ShowerStartFinderTool::IsClockwiseRotation(const CartesianVector &showerStartDirection, const CartesianVector &displacementVector) const
{
    // Clockwise with respect to +ve Z
    const float openingAngleFromCore(showerStartDirection.GetOpeningAngle(displacementVector));

    const CartesianVector clockwiseRotation(
        displacementVector.GetZ() * std::sin(openingAngleFromCore) + displacementVector.GetX() * std::cos(openingAngleFromCore), 0.f,
        displacementVector.GetZ() * std::cos(openingAngleFromCore) - displacementVector.GetX() * std::sin(openingAngleFromCore));

    const CartesianVector anticlockwiseRotation(
        displacementVector.GetZ() * std::sin(-1.f * openingAngleFromCore) + displacementVector.GetX() * std::cos(-1.f * openingAngleFromCore), 0.f,
        displacementVector.GetZ() * std::cos(-1.f * openingAngleFromCore) - displacementVector.GetX() * std::sin(-1.f * openingAngleFromCore));

    const float clockwiseT((clockwiseRotation - showerStartDirection).GetMagnitude());
    const float anticlockwiseT((anticlockwiseRotation - showerStartDirection).GetMagnitude());
    const bool isClockwise(clockwiseT < anticlockwiseT);

    return isClockwise;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerStartFinderTool::GetBoundaryExtremalPoints(const TwoDSlidingShowerFitResult &showerTwoDSlidingFit,
    const LayerFitResultMap &layerFitResultMap, const CartesianVector &showerStartPosition, const int showerStartLayer,
    const int showerEndLayer, CartesianVector &boundaryEdgeStart, CartesianVector &boundaryEdgeEnd) const
{
    int boundaryStartLayer(std::numeric_limits<int>::max());
    int boundaryEndLayer(std::numeric_limits<int>::max());

    for (auto &entry : layerFitResultMap)
    {
        const int bestStartSeparation(std::abs(showerStartLayer - boundaryStartLayer));
        const int thisStartSeparation(std::abs(showerStartLayer - entry.first));

        if (((bestStartSeparation - thisStartSeparation) > 0) && (thisStartSeparation < bestStartSeparation))
            boundaryStartLayer = entry.first;

        const int bestEndSeparation(std::abs(showerEndLayer - boundaryEndLayer));
        const int thisEndSeparation(std::abs(showerEndLayer - entry.first));

        if (((bestEndSeparation - thisEndSeparation) > 0) && (thisEndSeparation < bestEndSeparation))
            boundaryEndLayer = entry.first;
    }

    if (std::abs(showerStartLayer - boundaryStartLayer) > m_maxLayerSeparation)
        return STATUS_CODE_FAILURE;

    const float showerStartBoundaryLocalL(layerFitResultMap.at(boundaryStartLayer).GetL()),
        showerStartBoundaryLocalT(layerFitResultMap.at(boundaryStartLayer).GetFitT());
    const float showerEndBoundaryLocalL(layerFitResultMap.at(boundaryEndLayer).GetL()),
        showerEndBoundaryLocalT(layerFitResultMap.at(boundaryEndLayer).GetFitT());

    showerTwoDSlidingFit.GetShowerFitResult().GetGlobalPosition(showerStartBoundaryLocalL, showerStartBoundaryLocalT, boundaryEdgeStart);
    showerTwoDSlidingFit.GetShowerFitResult().GetGlobalPosition(showerEndBoundaryLocalL, showerEndBoundaryLocalT, boundaryEdgeEnd);

    // Swap start and end if needs be
    if ((showerStartPosition - boundaryEdgeEnd).GetMagnitudeSquared() < (showerStartPosition - boundaryEdgeStart).GetMagnitude())
    {
        const CartesianVector startTemp(boundaryEdgeStart), endTemp(boundaryEdgeEnd);

        boundaryEdgeStart = endTemp;
        boundaryEdgeEnd = startTemp;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerStartFinderTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SpineSlidingFitWindow", m_spineSlidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "LongitudinalCoordinateBinSize", m_longitudinalCoordinateBinSize));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "nInitialEnergyBins", m_nInitialEnergyBins));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinSigmaDeviation", m_minSigmaDeviation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxEdgeGap", m_maxEdgeGap));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "LongitudinalShowerFraction", m_longitudinalShowerFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinShowerOpeningAngle", m_minShowerOpeningAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MolliereRadius", m_molliereRadius));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ShowerSlidingFitWindow", m_showerSlidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxLayerSeparation", m_maxLayerSeparation));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
