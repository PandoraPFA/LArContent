/**
 *  @file   larpandoracontent/LArShowerRefinement/PeakDirectionFinderTool.cc
 *
 *  @brief  Implementation of the peak direction finder tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArShowerRefinement/PeakDirectionFinderTool.h"

using namespace pandora;

namespace lar_content
{

PeakDirectionFinderTool::PeakDirectionFinderTool() :
    m_pathwaySearchRegion(10.f),
    m_theta0XZBinSize(0.005f),
    m_smoothingWindow(1)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PeakDirectionFinderTool::Run(const ParticleFlowObject *const pShowerPfo, const CartesianVector &nuVertex3D, const CaloHitList *const pViewHitList, 
    const HitType hitType, CartesianPointVector &peakDirectionVector)
{
    CaloHitList viewShowerHitList;
    LArPfoHelper::GetCaloHits(pShowerPfo, hitType, viewShowerHitList);

    const CartesianVector nuVertex2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertex3D, hitType));

    // Get event hits in region of interest
    CaloHitList viewROIHits;
    this->CollectHitsWithinROI(viewShowerHitList, pViewHitList, nuVertex2D, viewROIHits);

    if (viewROIHits.empty())
        return STATUS_CODE_NOT_FOUND;

    // Fill angular decomposition map
    AngularDecompositionMap angularDecompositionMap;
    this->FillAngularDecompositionMap(viewROIHits, nuVertex2D, angularDecompositionMap);

    if (angularDecompositionMap.empty())
        return STATUS_CODE_NOT_FOUND;

    // Smooth angular decomposition map
    this->SmoothAngularDecompositionMap(angularDecompositionMap);
    
    // Get peak directions
    this->RetrievePeakDirections(angularDecompositionMap, peakDirectionVector);

    if (peakDirectionVector.empty())
        return STATUS_CODE_NOT_FOUND;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PeakDirectionFinderTool::CollectHitsWithinROI(const CaloHitList &showerHitList, const CaloHitList *const pViewHitList, 
    const CartesianVector &nuVertex2D, CaloHitList &viewROIHits) const
{
    float lowestTheta(std::numeric_limits<float>::max()), highestTheta((-1.f) * std::numeric_limits<float>::max());

    this->GetAngularExtrema(showerHitList, nuVertex2D, lowestTheta, highestTheta);
    this->CollectHitsWithinExtrema(pViewHitList, nuVertex2D, lowestTheta, highestTheta, viewROIHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PeakDirectionFinderTool::GetAngularExtrema(const CaloHitList &showerHitList, const CartesianVector &nuVertex2D, 
    float &lowestTheta, float &highestTheta) const
{
    lowestTheta = std::numeric_limits<float>::max();
    highestTheta = (-1.f) * std::numeric_limits<float>::max();

    const CartesianVector xAxis(1.f, 0.f, 0.f);

    for (const CaloHit *const pCaloHit : showerHitList)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        const CartesianVector displacementVector(hitPosition - nuVertex2D);

        float theta0XZ(xAxis.GetOpeningAngle(displacementVector));

        if (displacementVector.GetZ() < 0.f)
            theta0XZ = (2.0 * M_PI) - theta0XZ;

        if (theta0XZ < lowestTheta)
            lowestTheta = theta0XZ;

        if (theta0XZ > highestTheta)
            highestTheta = theta0XZ;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PeakDirectionFinderTool::CollectHitsWithinExtrema(const CaloHitList *const pViewHitList, const CartesianVector &nuVertex2D, 
    const float lowestTheta, const float highestTheta, CaloHitList &viewROIHits) const
{
    const CartesianVector xAxis(1.f, 0.f, 0.f);

    for (const CaloHit *const pCaloHit : *pViewHitList)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        const CartesianVector displacementVector(hitPosition - nuVertex2D);

        float theta0XZ(xAxis.GetOpeningAngle(displacementVector));

        if (displacementVector.GetZ() < 0.f)
            theta0XZ = (2.0 * M_PI) - theta0XZ;

        if ((theta0XZ < lowestTheta) || (theta0XZ > highestTheta))
            continue;

        viewROIHits.push_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PeakDirectionFinderTool::FillAngularDecompositionMap(const CaloHitList &viewShowerHitList, const CartesianVector &nuVertex2D, 
    AngularDecompositionMap &angularDecompositionMap) const
{
    const CartesianVector xAxis(1.f, 0.f, 0.f);

    for (const CaloHit *const pCaloHit : viewShowerHitList)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        const CartesianVector displacementVector(hitPosition - nuVertex2D);

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

void PeakDirectionFinderTool::SmoothAngularDecompositionMap(AngularDecompositionMap &angularDecompositionMap) const
{
    const int loopMin = (-1) * m_smoothingWindow;
    const int loopMax = m_smoothingWindow;

    const AngularDecompositionMap angularDecompositionMapTemp(angularDecompositionMap);

    angularDecompositionMap.clear();

    for (const auto &entry : angularDecompositionMapTemp)
    {
        float total(0.f);
        const int currentBin(entry.first);

        for (int binOffset = loopMin; binOffset <= loopMax; ++binOffset)
        {
            const int contributingBin = currentBin + binOffset;
            total += (angularDecompositionMapTemp.find(contributingBin) == angularDecompositionMapTemp.end()) ? 0.f : angularDecompositionMapTemp.at(contributingBin);
        }

        angularDecompositionMap[currentBin] = total / static_cast<float>((2.0 * m_smoothingWindow) + 1);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PeakDirectionFinderTool::RetrievePeakDirections(const AngularDecompositionMap &angularDecompositionMap, CartesianPointVector &peakDirectionVector) const
{
    IntVector orderedBinIndexVector;

    for (const auto &entry : angularDecompositionMap)
        orderedBinIndexVector.push_back(entry.first);

    // Order peak bin vector from highest to lowest bin height
    // Tie-break: highest index wins
    std::sort(orderedBinIndexVector.begin(), orderedBinIndexVector.end(),
        [&angularDecompositionMap](const int a, const int b) -> bool {

                  float aWeight = angularDecompositionMap.at(a);
                  float bWeight = angularDecompositionMap.at(b);

                  if (std::fabs(aWeight - bWeight) < std::numeric_limits<float>::epsilon())
                      return a > b;
                  else
                      return aWeight > bWeight;
              });


    for (int binIndex : orderedBinIndexVector)
    {
        const float binWeight(angularDecompositionMap.at(binIndex));

        int precedingBin(binIndex - 1);
        bool foundPreceeding(false);

        while (!foundPreceeding)
        {
            if ((angularDecompositionMap.find(precedingBin) == angularDecompositionMap.end()) ||
                (std::fabs(angularDecompositionMap.at(precedingBin) - angularDecompositionMap.at(binIndex)) > std::numeric_limits<float>::epsilon()))
            {
                foundPreceeding = true;
                break;
            }

            --precedingBin;
        }

        int followingBin(binIndex + 1);
        bool foundFollowing(false);

        while (!foundFollowing)
        {
            if ((angularDecompositionMap.find(followingBin) == angularDecompositionMap.end()) ||
                (std::fabs(angularDecompositionMap.at(followingBin) - angularDecompositionMap.at(binIndex)) > std::numeric_limits<float>::epsilon()))
            {
                foundFollowing = true;
                break;
            }

            ++followingBin;
        }

        const float precedingBinWeight(angularDecompositionMap.find(precedingBin) == angularDecompositionMap.end() ? 0.f : angularDecompositionMap.at(precedingBin));
        const float followingBinWeight(angularDecompositionMap.find(followingBin) == angularDecompositionMap.end() ? 0.f : angularDecompositionMap.at(followingBin));

        if ((binWeight < precedingBinWeight) || (binWeight < followingBinWeight))
            continue;

        // Convert to a direction
        const float theta0XZ(binIndex * m_theta0XZBinSize);
        const CartesianVector peakDirection(std::cos(theta0XZ), 0.f, std::sin(theta0XZ));

        peakDirectionVector.push_back(peakDirection);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------


StatusCode PeakDirectionFinderTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "PathwaySearchRegion", m_pathwaySearchRegion));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "Theta0XZBinSize", m_theta0XZBinSize));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "SmoothingWindow", m_smoothingWindow));

    return STATUS_CODE_SUCCESS;
}


} // namespace lar_content
