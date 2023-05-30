/**
 *  @file   larpandoracontent/LArShowerRefinement/ShowerSpineFinderTool.cc
 *
 *  @brief  Implementation of the shower spine finder tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArHelpers/LArHitWidthHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArConnectionPathwayHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArShowerRefinement/ShowerSpineFinderTool.h"

using namespace pandora;

namespace lar_content
{

ShowerSpineFinderTool::ShowerSpineFinderTool() :
    m_hitThresholdForSpine(7),
    m_growingFitInitialLength(10.f),
    m_initialFitDistanceToLine(0.5f),
    m_minInitialHitsFound(7),
    m_maxFittingHits(15),
    m_localSlidingFitWindow(10),
    m_growingFitSegmentLength(5.f),
    m_highResolutionSlidingFitWindow(5),
    m_distanceToLine(0.75f),
    m_hitConnectionDistance(0.75f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerSpineFinderTool::Run(const CartesianVector &nuVertex3D, const CaloHitList *const pViewHitList, const HitType hitType, 
    const CartesianVector &peakDirection, CaloHitList &unavailableHitList, CaloHitList &showerSpineHitList)
{
    const CartesianVector nuVertex2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertex3D, hitType));

    this->FindShowerSpine(pViewHitList, nuVertex2D, peakDirection, unavailableHitList, showerSpineHitList);

     // Demand that spine is significant, be lenient here as some have small stubs and a gap
    if (showerSpineHitList.size() < m_hitThresholdForSpine) 
        return STATUS_CODE_NOT_FOUND;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerSpineFinderTool::FindShowerSpine(const CaloHitList *const pViewHitList, const CartesianVector &nuVertex2D, 
    const CartesianVector &initialDirection, CaloHitList &unavailableHitList, CaloHitList &showerSpineHitList) const
{
    // Use initial direction to find seed hits for a starting fit
    float highestL(0.f);
    CartesianPointVector runningFitPositionVector;

    for (const CaloHit *const pCaloHit : *pViewHitList)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        const CartesianVector &displacementVector(hitPosition - nuVertex2D);
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

            runningFitPositionVector.push_back(hitPosition);
        }
    }

    // Require significant number of initial hits
    if (runningFitPositionVector.size() < m_minInitialHitsFound)
    {
        showerSpineHitList.clear();
        return;
    }

    // Perform a running fit to collect a pathway of hits
    unsigned int count(0);
    bool hitsCollected(true);
    const bool isEndDownstream(initialDirection.GetZ() > 0.f);
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    CartesianVector extrapolatedStartPosition(nuVertex2D), extrapolatedDirection(initialDirection);
    CartesianVector extrapolatedEndPosition(extrapolatedStartPosition + (extrapolatedDirection * highestL));

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

            const TwoDSlidingFitResult extrapolatedFit(&runningFitPositionVector, m_localSlidingFitWindow, slidingFitPitch);            

            extrapolatedStartPosition = count == 1 ? extrapolatedEndPosition : isEndDownstream ? extrapolatedFit.GetGlobalMaxLayerPosition() : extrapolatedFit.GetGlobalMinLayerPosition();
            extrapolatedDirection = isEndDownstream ? extrapolatedFit.GetGlobalMaxLayerDirection() : extrapolatedFit.GetGlobalMinLayerDirection() * (-1.f);
            extrapolatedEndPosition = extrapolatedStartPosition + (extrapolatedDirection * m_growingFitSegmentLength);


            hitsCollected = this->CollectSubsectionHits(extrapolatedFit, extrapolatedStartPosition, extrapolatedEndPosition, extrapolatedDirection,
                isEndDownstream, pViewHitList, runningFitPositionVector, unavailableHitList, showerSpineHitList);

            // If no hits found, as a final effort, reduce the sliding fit window
            if (!hitsCollected)
            {
                const TwoDSlidingFitResult microExtrapolatedFit(&runningFitPositionVector, m_highResolutionSlidingFitWindow, slidingFitPitch);
            
                extrapolatedStartPosition = count == 1 ? extrapolatedStartPosition : isEndDownstream ? microExtrapolatedFit.GetGlobalMaxLayerPosition() : 
                    microExtrapolatedFit.GetGlobalMinLayerPosition();
                extrapolatedDirection = isEndDownstream ? microExtrapolatedFit.GetGlobalMaxLayerDirection() : microExtrapolatedFit.GetGlobalMinLayerDirection() * (-1.f);
                extrapolatedEndPosition = extrapolatedStartPosition + (extrapolatedDirection * m_growingFitSegmentLength);

                hitsCollected = this->CollectSubsectionHits(microExtrapolatedFit, extrapolatedStartPosition, extrapolatedEndPosition, extrapolatedDirection,
                    isEndDownstream, pViewHitList, runningFitPositionVector, unavailableHitList, showerSpineHitList);
            }
        }
        catch (const StatusCodeException &)
        {
            return;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ShowerSpineFinderTool::CollectSubsectionHits(const TwoDSlidingFitResult &extrapolatedFit, const CartesianVector &extrapolatedStartPosition, 
    const CartesianVector &extrapolatedEndPosition, const CartesianVector &extrapolatedDirection, const bool isEndDownstream, 
    const CaloHitList *const pViewHitList, CartesianPointVector &runningFitPositionVector, CaloHitList &unavailableHitList, CaloHitList &showerSpineHitList) const
{ 
    float extrapolatedStartL(0.f), extrapolatedStartT(0.f);
    extrapolatedFit.GetLocalPosition(extrapolatedStartPosition, extrapolatedStartL, extrapolatedStartT);

    float extrapolatedEndL(0.f), extrapolatedEndT(0.f);
    extrapolatedFit.GetLocalPosition(extrapolatedEndPosition, extrapolatedEndL, extrapolatedEndT);

    CaloHitList collectedHits;

    for (const CaloHit *const pCaloHit : *pViewHitList)
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

    // Now find a continuous path of collected hits
    this->CollectConnectedHits(collectedHits, extrapolatedStartPosition, extrapolatedDirection, runningFitPositionVector, showerSpineHitList);

    const int nFinalHits(showerSpineHitList.size());

    return (nFinalHits != nInitialHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ShowerSpineFinderTool::IsCloseToLine(const CartesianVector &hitPosition, const CartesianVector &lineStart, 
    const CartesianVector &lineDirection, const float distanceToLine) const
{
    const float transverseDistanceFromLine(lineDirection.GetCrossProduct(hitPosition - lineStart).GetMagnitude());
    
    if (transverseDistanceFromLine > distanceToLine)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerSpineFinderTool::CollectConnectedHits(const CaloHitList &collectedHits, const CartesianVector &extrapolatedStartPosition, 
    const CartesianVector &extrapolatedDirection, CartesianPointVector &runningFitPositionVector, CaloHitList &showerSpineHitList) const
{
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
                if (LArHitWidthHelper::GetClosestDistance(pCaloHit, showerSpineHitList) > m_hitConnectionDistance)
                    continue;

                hitPosition = LArHitWidthHelper::GetClosestPointToLine2D(extrapolatedStartPosition, extrapolatedDirection, pCaloHit);
            }

            found = true;

            runningFitPositionVector.push_back(hitPosition);
            showerSpineHitList.push_back(pCaloHit);
         }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ShowerSpineFinderTool::GetClosestDistance(const CartesianVector &position, const CartesianPointVector &testPositions) const
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

StatusCode ShowerSpineFinderTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "HitThresholdForSpine", m_hitThresholdForSpine));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "GrowingFitInitialLength", m_growingFitInitialLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "InitialFitDistanceToLine", m_initialFitDistanceToLine));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinInitialHitsFound", m_minInitialHitsFound));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxFittingHits", m_maxFittingHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "LocalSlidingFitWindow", m_localSlidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "GrowingFitSegmentLength", m_growingFitSegmentLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "HighResolutionSlidingFitWindow", m_highResolutionSlidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "DistanceToLine", m_distanceToLine));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "HitConnectionDistance", m_hitConnectionDistance));

    return STATUS_CODE_SUCCESS;
}


} // namespace lar_content
