/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/LongitudinalTrackHitsBaseTool.cc
 *
 *  @brief  Implementation of the longitudinal track hits base tool.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArThreeDReco/LArHitCreation/LongitudinalTrackHitsBaseTool.h"
#include "larpandoracontent/LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

LongitudinalTrackHitsBaseTool::LongitudinalTrackHitsBaseTool() :
    m_vtxDisplacementCutSquared(5.f * 5.f),
    m_minTrackLengthSquared(7.5f * 7.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LongitudinalTrackHitsBaseTool::GetTrackHits3D(const CaloHitVector &inputTwoDHits, const MatchedSlidingFitMap &inputSlidingFitMap,
    ProtoHitVector &protoHitVector) const
{
    MatchedSlidingFitMap matchedSlidingFitMap;
    CartesianVector vtx3D(0.f, 0.f, 0.f), end3D(0.f, 0.f, 0.f);
    this->GetVertexAndEndPositions(inputSlidingFitMap, matchedSlidingFitMap, vtx3D, end3D);

    for (const CaloHit *const pCaloHit2D : inputTwoDHits)
    {
        try
        {
            ProtoHit protoHit(pCaloHit2D);
            this->GetLongitudinalTrackHit3D(matchedSlidingFitMap, vtx3D, end3D, protoHit);

            if (protoHit.IsPositionSet() && (protoHit.GetChi2() < m_chiSquaredCut))
                protoHitVector.push_back(protoHit);
        }
        catch (StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LongitudinalTrackHitsBaseTool::GetVertexAndEndPositions(const MatchedSlidingFitMap &inputSlidingFitMap,
    MatchedSlidingFitMap &outputSlidingFitMap, CartesianVector &bestVtx3D, CartesianVector &bestEnd3D) const
{
    // TODO Tidy up: The code below is quite repetitive...
    MatchedSlidingFitMap::const_iterator iterU = inputSlidingFitMap.find(TPC_VIEW_U);
    const bool foundU(inputSlidingFitMap.end() != iterU);

    MatchedSlidingFitMap::const_iterator iterV = inputSlidingFitMap.find(TPC_VIEW_V);
    const bool foundV(inputSlidingFitMap.end() != iterV);

    MatchedSlidingFitMap::const_iterator iterW = inputSlidingFitMap.find(TPC_VIEW_W);
    const bool foundW(inputSlidingFitMap.end() != iterW);

    bool useU(false), useV(false), useW(false);
    float bestChi2(std::numeric_limits<float>::max());

    for (unsigned int iPermutation = 0; iPermutation < 4; ++iPermutation)
    {
        const bool isForwardU((1 == iPermutation) ? false : true);
        const bool isForwardV((2 == iPermutation) ? false : true);
        const bool isForwardW((3 == iPermutation) ? false : true);

        CartesianVector vtxU(0.f, 0.f, 0.f), endU(0.f, 0.f, 0.f);
        CartesianVector vtxV(0.f, 0.f, 0.f), endV(0.f, 0.f, 0.f);
        CartesianVector vtxW(0.f, 0.f, 0.f), endW(0.f, 0.f, 0.f);

        if (foundU)
        {
            const TwoDSlidingFitResult &slidingFitResultU = iterU->second;
            vtxU = (isForwardU ? slidingFitResultU.GetGlobalMinLayerPosition() : slidingFitResultU.GetGlobalMaxLayerPosition());
            endU = (isForwardU ? slidingFitResultU.GetGlobalMaxLayerPosition() : slidingFitResultU.GetGlobalMinLayerPosition());
        }

        if (foundV)
        {
            const TwoDSlidingFitResult &slidingFitResultV = iterV->second;
            vtxV = (isForwardV ? slidingFitResultV.GetGlobalMinLayerPosition() : slidingFitResultV.GetGlobalMaxLayerPosition());
            endV = (isForwardV ? slidingFitResultV.GetGlobalMaxLayerPosition() : slidingFitResultV.GetGlobalMinLayerPosition());
        }

        if (foundW)
        {
            const TwoDSlidingFitResult &slidingFitResultW = iterW->second;
            vtxW = (isForwardW ? slidingFitResultW.GetGlobalMinLayerPosition() : slidingFitResultW.GetGlobalMaxLayerPosition());
            endW = (isForwardW ? slidingFitResultW.GetGlobalMaxLayerPosition() : slidingFitResultW.GetGlobalMinLayerPosition());
        }

        CartesianVector vtx3D(0.f, 0.f, 0.f), end3D(0.f, 0.f, 0.f);
        float vtxChi2(std::numeric_limits<float>::max()), endChi2(std::numeric_limits<float>::max());

        if (foundU && foundV)
        {
            this->UpdateBestPosition(TPC_VIEW_U, TPC_VIEW_V, vtxU, vtxV, vtx3D, vtxChi2);
            this->UpdateBestPosition(TPC_VIEW_U, TPC_VIEW_V, endU, endV, end3D, endChi2);
        }

        if (foundV && foundW)
        {
            this->UpdateBestPosition(TPC_VIEW_V, TPC_VIEW_W, vtxV, vtxW, vtx3D, vtxChi2);
            this->UpdateBestPosition(TPC_VIEW_V, TPC_VIEW_W, endV, endW, end3D, endChi2);
        }

        if (foundW && foundU)
        {
            this->UpdateBestPosition(TPC_VIEW_W, TPC_VIEW_U, vtxW, vtxU, vtx3D, vtxChi2);
            this->UpdateBestPosition(TPC_VIEW_W, TPC_VIEW_U, endW, endU, end3D, endChi2);
        }

        bool matchedU(false), matchedV(false), matchedW(false);
        unsigned int matchedViews(0);

        if (foundU)
        {
            const CartesianVector projVtxU(LArGeometryHelper::ProjectPosition(this->GetPandora(), vtx3D, TPC_VIEW_U));
            const CartesianVector projEndU(LArGeometryHelper::ProjectPosition(this->GetPandora(), end3D, TPC_VIEW_U));

            if((endU - vtxU).GetMagnitudeSquared() > m_minTrackLengthSquared &&
               (projVtxU - vtxU).GetMagnitudeSquared() < std::min(m_vtxDisplacementCutSquared, (projVtxU - endU).GetMagnitudeSquared()) &&
               (projEndU - endU).GetMagnitudeSquared() < std::min(m_vtxDisplacementCutSquared, (projEndU - vtxU).GetMagnitudeSquared()))
            {
                matchedU = true;  ++matchedViews;
            }
        }

        if (foundV)
        {
            const CartesianVector projVtxV(LArGeometryHelper::ProjectPosition(this->GetPandora(), vtx3D, TPC_VIEW_V));
            const CartesianVector projEndV(LArGeometryHelper::ProjectPosition(this->GetPandora(), end3D, TPC_VIEW_V));

            if((endV - vtxV).GetMagnitudeSquared() > m_minTrackLengthSquared &&
               (projVtxV - vtxV).GetMagnitudeSquared() < std::min(m_vtxDisplacementCutSquared, (projVtxV - endV).GetMagnitudeSquared()) &&
               (projEndV - endV).GetMagnitudeSquared() < std::min(m_vtxDisplacementCutSquared, (projEndV - vtxV).GetMagnitudeSquared()))
            {
                matchedV = true;  ++matchedViews;
            }
        }

        if (foundW)
        {
            const CartesianVector projVtxW(LArGeometryHelper::ProjectPosition(this->GetPandora(), vtx3D, TPC_VIEW_W));
            const CartesianVector projEndW(LArGeometryHelper::ProjectPosition(this->GetPandora(), end3D, TPC_VIEW_W));

            if((endW - vtxW).GetMagnitudeSquared() > m_minTrackLengthSquared &&
               (projVtxW - vtxW).GetMagnitudeSquared() < std::min(m_vtxDisplacementCutSquared, (projVtxW - endW).GetMagnitudeSquared()) &&
               (projEndW - endW).GetMagnitudeSquared() < std::min(m_vtxDisplacementCutSquared, (projEndW - vtxW).GetMagnitudeSquared()))
            {
                matchedW = true;  ++matchedViews;
            }
        }

        if (matchedViews < 2)
            continue;

        if (vtxChi2 + endChi2 < bestChi2)
        {
            useU = matchedU;
            useV = matchedV;
            useW = matchedW;

            bestVtx3D = vtx3D;
            bestEnd3D = end3D;
            bestChi2 = vtxChi2 + endChi2;
        }
    }

    if (useU)
        outputSlidingFitMap.insert(MatchedSlidingFitMap::value_type(iterU->first, iterU->second));

    if (useV)
        outputSlidingFitMap.insert(MatchedSlidingFitMap::value_type(iterV->first, iterV->second));

    if (useW)
        outputSlidingFitMap.insert(MatchedSlidingFitMap::value_type(iterW->first, iterW->second));

    if (outputSlidingFitMap.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LongitudinalTrackHitsBaseTool::UpdateBestPosition(const HitType hitType1, const HitType hitType2, const CartesianVector &vtx1,
    const CartesianVector &vtx2, CartesianVector &bestVtx, float &bestChi2) const
{
    CartesianVector mergedVtx(0.f, 0.f, 0.f);
    float mergedChi2(std::numeric_limits<float>::max());

    LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), hitType1, hitType2, vtx1, vtx2, mergedVtx, mergedChi2);

    if (mergedChi2 < bestChi2)
    {
        bestVtx = mergedVtx;
        bestChi2 = mergedChi2;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LongitudinalTrackHitsBaseTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    float vtxDisplacementCut = std::sqrt(m_vtxDisplacementCutSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VertexDisplacementCut", vtxDisplacementCut));
    m_vtxDisplacementCutSquared = vtxDisplacementCut * vtxDisplacementCut;

    float minTrackLength = std::sqrt(m_minTrackLengthSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinTrackLength", minTrackLength));
    m_minTrackLengthSquared =  minTrackLength * minTrackLength;

    return TrackHitsBaseTool::ReadSettings(xmlHandle);
}

} // namespace lar_content
