/**
 *  @file   larpandoracontent/LArThreeDReco/LArLongitudinalTrackMatching/ThreeViewLongitudinalTracksAlgorithm.cc
 *
 *  @brief  Implementation of the three view longitudinal tracks algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArThreeDReco/LArLongitudinalTrackMatching/ThreeViewLongitudinalTracksAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ThreeViewLongitudinalTracksAlgorithm::ThreeViewLongitudinalTracksAlgorithm() :
    m_nMaxTensorToolRepeats(1000),
    m_vertexChi2Cut(10.f),
    m_reducedChi2Cut(5.f),
    m_samplingPitch(1.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewLongitudinalTracksAlgorithm::CalculateOverlapResult(
    const Cluster *const pClusterU, const Cluster *const pClusterV, const Cluster *const pClusterW)
{
    LongitudinalOverlapResult overlapResult;
    this->CalculateOverlapResult(pClusterU, pClusterV, pClusterW, overlapResult);

    if (overlapResult.IsInitialized())
        this->GetMatchingControl().GetOverlapTensor().SetOverlapResult(pClusterU, pClusterV, pClusterW, overlapResult);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewLongitudinalTracksAlgorithm::CalculateOverlapResult(const Cluster *const pClusterU, const Cluster *const pClusterV,
    const Cluster *const pClusterW, LongitudinalOverlapResult &longitudinalOverlapResult)
{
    const TwoDSlidingFitResult &slidingFitResultU(this->GetCachedSlidingFitResult(pClusterU));
    const TwoDSlidingFitResult &slidingFitResultV(this->GetCachedSlidingFitResult(pClusterV));
    const TwoDSlidingFitResult &slidingFitResultW(this->GetCachedSlidingFitResult(pClusterW));

    // Loop over possible permutations of cluster direction
    TrackOverlapResult bestOverlapResult;

    for (unsigned int iPermutation = 0; iPermutation < 4; ++iPermutation)
    {
        const bool isForwardU((1 == iPermutation) ? false : true);
        const bool isForwardV((2 == iPermutation) ? false : true);
        const bool isForwardW((3 == iPermutation) ? false : true);

        // Get 2D start and end positions from each sliding fit for this permutation
        const CartesianVector vtxU((isForwardU) ? slidingFitResultU.GetGlobalMinLayerPosition() : slidingFitResultU.GetGlobalMaxLayerPosition());
        const CartesianVector endU((!isForwardU) ? slidingFitResultU.GetGlobalMinLayerPosition() : slidingFitResultU.GetGlobalMaxLayerPosition());

        const CartesianVector vtxV((isForwardV) ? slidingFitResultV.GetGlobalMinLayerPosition() : slidingFitResultV.GetGlobalMaxLayerPosition());
        const CartesianVector endV((!isForwardV) ? slidingFitResultV.GetGlobalMinLayerPosition() : slidingFitResultV.GetGlobalMaxLayerPosition());

        const CartesianVector vtxW((isForwardW) ? slidingFitResultW.GetGlobalMinLayerPosition() : slidingFitResultW.GetGlobalMaxLayerPosition());
        const CartesianVector endW((!isForwardW) ? slidingFitResultW.GetGlobalMinLayerPosition() : slidingFitResultW.GetGlobalMaxLayerPosition());

        // Merge start and end positions (three views)
        const float halfLengthSquaredU(0.25 * (vtxU - endU).GetMagnitudeSquared());
        const float halfLengthSquaredV(0.25 * (vtxV - endV).GetMagnitudeSquared());
        const float halfLengthSquaredW(0.25 * (vtxW - endW).GetMagnitudeSquared());

        CartesianVector vtxMergedU(0.f, 0.f, 0.f), vtxMergedV(0.f, 0.f, 0.f), vtxMergedW(0.f, 0.f, 0.f);
        CartesianVector endMergedU(0.f, 0.f, 0.f), endMergedV(0.f, 0.f, 0.f), endMergedW(0.f, 0.f, 0.f);

        float vtxChi2(std::numeric_limits<float>::max());
        float endChi2(std::numeric_limits<float>::max());

        LArGeometryHelper::MergeThreePositions(
            this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, vtxU, vtxV, vtxW, vtxMergedU, vtxMergedV, vtxMergedW, vtxChi2);

        if (vtxChi2 > m_vertexChi2Cut)
            continue;

        if (((vtxMergedU - vtxU).GetMagnitudeSquared() > std::min(halfLengthSquaredU, (vtxMergedU - endU).GetMagnitudeSquared())) ||
            ((vtxMergedV - vtxV).GetMagnitudeSquared() > std::min(halfLengthSquaredV, (vtxMergedV - endV).GetMagnitudeSquared())) ||
            ((vtxMergedW - vtxW).GetMagnitudeSquared() > std::min(halfLengthSquaredW, (vtxMergedW - endW).GetMagnitudeSquared())))
            continue;

        LArGeometryHelper::MergeThreePositions(
            this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, endU, endV, endW, endMergedU, endMergedV, endMergedW, endChi2);

        if (endChi2 > m_vertexChi2Cut)
            continue;

        if (((endMergedU - endU).GetMagnitudeSquared() > std::min(halfLengthSquaredU, (endMergedU - vtxU).GetMagnitudeSquared())) ||
            ((endMergedV - endV).GetMagnitudeSquared() > std::min(halfLengthSquaredV, (endMergedV - vtxV).GetMagnitudeSquared())) ||
            ((endMergedW - endW).GetMagnitudeSquared() > std::min(halfLengthSquaredW, (endMergedW - vtxW).GetMagnitudeSquared())))
            continue;

        // Merge start and end positions (two views)
        float chi2(0.f);
        CartesianVector position3D(0.f, 0.f, 0.f);
        CartesianPointVector vtxList3D, endList3D;

        LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, vtxU, vtxV, position3D, chi2);
        vtxList3D.push_back(position3D);

        LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), TPC_VIEW_V, TPC_VIEW_W, vtxV, vtxW, position3D, chi2);
        vtxList3D.push_back(position3D);

        LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), TPC_VIEW_W, TPC_VIEW_U, vtxW, vtxU, position3D, chi2);
        vtxList3D.push_back(position3D);

        LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, endU, endV, position3D, chi2);
        endList3D.push_back(position3D);

        LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), TPC_VIEW_V, TPC_VIEW_W, endV, endW, position3D, chi2);
        endList3D.push_back(position3D);

        LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), TPC_VIEW_W, TPC_VIEW_U, endW, endU, position3D, chi2);
        endList3D.push_back(position3D);

        // Find best matched 3D trajactory
        for (CartesianPointVector::const_iterator iterI = vtxList3D.begin(), iterEndI = vtxList3D.end(); iterI != iterEndI; ++iterI)
        {
            const CartesianVector &vtxMerged3D(*iterI);

            for (CartesianPointVector::const_iterator iterJ = endList3D.begin(), iterEndJ = endList3D.end(); iterJ != iterEndJ; ++iterJ)
            {
                const CartesianVector &endMerged3D(*iterJ);

                TrackOverlapResult overlapResult;
                this->CalculateOverlapResult(slidingFitResultU, slidingFitResultV, slidingFitResultW, vtxMerged3D, endMerged3D, overlapResult);

                if (overlapResult.IsInitialized() && (overlapResult.GetNMatchedSamplingPoints() > 0) && (overlapResult > bestOverlapResult))
                {
                    bestOverlapResult = overlapResult;
                    longitudinalOverlapResult = LongitudinalOverlapResult(overlapResult, vtxChi2, endChi2);
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewLongitudinalTracksAlgorithm::CalculateOverlapResult(const TwoDSlidingFitResult &slidingFitResultU,
    const TwoDSlidingFitResult &slidingFitResultV, const TwoDSlidingFitResult &slidingFitResultW, const CartesianVector &vtxMerged3D,
    const CartesianVector &endMerged3D, TrackOverlapResult &overlapResult) const
{
    // Calculate start and end positions of linear trajectory
    const CartesianVector vtxMergedU(LArGeometryHelper::ProjectPosition(this->GetPandora(), vtxMerged3D, TPC_VIEW_U));
    const CartesianVector vtxMergedV(LArGeometryHelper::ProjectPosition(this->GetPandora(), vtxMerged3D, TPC_VIEW_V));
    const CartesianVector vtxMergedW(LArGeometryHelper::ProjectPosition(this->GetPandora(), vtxMerged3D, TPC_VIEW_W));

    const CartesianVector endMergedU(LArGeometryHelper::ProjectPosition(this->GetPandora(), endMerged3D, TPC_VIEW_U));
    const CartesianVector endMergedV(LArGeometryHelper::ProjectPosition(this->GetPandora(), endMerged3D, TPC_VIEW_V));
    const CartesianVector endMergedW(LArGeometryHelper::ProjectPosition(this->GetPandora(), endMerged3D, TPC_VIEW_W));

    const unsigned int nSamplingPoints = static_cast<unsigned int>((endMerged3D - vtxMerged3D).GetMagnitude() / m_samplingPitch);

    if (0 == nSamplingPoints)
        return;

    // Loop over sampling points and calculate track overlap result
    float deltaChi2(0.f), totalChi2(0.f);
    unsigned int nMatchedSamplingPoints(0);

    for (unsigned int n = 0; n < nSamplingPoints; ++n)
    {
        const float alpha((0.5f + static_cast<float>(n)) / static_cast<float>(nSamplingPoints));
        const CartesianVector linearU(vtxMergedU + (endMergedU - vtxMergedU) * alpha);
        const CartesianVector linearV(vtxMergedV + (endMergedV - vtxMergedV) * alpha);
        const CartesianVector linearW(vtxMergedW + (endMergedW - vtxMergedW) * alpha);

        CartesianVector posU(0.f, 0.f, 0.f), posV(0.f, 0.f, 0.f), posW(0.f, 0.f, 0.f);
        if ((STATUS_CODE_SUCCESS != slidingFitResultU.GetGlobalFitProjection(linearU, posU)) ||
            (STATUS_CODE_SUCCESS != slidingFitResultV.GetGlobalFitProjection(linearV, posV)) ||
            (STATUS_CODE_SUCCESS != slidingFitResultW.GetGlobalFitProjection(linearW, posW)))
        {
            continue;
        }

        CartesianVector mergedU(0.f, 0.f, 0.f), mergedV(0.f, 0.f, 0.f), mergedW(0.f, 0.f, 0.f);
        LArGeometryHelper::MergeThreePositions(this->GetPandora(), posU, posV, posW, mergedU, mergedV, mergedW, deltaChi2);

        if (deltaChi2 < m_reducedChi2Cut)
            ++nMatchedSamplingPoints;

        totalChi2 += deltaChi2;
    }

    if (nMatchedSamplingPoints > 0)
        overlapResult = TrackOverlapResult(nMatchedSamplingPoints, nSamplingPoints, totalChi2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeViewLongitudinalTracksAlgorithm::ExamineOverlapContainer()
{
    unsigned int repeatCounter(0);

    for (TensorToolVector::const_iterator iter = m_algorithmToolVector.begin(), iterEnd = m_algorithmToolVector.end(); iter != iterEnd;)
    {
        if ((*iter)->Run(this, this->GetMatchingControl().GetOverlapTensor()))
        {
            iter = m_algorithmToolVector.begin();

            if (++repeatCounter > m_nMaxTensorToolRepeats)
                break;
        }
        else
        {
            ++iter;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeViewLongitudinalTracksAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "TrackTools", algorithmToolVector));

    for (AlgorithmToolVector::const_iterator iter = algorithmToolVector.begin(), iterEnd = algorithmToolVector.end(); iter != iterEnd; ++iter)
    {
        LongitudinalTensorTool *const pLongitudinalTensorTool(dynamic_cast<LongitudinalTensorTool *>(*iter));

        if (!pLongitudinalTensorTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolVector.push_back(pLongitudinalTensorTool);
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NMaxTensorToolRepeats", m_nMaxTensorToolRepeats));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexChi2Cut", m_vertexChi2Cut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ReducedChi2Cut", m_reducedChi2Cut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SamplingPitch", m_samplingPitch));

    if (m_samplingPitch < std::numeric_limits<float>::epsilon())
        return STATUS_CODE_INVALID_PARAMETER;

    return BaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
