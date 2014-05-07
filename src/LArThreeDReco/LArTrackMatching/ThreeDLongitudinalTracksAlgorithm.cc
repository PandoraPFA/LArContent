/**
 *  @file   LArContent/src/LArThreeDReco/LArTrackMatching/ThreeDLongitudinalTracksAlgorithm.cc
 *
 *  @brief  Implementation of the three dimensional longitudinal tracks algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"

#include "LArThreeDReco/LArTrackMatching/ThreeDLongitudinalTracksAlgorithm.h"

using namespace pandora;

namespace lar
{

void ThreeDLongitudinalTracksAlgorithm::CalculateOverlapResult(Cluster *pClusterU, Cluster *pClusterV, Cluster *pClusterW)
{
    TwoDSlidingFitResult slidingFitResultU, slidingFitResultV, slidingFitResultW;
    LArClusterHelper::LArTwoDSlidingFit(pClusterU, 20, slidingFitResultU);
    LArClusterHelper::LArTwoDSlidingFit(pClusterV, 20, slidingFitResultV);
    LArClusterHelper::LArTwoDSlidingFit(pClusterW, 20, slidingFitResultW);

    LongitudinalOverlapResult bestOverlapResult(0, 1, m_reducedChi2Cut, 0.f, 0.f);

    for (unsigned int iPermutation = 0; iPermutation < 8; ++iPermutation)
    {
        const bool isForwardU((iPermutation >> 0) & 0x1);
        const bool isForwardV((iPermutation >> 1) & 0x1);
        const bool isForwardW((iPermutation >> 2) & 0x1);

        // Get 2D start and end positions for each sliding window fit
        const CartesianVector vtxU((isForwardU)  ? slidingFitResultU.GetGlobalMinLayerPosition() : slidingFitResultU.GetGlobalMaxLayerPosition());
        const CartesianVector endU((!isForwardU) ? slidingFitResultU.GetGlobalMinLayerPosition() : slidingFitResultU.GetGlobalMaxLayerPosition());

        const CartesianVector vtxV((isForwardV)  ? slidingFitResultV.GetGlobalMinLayerPosition() : slidingFitResultV.GetGlobalMaxLayerPosition());
        const CartesianVector endV((!isForwardV) ? slidingFitResultV.GetGlobalMinLayerPosition() : slidingFitResultV.GetGlobalMaxLayerPosition());

        const CartesianVector vtxW((isForwardW)  ? slidingFitResultW.GetGlobalMinLayerPosition() : slidingFitResultW.GetGlobalMaxLayerPosition());
        const CartesianVector endW((!isForwardW) ? slidingFitResultW.GetGlobalMinLayerPosition() : slidingFitResultW.GetGlobalMaxLayerPosition());

        float chi2(0.f);
        CartesianVector position3D(0.f,0.f,0.f);
        CartesianPointList vtxList3D, endList3D;

        // Calculate possible 3D start positions
        LArGeometryHelper::MergeTwoPositions3D(TPC_VIEW_U, TPC_VIEW_V, vtxU, vtxV, position3D, chi2);
        if (chi2 < m_vertexChi2Cut)
            vtxList3D.push_back(position3D);

        LArGeometryHelper::MergeTwoPositions3D(TPC_VIEW_V, TPC_VIEW_W, vtxV, vtxW, position3D, chi2);
        if (chi2 < m_vertexChi2Cut)
            vtxList3D.push_back(position3D);

        LArGeometryHelper::MergeTwoPositions3D(TPC_VIEW_W, TPC_VIEW_U, vtxW, vtxU, position3D, chi2);
        if (chi2 < m_vertexChi2Cut)
            vtxList3D.push_back(position3D);

        // Calculate possible 3D end positions
        LArGeometryHelper::MergeTwoPositions3D(TPC_VIEW_U, TPC_VIEW_V, endU, endV, position3D, chi2);
        if (chi2 < m_vertexChi2Cut)
            endList3D.push_back(position3D);

        LArGeometryHelper::MergeTwoPositions3D(TPC_VIEW_V, TPC_VIEW_W, endV, endW, position3D, chi2);
        if (chi2 < m_vertexChi2Cut)
            endList3D.push_back(position3D);

        LArGeometryHelper::MergeTwoPositions3D(TPC_VIEW_W, TPC_VIEW_U, endW, endU, position3D, chi2);
        if (chi2 < m_vertexChi2Cut)
            endList3D.push_back(position3D);

        // Find best matched 3D trajactory
        for (CartesianPointList::const_iterator iterI = vtxList3D.begin(), iterEndI = vtxList3D.end(); iterI != iterEndI; ++iterI)
        {
            const CartesianVector &vtxMerged3D(*iterI);

            const CartesianVector vtxMergedU(LArGeometryHelper::ProjectPosition(vtxMerged3D, TPC_VIEW_U));
            const CartesianVector vtxMergedV(LArGeometryHelper::ProjectPosition(vtxMerged3D, TPC_VIEW_V));
            const CartesianVector vtxMergedW(LArGeometryHelper::ProjectPosition(vtxMerged3D, TPC_VIEW_W));

            for (CartesianPointList::const_iterator iterJ = endList3D.begin(), iterEndJ = endList3D.end(); iterJ != iterEndJ; ++iterJ)
            {
                const CartesianVector &endMerged3D(*iterJ);

                const CartesianVector endMergedU(LArGeometryHelper::ProjectPosition(endMerged3D, TPC_VIEW_U));
                const CartesianVector endMergedV(LArGeometryHelper::ProjectPosition(endMerged3D, TPC_VIEW_V));
                const CartesianVector endMergedW(LArGeometryHelper::ProjectPosition(endMerged3D, TPC_VIEW_W));

                if ( ((endMergedU - vtxMergedU).GetCosOpeningAngle(endU - vtxU) < m_cosOpeningAngleCut) ||
                     ((endMergedV - vtxMergedV).GetCosOpeningAngle(endV - vtxV) < m_cosOpeningAngleCut) ||
                     ((endMergedW - vtxMergedW).GetCosOpeningAngle(endW - vtxW) < m_cosOpeningAngleCut) )
                {
                    continue;
                }

                if ( ((vtxMergedU - vtxU).GetMagnitudeSquared() > (vtxMergedU - endU).GetMagnitudeSquared()) ||
                     ((vtxMergedV - vtxV).GetMagnitudeSquared() > (vtxMergedV - endV).GetMagnitudeSquared()) ||
                     ((vtxMergedW - vtxW).GetMagnitudeSquared() > (vtxMergedW - endW).GetMagnitudeSquared()) ||
                     ((endMergedU - endU).GetMagnitudeSquared() > (endMergedU - vtxU).GetMagnitudeSquared()) ||
                     ((endMergedV - endV).GetMagnitudeSquared() > (endMergedV - vtxV).GetMagnitudeSquared()) ||
                     ((endMergedW - endW).GetMagnitudeSquared() > (endMergedW - vtxW).GetMagnitudeSquared()) ||
                     ((vtxMergedU - vtxU).GetMagnitudeSquared() > (endMergedU - vtxU).GetMagnitudeSquared()) ||
                     ((vtxMergedV - vtxV).GetMagnitudeSquared() > (endMergedV - vtxV).GetMagnitudeSquared()) ||
                     ((vtxMergedW - vtxW).GetMagnitudeSquared() > (endMergedW - vtxW).GetMagnitudeSquared()) ||
                     ((endMergedU - endU).GetMagnitudeSquared() > (vtxMergedU - endU).GetMagnitudeSquared()) ||
                     ((endMergedV - endV).GetMagnitudeSquared() > (vtxMergedV - endV).GetMagnitudeSquared()) ||
                     ((endMergedW - endW).GetMagnitudeSquared() > (vtxMergedW - endW).GetMagnitudeSquared()) )
                {
                    continue;
                }

                LongitudinalOverlapResult thisOverlapResult(0, 1, m_reducedChi2Cut, 0.f, 0.f);
                this->CalculateOverlapResult(slidingFitResultU, slidingFitResultV, slidingFitResultW, vtxMerged3D, endMerged3D, thisOverlapResult);

                if (thisOverlapResult.GetNMatchedSamplingPoints() > 0 && thisOverlapResult.GetReducedChi2() < bestOverlapResult.GetReducedChi2())
                    bestOverlapResult = thisOverlapResult;
            }
        }
    }

    if (bestOverlapResult.GetNMatchedSamplingPoints() > 0)
        m_overlapTensor.SetOverlapResult(pClusterU, pClusterV, pClusterW, bestOverlapResult);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDLongitudinalTracksAlgorithm::CalculateOverlapResult(const TwoDSlidingFitResult &slidingFitResultU, const TwoDSlidingFitResult &slidingFitResultV,
    const TwoDSlidingFitResult &slidingFitResultW, const CartesianVector &vtxMerged3D, const CartesianVector &endMerged3D, LongitudinalOverlapResult &overlapResult) const
{
    // Calculate start and end positions of linear trajectory
    const CartesianVector vtxMergedU(LArGeometryHelper::ProjectPosition(vtxMerged3D, TPC_VIEW_U));
    const CartesianVector vtxMergedV(LArGeometryHelper::ProjectPosition(vtxMerged3D, TPC_VIEW_V));
    const CartesianVector vtxMergedW(LArGeometryHelper::ProjectPosition(vtxMerged3D, TPC_VIEW_W));

    const CartesianVector endMergedU(LArGeometryHelper::ProjectPosition(endMerged3D, TPC_VIEW_U));
    const CartesianVector endMergedV(LArGeometryHelper::ProjectPosition(endMerged3D, TPC_VIEW_V));
    const CartesianVector endMergedW(LArGeometryHelper::ProjectPosition(endMerged3D, TPC_VIEW_W));

    const unsigned int nTotalSamplingPoints = static_cast<unsigned int>((endMerged3D - vtxMerged3D).GetMagnitude()/ m_samplingPitch);

    if (0 == nTotalSamplingPoints)
        return;

    float deltaChi2(0.f), innerChi2(0.f), outerChi2(0.f), totalChi2(0.f);
    unsigned int nSamplingPoints(0), nMatchedSamplingPoints(0);
    bool foundChi2(false);

    for (unsigned int n = 0; n < nTotalSamplingPoints; ++n)
    {
        const float alpha((0.5f + static_cast<float>(n)) / static_cast<float>(nTotalSamplingPoints));
        const CartesianVector linearU(vtxMergedU + (endMergedU - vtxMergedU) * alpha);
        const CartesianVector linearV(vtxMergedV + (endMergedV - vtxMergedV) * alpha);
        const CartesianVector linearW(vtxMergedW + (endMergedW - vtxMergedW) * alpha);

        try
        {
            CartesianVector posU(0.f,0.f,0.f), posV(0.f,0.f,0.f), posW(0.f,0.f,0.f);
            slidingFitResultU.GetGlobalFitProjection(linearU, posU);
            slidingFitResultV.GetGlobalFitProjection(linearV, posV);
            slidingFitResultW.GetGlobalFitProjection(linearW, posW);

            CartesianVector mergedU(0.f,0.f,0.f), mergedV(0.f,0.f,0.f), mergedW(0.f,0.f,0.f);
            LArGeometryHelper::MergeThreePositions(posU, posV, posW, mergedU, mergedV, mergedW, deltaChi2);

            if (deltaChi2 < m_reducedChi2Cut)
                ++nMatchedSamplingPoints;

            ++nSamplingPoints;

            if (!foundChi2)
                innerChi2 = deltaChi2;

            outerChi2 = deltaChi2;
            foundChi2 = true;

            totalChi2 += deltaChi2;
        }
        catch (StatusCodeException &)
        {
        }
    }

    if (nSamplingPoints > 0)
    {
        overlapResult = LongitudinalOverlapResult(nMatchedSamplingPoints, nSamplingPoints, totalChi2, innerChi2, outerChi2);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDLongitudinalTracksAlgorithm::ExamineTensorOld()
{
    while (true)
    {
        float bestReducedChi2(m_reducedChi2Cut);
        Cluster *pBestClusterU(NULL), *pBestClusterV(NULL), *pBestClusterW(NULL);

        for (TensorType::const_iterator iterU = m_overlapTensor.begin(), iterUEnd = m_overlapTensor.end(); iterU != iterUEnd; ++iterU)
        {
            for (TensorType::OverlapMatrix::const_iterator iterV = iterU->second.begin(), iterVEnd = iterU->second.end(); iterV != iterVEnd; ++iterV)
            {
                for (TensorType::OverlapList::const_iterator iterW = iterV->second.begin(), iterWEnd = iterV->second.end(); iterW != iterWEnd; ++iterW)
                {
                    const LongitudinalOverlapResult &overlapResult(iterW->second);

                    if (overlapResult.GetReducedChi2() < bestReducedChi2)
                    {
                        bestReducedChi2 = overlapResult.GetReducedChi2();
                        pBestClusterU = iterU->first;
                        pBestClusterV = iterV->first;
                        pBestClusterW = iterW->first;
                    }
                }
            }
        }

        if (!pBestClusterU || !pBestClusterV || !pBestClusterW)
            break;

        ProtoParticle protoParticle;
        protoParticle.m_clusterListU.insert(pBestClusterU);
        protoParticle.m_clusterListV.insert(pBestClusterV);
        protoParticle.m_clusterListW.insert(pBestClusterW);

        ProtoParticleVector protoParticleVector;
        protoParticleVector.push_back(protoParticle);

        this->CreateThreeDParticles(protoParticleVector);
        this->RemoveUnavailableTensorElements();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDLongitudinalTracksAlgorithm::ExamineTensor()
{
    unsigned int repeatCounter(0);

    for (LongitudinalTensorToolList::const_iterator iter = m_algorithmToolList.begin(), iterEnd = m_algorithmToolList.end(); iter != iterEnd; )
    {
        if ((*iter)->Run(this, m_overlapTensor))
        {
            iter = m_algorithmToolList.begin();

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
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDLongitudinalTracksAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    AlgorithmToolList algorithmToolList;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle,
        "TrackTools", algorithmToolList));

    for (AlgorithmToolList::const_iterator iter = algorithmToolList.begin(), iterEnd = algorithmToolList.end(); iter != iterEnd; ++iter)
    {
        LongitudinalTensorTool *pTensorManipulationTool(dynamic_cast<LongitudinalTensorTool*>(*iter));

        if (NULL == pTensorManipulationTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolList.push_back(pTensorManipulationTool);
    }

    m_nMaxTensorToolRepeats = 5000;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NMaxTensorToolRepeats", m_nMaxTensorToolRepeats));

    m_vertexChi2Cut = 5.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VertexChi2Cut", m_vertexChi2Cut));

    m_reducedChi2Cut = 5.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ReducedChi2Cut", m_reducedChi2Cut));

    m_cosOpeningAngleCut = 0.5;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CosOpeningAngleCut", m_cosOpeningAngleCut));

    m_samplingPitch = 1.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SamplingPitch", m_samplingPitch));

    if(m_samplingPitch < std::numeric_limits<float>::epsilon())
        return STATUS_CODE_INVALID_PARAMETER;

    return ThreeDBaseAlgorithm<LongitudinalOverlapResult>::ReadSettings(xmlHandle);
}

} // namespace lar
