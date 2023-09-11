/**
 *  @file   larpandoracontent/LArTrackShowerId/CutPfoCharacterisationAlgorithm.cc
 *
 *  @brief  Implementation of the cut based pfo characterisation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArTrackShowerId/CutClusterCharacterisationAlgorithm.h"
#include "larpandoracontent/LArTrackShowerId/CutPfoCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CutPfoCharacterisationAlgorithm::CutPfoCharacterisationAlgorithm() :
    m_postBranchAddition(false),
    m_slidingFitWindow(5),
    m_slidingShowerFitWindow(10),
    m_maxShowerLengthCut(80.f),
    m_dTdLWidthRatioCut(0.045f),
    m_vertexDistanceRatioCut(0.6f),
    m_showerWidthRatioCut(0.2f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CutPfoCharacterisationAlgorithm::IsClearTrack(const Cluster *const pCluster) const
{
    float straightLineLength(-1.f);
    float dTdLMin(+std::numeric_limits<float>::max()), dTdLMax(-std::numeric_limits<float>::max());

    try
    {
        const HitType view{LArClusterHelper::GetClusterHitType(pCluster)};
        const TwoDSlidingFitResult slidingFitResult(pCluster, m_slidingFitWindow, LArGeometryHelper::GetWirePitch(this->GetPandora(), view));
        straightLineLength = (slidingFitResult.GetGlobalMaxLayerPosition() - slidingFitResult.GetGlobalMinLayerPosition()).GetMagnitude();

        for (const auto &mapEntry : slidingFitResult.GetLayerFitResultMap())
        {
            dTdLMin = std::min(dTdLMin, static_cast<float>(mapEntry.second.GetGradient()));
            dTdLMax = std::max(dTdLMax, static_cast<float>(mapEntry.second.GetGradient()));
        }
    }
    catch (const StatusCodeException &)
    {
    }

    if (straightLineLength < std::numeric_limits<float>::epsilon())
        return false;

    if (straightLineLength > m_maxShowerLengthCut)
        return true;

    if ((dTdLMax - dTdLMin) / straightLineLength > m_dTdLWidthRatioCut)
        return false;

    const float vertexDistance(CutClusterCharacterisationAlgorithm::GetVertexDistance(this, pCluster));

    if ((vertexDistance > std::numeric_limits<float>::epsilon()) && ((vertexDistance / straightLineLength) > m_vertexDistanceRatioCut))
        return false;

    const float showerFitWidth(CutClusterCharacterisationAlgorithm::GetShowerFitWidth(this, pCluster, m_slidingShowerFitWindow));

    if ((showerFitWidth < std::numeric_limits<float>::epsilon()) || ((showerFitWidth / straightLineLength) > m_showerWidthRatioCut))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CutPfoCharacterisationAlgorithm::IsClearTrack(const pandora::ParticleFlowObject *const /*pPfo*/) const
{
    throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CutPfoCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PostBranchAddition", m_postBranchAddition));

    // Allow change in default values via a single xml tag, can subsequently override all individual values below, if required
    if (m_postBranchAddition)
    {
        m_maxShowerLengthCut = 80.f;
        m_dTdLWidthRatioCut = 0.03f;
        m_vertexDistanceRatioCut = 1.f;
        m_showerWidthRatioCut = 0.3f;
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingShowerFitWindow", m_slidingShowerFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxShowerLengthCut", m_maxShowerLengthCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "DTDLWidthRatioCut", m_dTdLWidthRatioCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexDistanceRatioCut", m_vertexDistanceRatioCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ShowerWidthRatioCut", m_showerWidthRatioCut));

    return PfoCharacterisationBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
