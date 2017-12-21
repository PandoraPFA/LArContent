/**
 *  @file   larpandoracontent/LArTrackShowerId/CutClusterCharacterisationAlgorithm.cc
 *
 *  @brief  Implementation of the cut based cluster characterisation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingShowerFitResult.h"

#include "larpandoracontent/LArTrackShowerId/CutClusterCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CutClusterCharacterisationAlgorithm::CutClusterCharacterisationAlgorithm() :
    m_slidingFitWindow(5),
    m_slidingShowerFitWindow(10),
    m_minCaloHitsCut(6),
    m_maxShowerLengthCut(80.f),
    m_pathLengthRatioCut(1.005f),
    m_rTWidthRatioCut(0.05f),
    m_vertexDistanceRatioCut(0.5f),
    m_showerWidthRatioCut(0.35f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

float CutClusterCharacterisationAlgorithm::GetVertexDistance(const Algorithm *const pAlgorithm, const Cluster *const pCluster)
{
    const VertexList *pVertexList = nullptr;
    (void) PandoraContentApi::GetCurrentList(*pAlgorithm, pVertexList);

    if (!pVertexList || (pVertexList->size() != 1) || (VERTEX_3D != pVertexList->front()->GetVertexType()))
        return -1.f;

    const Vertex *const pVertex(pVertexList->front());
    const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));

    const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), pVertex->GetPosition(), hitType));
    return LArClusterHelper::GetClosestDistance(vertexPosition2D, pCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float CutClusterCharacterisationAlgorithm::GetShowerFitWidth(const Algorithm *const pAlgorithm, const Cluster *const pCluster,
    const unsigned int showerFitWindow)
{
    try
    {
        const TwoDSlidingShowerFitResult showerFitResult(pCluster, showerFitWindow, LArGeometryHelper::GetWireZPitch(pAlgorithm->GetPandora()));
        const LayerFitResultMap &layerFitResultMapS(showerFitResult.GetShowerFitResult().GetLayerFitResultMap());
        const LayerFitResultMap &layerFitResultMapP(showerFitResult.GetPositiveEdgeFitResult().GetLayerFitResultMap());
        const LayerFitResultMap &layerFitResultMapN(showerFitResult.GetNegativeEdgeFitResult().GetLayerFitResultMap());

        if (!layerFitResultMapS.empty())
        {
            float showerFitWidth(0.f);

            for (const auto &mapEntryS : layerFitResultMapS)
            {
                LayerFitResultMap::const_iterator iterP = layerFitResultMapP.find(mapEntryS.first);
                LayerFitResultMap::const_iterator iterN = layerFitResultMapN.find(mapEntryS.first);

                if ((layerFitResultMapP.end() != iterP) && (layerFitResultMapN.end() != iterN))
                    showerFitWidth += std::fabs(iterP->second.GetFitT() - iterN->second.GetFitT());
            }

            return showerFitWidth;
        }
    }
    catch (const StatusCodeException &)
    {
    }

    return -1.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CutClusterCharacterisationAlgorithm::IsClearTrack(const Cluster *const pCluster) const
{
    if (pCluster->GetNCaloHits() < m_minCaloHitsCut)
        return false;

    float straightLineLength(-1.f), integratedPathLength(-1.f);
    float rTMin(+std::numeric_limits<float>::max()), rTMax(-std::numeric_limits<float>::max());

    try
    {
        const TwoDSlidingFitResult slidingFitResult(pCluster, m_slidingFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const CartesianVector globalMinLayerPosition(slidingFitResult.GetGlobalMinLayerPosition());
        straightLineLength = (slidingFitResult.GetGlobalMaxLayerPosition() - globalMinLayerPosition).GetMagnitude();

        integratedPathLength = 0.f;
        CartesianVector previousFitPosition(globalMinLayerPosition);

        for (const auto &mapEntry : slidingFitResult.GetLayerFitResultMap())
        {
            rTMin = std::min(rTMin, static_cast<float>(mapEntry.second.GetFitT()));
            rTMax = std::max(rTMax, static_cast<float>(mapEntry.second.GetFitT()));

            CartesianVector thisFitPosition(0.f, 0.f, 0.f);
            slidingFitResult.GetGlobalPosition(mapEntry.second.GetL(), mapEntry.second.GetFitT(), thisFitPosition);
            integratedPathLength += (thisFitPosition - previousFitPosition).GetMagnitude();
            previousFitPosition = thisFitPosition;
        }
    }
    catch (const StatusCodeException &)
    {
    }

    if (straightLineLength < std::numeric_limits<float>::epsilon())
        return false;

    if (straightLineLength > m_maxShowerLengthCut)
        return true;

    if ((integratedPathLength < std::numeric_limits<float>::epsilon()) || (integratedPathLength / straightLineLength > m_pathLengthRatioCut))
        return false;

    if ((rTMax - rTMin) / straightLineLength > m_rTWidthRatioCut)
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

StatusCode CutClusterCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingShowerFitWindow", m_slidingShowerFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCaloHitsCut", m_minCaloHitsCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxShowerLengthCut", m_maxShowerLengthCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PathLengthRatioCut", m_pathLengthRatioCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "RTWidthRatioCut", m_rTWidthRatioCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VertexDistanceRatioCut", m_vertexDistanceRatioCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowerWidthRatioCut", m_showerWidthRatioCut));

    return ClusterCharacterisationBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
