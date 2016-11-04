/**
 *  @file   larpandoracontent/LArTrackShowerId/ClusterCharacterisationAlgorithm.cc
 *
 *  @brief  Implementation of the cluster characterisation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingShowerFitResult.h"

#include "larpandoracontent/LArTrackShowerId/ClusterCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ClusterCharacterisationAlgorithm::ClusterCharacterisationAlgorithm() :
    m_overwriteExistingId(false),
    m_useUnavailableClusters(false),
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

StatusCode ClusterCharacterisationAlgorithm::Run()
{
    for (const std::string &clusterListName : m_inputClusterListNames)
    {
        const ClusterList *pClusterList = NULL;
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));

        if (!pClusterList || pClusterList->empty())
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "ClusterCharacterisationAlgorithm: unable to find cluster list " << clusterListName << std::endl;

            continue;
        }

        for (const Cluster *const pCluster : *pClusterList)
        {
            if (!m_overwriteExistingId && (UNKNOWN_PARTICLE_TYPE != pCluster->GetParticleId()))
                continue;

            if (!m_useUnavailableClusters && !PandoraContentApi::IsAvailable(*this, pCluster))
                continue;

            PandoraContentApi::Cluster::Metadata metadata;

            if (this->IsClearTrack(pCluster))
            {
                metadata.m_particleId = MU_MINUS;
            }
            else
            {
                metadata.m_particleId = E_MINUS;
            }

            if (pCluster->GetParticleId() != metadata.m_particleId.Get())
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::AlterMetadata(*this, pCluster, metadata));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ClusterCharacterisationAlgorithm::IsClearTrack(const Cluster *const pCluster) const
{
    if (pCluster->GetNCaloHits() < m_minCaloHitsCut)
        return false;

    // Quantities related to sliding linear fit
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

    if (integratedPathLength / straightLineLength > m_pathLengthRatioCut)
        return false;

    if ((rTMax - rTMin) / straightLineLength > m_rTWidthRatioCut)
        return false;

    // Distance to interaction vertex
    const VertexList *pVertexList = nullptr;
    (void) PandoraContentApi::GetCurrentList(*this, pVertexList);

    if (pVertexList && (pVertexList->size() == 1) && (VERTEX_3D == pVertexList->front()->GetVertexType()))
    {
        const Vertex *const pVertex(pVertexList->front());
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));

        const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), hitType));
        const float vertexDistance(LArClusterHelper::GetClosestDistance(vertexPosition2D, pCluster));

        if (vertexDistance / straightLineLength > m_vertexDistanceRatioCut)
            return false;
    }

    // Shower fit width
    float showerFitWidth(-1.f);

    try
    {
        const TwoDSlidingShowerFitResult showerFitResult(pCluster, m_slidingShowerFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const LayerFitResultMap &layerFitResultMapS(showerFitResult.GetShowerFitResult().GetLayerFitResultMap());
        const LayerFitResultMap &layerFitResultMapP(showerFitResult.GetPositiveEdgeFitResult().GetLayerFitResultMap());
        const LayerFitResultMap &layerFitResultMapN(showerFitResult.GetNegativeEdgeFitResult().GetLayerFitResultMap());

        if (!layerFitResultMapS.empty())
        {
            showerFitWidth = 0.f;

            for (LayerFitResultMap::const_iterator iterS = layerFitResultMapS.begin(); iterS != layerFitResultMapS.end(); ++iterS)
            {
                LayerFitResultMap::const_iterator iterP = layerFitResultMapP.find(iterS->first);
                LayerFitResultMap::const_iterator iterN = layerFitResultMapN.find(iterS->first);

                if ((layerFitResultMapP.end() != iterP) && (layerFitResultMapN.end() != iterN))
                    showerFitWidth += std::fabs(iterP->second.GetFitT() - iterN->second.GetFitT());
            }
        }
    }
    catch (const StatusCodeException &)
    {
    }

    if (showerFitWidth < std::numeric_limits<float>::epsilon())
        return false;

    if (showerFitWidth / straightLineLength > m_showerWidthRatioCut)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "InputClusterListNames", m_inputClusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OverwriteExistingId", m_overwriteExistingId));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseUnavailableClusters", m_useUnavailableClusters));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseUnavailableClusters", m_useUnavailableClusters));

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

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
