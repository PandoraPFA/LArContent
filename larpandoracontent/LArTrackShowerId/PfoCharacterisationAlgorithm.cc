/**
 *  @file   larpandoracontent/LArTrackShowerId/PfoCharacterisationAlgorithm.cc
 * 
 *  @brief  Implementation of the pfo characterisation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingShowerFitResult.h"

#include "larpandoracontent/LArTrackShowerId/PfoCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

PfoCharacterisationAlgorithm::PfoCharacterisationAlgorithm() :
    m_updateClusterIds(true),
    m_postBranchAddition(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoCharacterisationAlgorithm::Run()
{
    PfoList tracksToShowers, showersToTracks;

    for (const std::string &pfoListName : m_inputPfoListNames)
    {
        const PfoList *pPfoList(nullptr);
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, pfoListName, pPfoList));

        if (!pPfoList || pPfoList->empty())
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "PfoCharacterisationAlgorithm: unable to find pfo list " << pfoListName << std::endl;

            continue;
        }

        for (const ParticleFlowObject *const pPfo : *pPfoList)
        {
            PandoraContentApi::ParticleFlowObject::Metadata pfoMetadata;

            if (this->IsClearTrack(pPfo))
            {
                pfoMetadata.m_particleId = MU_MINUS;

                if (m_showerPfoListName == pfoListName)
                    showersToTracks.push_back(pPfo);
            }
            else
            {
                pfoMetadata.m_particleId = E_MINUS;

                if (m_trackPfoListName == pfoListName)
                    tracksToShowers.push_back(pPfo);
            }

            if (pPfo->GetParticleId() != pfoMetadata.m_particleId.Get())
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, pfoMetadata));

            if (!m_updateClusterIds)
                continue;

            ClusterList twoDClusterList;
            LArPfoHelper::GetTwoDClusterList(pPfo, twoDClusterList);

            for (const Cluster *const pCluster : twoDClusterList)
            {
                if (pCluster->GetParticleId() == pfoMetadata.m_particleId.Get())
                    continue;

                PandoraContentApi::Cluster::Metadata clusterMetadata;
                clusterMetadata.m_particleId = pfoMetadata.m_particleId.Get();
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::AlterMetadata(*this, pCluster, clusterMetadata));
            }
        }
    }

    if (!tracksToShowers.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_trackPfoListName, m_showerPfoListName, tracksToShowers));

    if (!showersToTracks.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_showerPfoListName, m_trackPfoListName, showersToTracks));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool PfoCharacterisationAlgorithm::IsClearTrack(const ParticleFlowObject *const pPfo) const
{
    // TODO, methodise
    const Cluster *pCluster(nullptr);

    // Prefer w cluster, if available, then u, followed by v
    for (const HitType hitType : {TPC_VIEW_W, TPC_VIEW_U, TPC_VIEW_V})
    {
        ClusterList clusterList;
        LArPfoHelper::GetClusters(pPfo, hitType, clusterList);

        if (clusterList.size() > 1)
            return false;

        if (!clusterList.empty())
        {
            pCluster = clusterList.front();
            break;
        }
    }
    
    if (!pCluster)
        return false;

    const float straightLineLengthCut(m_postBranchAddition ? 80.f : 40.f);
    const float dTdLWidthCut(m_postBranchAddition ? 0.030f : 0.045f);
    const float vertexDistanceCut(m_postBranchAddition ? 1.0f : 0.6f);
    const float showerFitWidthCut(m_postBranchAddition ? 0.3f : 0.2f);

    // Quantities related to sliding linear fit
    float straightLineLength(-1.f);
    float dTdLMin(+std::numeric_limits<float>::max()), dTdLMax(-std::numeric_limits<float>::max());

    try
    {
        const TwoDSlidingFitResult slidingFitResult(pCluster, 5, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
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

    if (straightLineLength > straightLineLengthCut)
        return true;

    const float dTdLWidth(dTdLMax - dTdLMin);

    if (dTdLWidth / straightLineLength > dTdLWidthCut)
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

        if (vertexDistance / straightLineLength > vertexDistanceCut)
            return false;
    }

    // Shower fit width
    float showerFitWidth(-1.f);

    try
    {
        const TwoDSlidingShowerFitResult showerFitResult(pCluster, 10, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
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

    if (showerFitWidth / straightLineLength > showerFitWidthCut)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrackPfoListName", m_trackPfoListName));
    m_inputPfoListNames.push_back(m_trackPfoListName);

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ShowerPfoListName", m_showerPfoListName));
    m_inputPfoListNames.push_back(m_showerPfoListName);

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UpdateClusterIds", m_updateClusterIds));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PostBranchAddition", m_postBranchAddition));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
