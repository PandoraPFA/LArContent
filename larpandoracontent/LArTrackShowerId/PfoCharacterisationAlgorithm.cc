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
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArTrackShowerId/ClusterCharacterisationAlgorithm.h"
#include "larpandoracontent/LArTrackShowerId/PfoCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

PfoCharacterisationAlgorithm::PfoCharacterisationAlgorithm() :
    m_updateClusterIds(true),
    m_postBranchAddition(false),
    m_minTrackLikeViews(2),
    m_slidingFitWindow(5),
    m_slidingShowerFitWindow(10),
    m_maxShowerLengthCut(80.f),
    m_dTdLWidthRatioCut(0.045f),
    m_vertexDistanceRatioCut(0.6f),
    m_showerWidthRatioCut(0.2f)
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
    ClusterList twoDClusterList;
    LArPfoHelper::GetTwoDClusterList(pPfo, twoDClusterList);

    unsigned int nTrackLikeViews(0);
    HitTypeSet hitTypeSet;

    for (const Cluster *const pCluster : twoDClusterList)
    {
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));

        if (!hitTypeSet.insert(hitType).second)
            continue;

        if (this->IsClearTrack(pCluster))
            ++nTrackLikeViews;

        if (nTrackLikeViews >= m_minTrackLikeViews)
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool PfoCharacterisationAlgorithm::IsClearTrack(const Cluster *const pCluster) const
{
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

    if (straightLineLength > m_maxShowerLengthCut)
        return true;

    if ((dTdLMax - dTdLMin) / straightLineLength > m_dTdLWidthRatioCut)
        return false;

    const float vertexDistance(ClusterCharacterisationAlgorithm::GetVertexDistance(this, pCluster));

    if ((vertexDistance > std::numeric_limits<float>::epsilon()) && ((vertexDistance / straightLineLength) > m_vertexDistanceRatioCut))
        return false;

    const float showerFitWidth(ClusterCharacterisationAlgorithm::GetShowerFitWidth(this, pCluster, m_slidingShowerFitWindow));

    if ((showerFitWidth < std::numeric_limits<float>::epsilon()) || ((showerFitWidth / straightLineLength) > m_showerWidthRatioCut))
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

    // Allow change in default values via a single xml tag, can subsequently override all individual values below, if required
    if (m_postBranchAddition)
    {
        m_maxShowerLengthCut = 80.f;
        m_dTdLWidthRatioCut = 0.03f;
        m_vertexDistanceRatioCut = 1.f;
        m_showerWidthRatioCut = 0.3f;
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinTrackLikeViews", m_minTrackLikeViews));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingShowerFitWindow", m_slidingShowerFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxShowerLengthCut", m_maxShowerLengthCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DTDLWidthRatioCut", m_dTdLWidthRatioCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VertexDistanceRatioCut", m_vertexDistanceRatioCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowerWidthRatioCut", m_showerWidthRatioCut));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
