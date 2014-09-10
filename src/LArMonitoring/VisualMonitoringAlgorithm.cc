/**
 *  @file   LArContent/src/LArMonitoring/VisualMonitoringAlgorithm.cc
 * 
 *  @brief  Implementation of the visual monitoring algorithm class 
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArMonitoring/VisualMonitoringAlgorithm.h"

#include <algorithm>
#include <string>

namespace lar_content
{

using namespace pandora;

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VisualMonitoringAlgorithm::Run()
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), m_showDetector, (m_detectorView.find("xz") != std::string::npos) ? DETECTOR_VIEW_XZ :
        (m_detectorView.find("xy") != std::string::npos) ? DETECTOR_VIEW_XY : DETECTOR_VIEW_DEFAULT, m_transparencyThresholdE, m_energyScaleThresholdE));

    // Show current mc particles
    if (m_showCurrentMCParticles)
    {
        this->VisualizeMCParticleList(std::string());
    }

    // Show specified lists of mc particles
    for (StringVector::const_iterator iter = m_mcParticleListNames.begin(), iterEnd = m_mcParticleListNames.end(); iter != iterEnd; ++iter)
    {
        this->VisualizeMCParticleList(*iter);
    }

    // Show current calo hit list
    if (m_showCurrentCaloHits)
    {
        this->VisualizeCaloHitList(std::string());
    }

    // Show specified lists of calo hits
    for (StringVector::const_iterator iter = m_caloHitListNames.begin(), iterEnd = m_caloHitListNames.end(); iter != iterEnd; ++iter)
    {
        this->VisualizeCaloHitList(*iter);
    }

    // Show current cluster list
    if (m_showCurrentClusters)
    {
        this->VisualizeClusterList(std::string());
    }

    // Show specified lists of clusters
    for (StringVector::const_iterator iter = m_clusterListNames.begin(), iterEnd = m_clusterListNames.end(); iter != iterEnd; ++iter)
    {
        this->VisualizeClusterList(*iter);
    }

    // Show current track list
    if (m_showCurrentTracks)
    {
        this->VisualizeTrackList(std::string());
    }

    // Show specified lists of tracks
    for (StringVector::const_iterator iter = m_trackListNames.begin(), iterEnd = m_trackListNames.end(); iter != iterEnd; ++iter)
    {
        this->VisualizeTrackList(*iter);
    }

    // Show current particle flow objects
    if (m_showCurrentPfos)
    {
        this->VisualizeParticleFlowList(std::string());
    }

    // Show specified lists of pfo
    for (StringVector::const_iterator iter = m_pfoListNames.begin(), iterEnd = m_pfoListNames.end(); iter != iterEnd; ++iter)
    {
        this->VisualizeParticleFlowList(*iter);
    }

    // Show current vertex objects
    if (m_showCurrentVertices)
    {
        this->VisualizeVertexList(std::string());
    }

    // Show specified lists of vertices
    for (StringVector::const_iterator iter = m_vertexListNames.begin(), iterEnd = m_vertexListNames.end(); iter != iterEnd; ++iter)
    {
        this->VisualizeVertexList(*iter);
    }

    // Finally, display the event and pause application
    if (m_displayEvent)
    {
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VisualMonitoringAlgorithm::VisualizeMCParticleList(const std::string &listName) const
{
    const MCParticleList *pMCParticleList = NULL;

    if (listName.empty())
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetCurrentList(*this, pMCParticleList))
        {
            std::cout << "VisualMonitoringAlgorithm: mc particle list unavailable." << std::endl;
            return;
        }
    }
    else
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, listName, pMCParticleList))
        {
            std::cout << "VisualMonitoringAlgorithm: mc particle list unavailable." << std::endl;
            return;
        }
    }

    PANDORA_MONITORING_API(VisualizeMCParticles(this->GetPandora(), pMCParticleList, listName.empty() ? "currentMCParticles" : listName.c_str(),
        AUTO, &m_particleSuppressionMap));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VisualMonitoringAlgorithm::VisualizeCaloHitList(const std::string &listName) const
{
    const CaloHitList *pCaloHitList = NULL;

    if (listName.empty())
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetCurrentList(*this, pCaloHitList))
        {
            std::cout << "VisualMonitoringAlgorithm: current calo hit list unavailable." << std::endl;
            return;
        }
    }
    else
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, listName, pCaloHitList))
        {
            std::cout << "VisualMonitoringAlgorithm: calo hit list " << listName << " unavailable." << std::endl;
            return;
        }
    }

    CaloHitList caloHitList(*pCaloHitList);

    // Filter calo hit list
    for (CaloHitList::const_iterator hitIter = caloHitList.begin(), hitIterEnd = caloHitList.end(); hitIter != hitIterEnd; )
    {
        if (((*hitIter)->GetElectromagneticEnergy() < m_thresholdEnergy))
        {
            caloHitList.erase(hitIter++);
        }
        else if (m_showOnlyAvailable && !PandoraContentApi::IsAvailable(*this, *hitIter))
        {
            caloHitList.erase(hitIter++);
        }
        else
        {
            hitIter++;
        }
    }

    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitList, listName.empty() ? "currentCaloHits" : listName.c_str(),
        (m_hitColors.find("energy") != std::string::npos ? AUTOENERGY : GRAY)));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VisualMonitoringAlgorithm::VisualizeTrackList(const std::string &listName) const
{
    const TrackList *pTrackList = NULL;

    if (listName.empty())
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetCurrentList(*this, pTrackList))
        {
            std::cout << "VisualMonitoringAlgorithm: current track list unavailable." << std::endl;
            return;
        }
    }
    else
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, listName, pTrackList))
        {
            std::cout << "VisualMonitoringAlgorithm: track list " << listName << " unavailable." << std::endl;
            return;
        }
    }

    // Filter track list
    TrackList trackList;

    for (TrackList::const_iterator iter = pTrackList->begin(), iterEnd = pTrackList->end(); iter != iterEnd; ++iter)
    {
        Track *pTrack = *iter;

        if (!m_showOnlyAvailable || PandoraContentApi::IsAvailable(*this, pTrack))
        {
            trackList.insert(pTrack);
        }
    }

    PANDORA_MONITORING_API(VisualizeTracks(this->GetPandora(), &trackList, listName.empty() ? "currentTracks" : listName.c_str(), GRAY));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VisualMonitoringAlgorithm::VisualizeClusterList(const std::string &listName) const
{
    const ClusterList *pClusterList = NULL;

    if (listName.empty())
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetCurrentList(*this, pClusterList))
        {
            std::cout << "VisualMonitoringAlgorithm: current cluster list unavailable." << std::endl;
            return;
        }
    }
    else
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, listName, pClusterList))
        {
            std::cout << "VisualMonitoringAlgorithm: cluster list " << listName << " unavailable." << std::endl;
            return;
        }
    }

    // Filter cluster list
    ClusterList clusterList;

    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (!m_showOnlyAvailable || PandoraContentApi::IsAvailable(*this, pCluster))
        {
            clusterList.insert(pCluster);
        }
    }

    PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clusterList, listName.empty() ? "currentClusters" : listName.c_str(),
        (m_hitColors.find("particleid") != std::string::npos) ? AUTOID :
        (m_hitColors.find("iterate") != std::string::npos ? AUTOITER :
        (m_hitColors.find("energy") != std::string::npos ? AUTOENERGY :
        AUTO))));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VisualMonitoringAlgorithm::VisualizeParticleFlowList(const std::string &listName) const
{
    const PfoList *pPfoList = NULL;

    if (listName.empty())
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetCurrentList(*this, pPfoList))
        {
            std::cout << "VisualMonitoringAlgorithm: current pfo list unavailable." << std::endl;
            return;
        }
    }
    else
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, listName, pPfoList))
        {
            std::cout << "VisualMonitoringAlgorithm: pfo list " << listName << " unavailable." << std::endl;
            return;
        }
    }

    PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), pPfoList, listName.empty() ? "currentPfos" : listName.c_str(),
        (m_hitColors.find("particleid") != std::string::npos) ? AUTOID :
        (m_hitColors.find("iterate") != std::string::npos ? AUTOITER :
        (m_hitColors.find("energy") != std::string::npos ? AUTOENERGY :
        AUTO)), m_showPfoVertices, m_showPfoHierarchy));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VisualMonitoringAlgorithm::VisualizeVertexList(const std::string &listName) const
{
    const VertexList *pVertexList = NULL;

    if (listName.empty())
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetCurrentList(*this, pVertexList))
        {
            std::cout << "VisualMonitoringAlgorithm: current vertex list unavailable." << std::endl;
            return;
        }
    }
    else
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, listName, pVertexList))
        {
            std::cout << "VisualMonitoringAlgorithm: vertex list " << listName << " unavailable." << std::endl;
            return;
        }
    }

    PANDORA_MONITORING_API(VisualizeVertices(this->GetPandora(), pVertexList, listName.empty() ? "currentVertices" : listName.c_str(), AUTO));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VisualMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_showCurrentMCParticles = false;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowCurrentMCParticles", m_showCurrentMCParticles));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "MCParticleListNames", m_mcParticleListNames));

    m_showCurrentCaloHits = false;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowCurrentCaloHits", m_showCurrentCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "CaloHitListNames", m_caloHitListNames));

    m_showCurrentTracks = false;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowCurrentTracks", m_showCurrentTracks));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "TrackListNames", m_trackListNames));

    m_showCurrentClusters = true;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowCurrentClusters", m_showCurrentClusters));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "ClusterListNames", m_clusterListNames));

    m_showCurrentPfos = true;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowCurrentPfos", m_showCurrentPfos));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "PfoListNames", m_pfoListNames));

    m_showCurrentVertices = false;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowCurrentVertices", m_showCurrentVertices));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "VertexListNames", m_vertexListNames));

    m_displayEvent = true;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DisplayEvent", m_displayEvent));

    m_showDetector = true;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowDetector", m_showDetector));

    m_detectorView = "xz";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DetectorView", m_detectorView));
    std::transform(m_detectorView.begin(), m_detectorView.end(), m_detectorView.begin(), ::tolower);

    m_showOnlyAvailable = false;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowOnlyAvailable", m_showOnlyAvailable));

    m_hitColors = "pfo";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "HitColors", m_hitColors));
    std::transform(m_hitColors.begin(), m_hitColors.end(), m_hitColors.begin(), ::tolower);

    m_thresholdEnergy = -1.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ThresholdEnergy", m_thresholdEnergy));

    m_transparencyThresholdE = -1.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TransparencyThresholdE", m_transparencyThresholdE));

    m_energyScaleThresholdE = 1.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "EnergyScaleThresholdE", m_energyScaleThresholdE));

    m_showPfoVertices = true;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowPfoVertices", m_showPfoVertices));

    m_showPfoHierarchy = true;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowPfoHierarchy", m_showPfoHierarchy));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "SuppressMCParticles", m_suppressMCParticles));

    for (StringVector::iterator iter = m_suppressMCParticles.begin(), iterEnd = m_suppressMCParticles.end(); iter != iterEnd; ++iter)
    {
        const std::string pdgEnergy(*iter);
        StringVector pdgEnergySeparated;
        const std::string delimiter = ":";
        XmlHelper::TokenizeString(pdgEnergy, pdgEnergySeparated, delimiter);

        if (pdgEnergySeparated.size() != 2)
            return STATUS_CODE_INVALID_PARAMETER;

        int pdgCode(0);
        float energy(0.f);

        if (!StringToType(pdgEnergySeparated.at(0), pdgCode) || !StringToType(pdgEnergySeparated.at(1), energy))
            return STATUS_CODE_INVALID_PARAMETER;

        m_particleSuppressionMap.insert(PdgCodeToEnergyMap::value_type(pdgCode, energy));
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
