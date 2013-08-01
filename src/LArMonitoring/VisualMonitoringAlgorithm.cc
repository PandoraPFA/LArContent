/**
 *  @file   LArContent/src/LArMonitoring/VisualMonitoringAlgorithm.cc
 * 
 *  @brief  Implementation of the visual monitoring algorithm class 
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArVertexHelper.h"

#include "LArMonitoring/VisualMonitoringAlgorithm.h"

#include <algorithm>
#include <string>

namespace lar
{

using namespace pandora;

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VisualMonitoringAlgorithm::Run()
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(m_blackBackground, m_showDetector, m_transparencyThresholdE, m_energyScaleThresholdE));

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

    // Show current vertex
    if (m_showCurrentVertex)
    {
        this->VisualizeVertex(std::string());
    }

    // Show specified lists of vertices
    for (StringVector::const_iterator iter = m_vertexNames.begin(), iterEnd = m_vertexNames.end(); iter != iterEnd; ++iter)
    {
        this->VisualizeVertex(*iter);
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

    // Finally, display the event and pause application
    if (m_displayEvent)
    {
        PANDORA_MONITORING_API(ViewEvent());
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VisualMonitoringAlgorithm::VisualizeMCParticleList(const std::string &listName) const
{
    const MCParticleList *pMCParticleList = NULL;

    if (listName.empty())
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetCurrentMCParticleList(*this, pMCParticleList))
        {
            std::cout << "VisualMonitoringAlgorithm: mc particle list unavailable." << std::endl;
            return;
        }
    }
    else
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetMCParticleList(*this, listName, pMCParticleList))
        {
            std::cout << "VisualMonitoringAlgorithm: mc particle list unavailable." << std::endl;
            return;
        }
    }

    PANDORA_MONITORING_API(VisualizeMCParticles(pMCParticleList, listName.empty() ? "currentMCParticles" : listName.c_str(),
        AUTO, &m_particleSuppressionMap));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VisualMonitoringAlgorithm::VisualizeVertex(const std::string &vertexName) const
{
    try
    {
        if (vertexName.empty())
        {
            const CartesianVector &currentVertex(LArVertexHelper::GetCurrentVertex());
            PANDORA_MONITORING_API(AddMarkerToVisualization(&currentVertex, "CurrentVertex", AUTO, 1.));
        }
        else
        {
            const CartesianVector &vertex(LArVertexHelper::GetVertex(vertexName));
            PANDORA_MONITORING_API(AddMarkerToVisualization(&vertex, vertexName, AUTO, 1.));
        }
    }
    catch (StatusCodeException &)
    {
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VisualMonitoringAlgorithm::VisualizeCaloHitList(const std::string &listName) const
{
    const CaloHitList *pCaloHitList = NULL;

    if (listName.empty())
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetCurrentCaloHitList(*this, pCaloHitList))
        {
            std::cout << "VisualMonitoringAlgorithm: current calo hit list unavailable." << std::endl;
            return;
        }
    }
    else
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetCaloHitList(*this, listName, pCaloHitList))
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
        else if (m_showOnlyAvailable && !PandoraContentApi::IsCaloHitAvailable(*this, *hitIter))
        {
            caloHitList.erase(hitIter++);
        }
        else
        {
            hitIter++;
        }
    }

    PANDORA_MONITORING_API(VisualizeCaloHits(&caloHitList, listName.empty() ? "currentCaloHits" : listName.c_str(),
        (m_hitColors.find("energy") != std::string::npos ? AUTOENERGY : GRAY)));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VisualMonitoringAlgorithm::VisualizeClusterList(const std::string &listName) const
{
    const ClusterList *pClusterList = NULL;

    if (listName.empty())
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetCurrentClusterList(*this, pClusterList))
        {
            std::cout << "VisualMonitoringAlgorithm: current cluster list unavailable." << std::endl;
            return;
        }
    }
    else
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetClusterList(*this, listName, pClusterList))
        {
            std::cout << "VisualMonitoringAlgorithm: cluster list " << listName << " unavailable." << std::endl;
            return;
        }
    }

    ClusterList clusterList(*pClusterList);

    // Filter cluster list
    for (ClusterList::const_iterator clsIter = clusterList.begin(), clsIterEnd = clusterList.end(); clsIter != clsIterEnd; )
    {
        if (m_showOnlyAvailable && !(*clsIter)->IsAvailable())
        {
            clusterList.erase(clsIter++);
        }
        else
        {
            clsIter++;
        }
    }

    PANDORA_MONITORING_API(VisualizeClusters(&clusterList, listName.empty() ? "currentClusters" : listName.c_str(),
        (m_hitColors.find("particleid") != std::string::npos) ? AUTOID :
        (m_hitColors.find("particletype") != std::string::npos) ? AUTOTYPE :
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
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetCurrentPfoList(*this, pPfoList))
        {
            std::cout << "VisualMonitoringAlgorithm: current pfo list unavailable." << std::endl;
            return;
        }
    }
    else
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetPfoList(*this, listName, pPfoList))
        {
            std::cout << "VisualMonitoringAlgorithm: pfo list " << listName << " unavailable." << std::endl;
            return;
        }
    }

    PANDORA_MONITORING_API(VisualizeParticleFlowObjects(pPfoList, listName.empty() ? "currentPfos" : listName.c_str(),
        (m_hitColors.find("particleid") != std::string::npos) ? AUTOID :
        (m_hitColors.find("particletype") != std::string::npos) ? AUTOTYPE :
        (m_hitColors.find("iterate") != std::string::npos ? AUTOITER :
        (m_hitColors.find("energy") != std::string::npos ? AUTOENERGY :
        AUTO))));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VisualMonitoringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    m_showCurrentMCParticles = false;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowCurrentMCParticles", m_showCurrentMCParticles));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "MCParticleListNames", m_mcParticleListNames));

    m_showCurrentVertex = true;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowCurrentVertex", m_showCurrentVertex));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "VertexNames", m_vertexNames));

    m_showCurrentCaloHits = false;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowCurrentCaloHits", m_showCurrentCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "CaloHitListNames", m_caloHitListNames));

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

    m_hitColors = "pfo";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "HitColors", m_hitColors));
    std::transform(m_hitColors.begin(), m_hitColors.end(), m_hitColors.begin(), ::tolower);

    m_thresholdEnergy = -1.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ThresholdEnergy", m_thresholdEnergy));

    m_showOnlyAvailable = false;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowOnlyAvailable", m_showOnlyAvailable));

    m_displayEvent = true;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DisplayEvent", m_displayEvent));

    m_transparencyThresholdE = -1.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TransparencyThresholdE", m_transparencyThresholdE));

    m_energyScaleThresholdE = 1.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "EnergyScaleThresholdE", m_energyScaleThresholdE));

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

    m_blackBackground = false;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "BlackBackground", m_blackBackground));

    m_showDetector = true;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShowDetector", m_showDetector));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
