/**
 *  @file   LArContent/src/LArReclustering/LowEnergyTrackSplittingAlgorithm.cc
 * 
 *  @brief  Implementation of the low energy track splitting algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArReclustering/LowEnergyTrackSplittingAlgorithm.h"

using namespace pandora;

namespace lar
{

LowEnergyTrackSplittingAlgorithm::LowEnergyTrackSplittingAlgorithm()
{
    PANDORA_MONITORING_API(Create());
}

//------------------------------------------------------------------------------------------------------------------------------------------

LowEnergyTrackSplittingAlgorithm::~LowEnergyTrackSplittingAlgorithm()
{
    PANDORA_MONITORING_API(SaveTree(m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LowEnergyTrackSplittingAlgorithm::Run()
{
    const ClusterList *pInputClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentClusterList(*this, m_inputSeedClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentClusterList(*this, pInputClusterList));

    ClusterVector clusterVector(pInputClusterList->begin(), pInputClusterList->end());
    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNOccupiedLayers);

    // Monitoring containers
    FloatVector trackLengths, trackWidths;
    IntVector pdgCodes;

    // The target lists
    ClusterList allShowerSeedClusters, allTrackSeedClusters, allNonSeedClusters;

    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        // Store monitoring information
        int pdgCode(-std::numeric_limits<int>::max());

        try
        {
            const MCParticle *pMCParticle(MCParticleHelper::GetMainMCParticle(pCluster));
            pdgCode = pMCParticle->GetParticleId();
        }
        catch (StatusCodeException &)
        {
        }

        float trackLength(-std::numeric_limits<float>::max());

        try
        {
            trackLength = LArClusterHelper::GetLength(pCluster);
        }
        catch (StatusCodeException &)
        {
        }

        float trackWidth(-std::numeric_limits<float>::max());

        try
        {
            trackWidth = LArClusterHelper::LArTrackWidth(pCluster);
        }
        catch (StatusCodeException &)
        {
        }

        pdgCodes.push_back(pdgCode);
        trackLengths.push_back(trackLength);
        trackWidths.push_back(trackWidth);

        if ((std::abs(pdgCode) == PROTON) || (std::abs(pdgCode) == MU_MINUS) || (std::abs(pdgCode) == PI_PLUS))
        {
            allTrackSeedClusters.insert(pCluster);
        }
        else if ((std::abs(pdgCode) == E_MINUS) || (std::abs(pdgCode) == PHOTON))
        {
            allShowerSeedClusters.insert(pCluster);
        }
        else
        {
            allNonSeedClusters.insert(pCluster);
        }
    }

//std::cout << "Cheated Track/Shower Separation " << std::endl;
//PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f); 
//PandoraMonitoringApi::VisualizeClusters(&allTrackSeedClusters, "allTrackSeedClusters", RED);
//PandoraMonitoringApi::VisualizeClusters(&allShowerSeedClusters, "allShowerSeedClusters", GREEN);
//PandoraMonitoringApi::VisualizeClusters(&allNonSeedClusters, "allNonSeedClusters", BLUE);
//PandoraMonitoringApi::ViewEvent();

    // Cluster output
    if (!allShowerSeedClusters.empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveClusterList(*this, m_inputSeedClusterListName,
            m_showerSeedClusterListName, allShowerSeedClusters));
    }

    if (!allTrackSeedClusters.empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveClusterList(*this, m_inputSeedClusterListName,
            m_trackSeedClusterListName, allTrackSeedClusters));
    }

    if (!allNonSeedClusters.empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveClusterList(*this, m_inputSeedClusterListName,
            m_nonSeedClusterListName, allNonSeedClusters));
    }

    if (!pInputClusterList->empty())
        return STATUS_CODE_FAILURE;

    // Write monitoring information
    PANDORA_MONITORING_API(SetTreeVariable(m_treeName.c_str(), "pdgCodes", &pdgCodes));
    PANDORA_MONITORING_API(SetTreeVariable(m_treeName.c_str(), "trackLengths", &trackLengths));
    PANDORA_MONITORING_API(SetTreeVariable(m_treeName.c_str(), "trackWidths", &trackWidths));
    PANDORA_MONITORING_API(FillTree(m_treeName.c_str()));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LowEnergyTrackSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputSeedClusterListName",
        m_inputSeedClusterListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ShowerSeedClusterListName",
        m_showerSeedClusterListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrackSeedClusterListName",
        m_trackSeedClusterListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NonSeedClusterListName",
        m_nonSeedClusterListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputTree", m_treeName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputFile", m_fileName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
