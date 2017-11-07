/**
 *  @file   larpandoracontent/LArTrackShowerId/SvmPfoCharacterisationAlgorithm.cc
 * 
 *  @brief  Implementation of the svm pfo characterisation algorithm class.
 * 
 *  $Log: $
 */
#ifdef CETLIB_AVAILABLE
#include "cetlib/search_path.h"
#endif

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArSvmHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArTrackShowerId/SvmPfoCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

SvmPfoCharacterisationAlgorithm::SvmPfoCharacterisationAlgorithm() :
    m_trainingSetMode(false),
    m_enableProbability(false),
    m_useThreeDInformation(false),
    m_minProbabilityCut(0.5f),
    m_minCaloHitsCut(5)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SvmPfoCharacterisationAlgorithm::IsClearTrack(const Cluster *const pCluster) const
{
    if (pCluster->GetNCaloHits() < m_minCaloHitsCut)
        return false;

    SupportVectorMachine::DoubleVector featureVector(LArSvmHelper::CalculateFeatures(m_featureToolVector, this, pCluster));

    if (m_trainingSetMode)
    {
        bool isTrueTrack(false);

        try
        {
            const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCluster));
            isTrueTrack = ((PHOTON != pMCParticle->GetParticleId()) && (E_MINUS != std::abs(pMCParticle->GetParticleId())));
        }
        catch (const StatusCodeException &) {}

        LArSvmHelper::ProduceTrainingExample(m_trainingOutputFile, isTrueTrack, featureVector);
        return isTrueTrack;
    }

    if (!m_enableProbability)
    {
        return LArSvmHelper::Classify(m_supportVectorMachine, featureVector);
    }
    else
    {
        return (LArSvmHelper::CalculateProbability(m_supportVectorMachine, featureVector) > m_minProbabilityCut);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SvmPfoCharacterisationAlgorithm::IsClearTrack(const pandora::ParticleFlowObject *const pPfo) const
{
    //charge related features are only calculated using hits in W view
    ClusterList wClusterList;
    LArPfoHelper::GetClusters(pPfo, TPC_VIEW_W, wClusterList);

    PfoCharacterisationFeatureTool::FeatureToolVector featureToolVector(wClusterList.empty() ? m_featureToolVectorMissingW : m_featureToolVectorThreeD);
    SupportVectorMachine::DoubleVector featureVector(LArSvmHelper::CalculateFeatures(featureToolVector, this, pPfo));
    SupportVectorMachine supportVectorMachine(wClusterList.empty() ? m_supportVectorMachineMissingW : m_supportVectorMachine);

    if (m_trainingSetMode)
    {
        bool isTrueTrack(false);
        int trueMCPdg(0);

        try
        {
            const MCParticle *const pMCParticle(LArMCParticleHelper::GetMainMCParticle(pPfo));
            isTrueTrack = ((PHOTON != pMCParticle->GetParticleId()) && (E_MINUS != std::abs(pMCParticle->GetParticleId())));
            trueMCPdg = pMCParticle->GetParticleId();
        }
        catch (const StatusCodeException &) {}

        if (trueMCPdg!=0)
        {
	  std::string outputFile;
            outputFile.append(m_trainingOutputFile);
            const char* end=((wClusterList.empty()) ? "MissingW.txt" : ".txt");
            outputFile.append(end);
            LArSvmHelper::ProduceTrainingExample(outputFile, isTrueTrack, featureVector);
        }
        return isTrueTrack;
    }// training mode

    if (!m_enableProbability)
    {
        return LArSvmHelper::Classify(supportVectorMachine, featureVector);
    }
    else
    {
        return (m_minProbabilityCut <= LArSvmHelper::CalculateProbability(supportVectorMachine, featureVector));
    }
 
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SvmPfoCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TrainingSetMode", m_trainingSetMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCaloHitsCut", m_minCaloHitsCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseThreeDInformation", m_useThreeDInformation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SvmFileName", m_svmFileName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SvmName", m_svmName));

    if (m_useThreeDInformation)
    {
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
            "SvmFileNameMissingW", m_svmFileNameMissingW));

        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
            "SvmNameMissingW", m_svmNameMissingW));
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "EnableProbability",  m_enableProbability));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinProbabilityCut", m_minProbabilityCut));

    if (m_trainingSetMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrainingOutputFileName", m_trainingOutputFile));
    }
    else
    {
        if (m_svmFileName.empty() || m_svmName.empty())
        {
            std::cout << "SvmPfoCharacterisationAlgorithm: SvmFileName and SvmName must be set if in classification mode " << std::endl;
            return STATUS_CODE_INVALID_PARAMETER;
        }

        std::string fullSvmFileName(m_svmFileName), fullSvmFileNameMissingW(m_svmFileNameMissingW);;
#ifdef CETLIB_AVAILABLE
        cet::search_path sp("FW_SEARCH_PATH");

        if (!sp.find_file(m_svmFileName, fullSvmFileName))
        {
            std::cout << " SvmPfoCharacterisationAlgorithm::ReadSettings - Failed to find svm file " << m_svmFileName << " in FW search path" << std::endl;
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
        }
        if (m_useThreeDInformation && (!sp.find_file(m_svmFileNameMissingW, fullSvmFileNameMissingW)))
        {
            std::cout << " SvmPfoCharacterisationAlgorithm::ReadSettings - Failed to find svm file " << m_svmFileNameMissingW << " in FW search path" << std::endl;
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
        }
#endif
        m_supportVectorMachine.Initialize(fullSvmFileName, m_svmName);
        if (m_useThreeDInformation)
            m_supportVectorMachineMissingW.Initialize(fullSvmFileNameMissingW, m_svmNameMissingW);
    }

    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "FeatureTools", algorithmToolVector));

    if (m_useThreeDInformation)
    {
        AlgorithmToolVector algorithmToolVectorMissingW;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "FeatureToolsMissingW", algorithmToolVectorMissingW));
        for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArSvmHelper::AddFeatureToolToVector(pAlgorithmTool, m_featureToolVectorThreeD));
        for (AlgorithmTool *const pAlgorithmTool : algorithmToolVectorMissingW)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArSvmHelper::AddFeatureToolToVector(pAlgorithmTool, m_featureToolVectorMissingW));
    }
    else
    {
        for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArSvmHelper::AddFeatureToolToVector(pAlgorithmTool, m_featureToolVector));x
    }

    return PfoCharacterisationBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
