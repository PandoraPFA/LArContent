/**
 *  @file   larpandoracontent/LArTrackShowerId/SvmPfoCharacterisationAlgorithm.cc
 * 
 *  @brief  Implementation of the svm pfo characterisation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArSvmHelper.h"

#include "larpandoracontent/LArTrackShowerId/SvmPfoCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

SvmPfoCharacterisationAlgorithm::SvmPfoCharacterisationAlgorithm() :
    m_trainingSetMode(false),
    m_ratioVariables(true),
    m_enableProbability(false),
    m_minProbabilityCut(0.5f),
    m_minCaloHitsCut(5),
    m_filePathEnvironmentVariable("FW_SEARCH_PATH")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SvmPfoCharacterisationAlgorithm::IsClearTrack(const Cluster *const pCluster) const
{
    if (pCluster->GetNCaloHits() < m_minCaloHitsCut)
        return false;

    SupportVectorMachine::DoubleVector featureVector(LArSvmHelper::CalculateFeatures(m_featureToolVector, this, pCluster));

    if (m_ratioVariables)
    {
        // TODO This assumption is very bad - remove
        const double straightLineLength(featureVector.at(0));

        if (straightLineLength > std::numeric_limits<double>::epsilon())
        {
            for (unsigned int i = 1; i < featureVector.size(); ++i)
                featureVector[i] /= straightLineLength;
        }
    }

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

StatusCode SvmPfoCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TrainingSetMode", m_trainingSetMode));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "RatioVariables", m_ratioVariables));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCaloHitsCut", m_minCaloHitsCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FilePathEnvironmentVariable", m_filePathEnvironmentVariable));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SvmFileName", m_svmFileName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SvmName", m_svmName));

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

        const std::string fullSvmFileName(LArFileHelper::FindFileInPath(m_svmFileName, m_filePathEnvironmentVariable));
        m_supportVectorMachine.Initialize(fullSvmFileName, m_svmName);
    }

    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "FeatureTools", algorithmToolVector));

    for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArSvmHelper::AddFeatureToolToVector(pAlgorithmTool, m_featureToolVector));

    return PfoCharacterisationBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
