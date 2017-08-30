/**
 *  @file   PandoraSDK/src/Persistency/EventReadingAlgorithm.cc
 * 
 *  @brief  Implementation of the event reading algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "Persistency/BinaryFileReader.h"
#include "Persistency/XmlFileReader.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"

#include "larpandoracontent/LArPersistency/EventReadingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

EventReadingAlgorithm::EventReadingAlgorithm() :
    m_geometryFileType(UNKNOWN_FILE_TYPE),
    m_eventFileType(UNKNOWN_FILE_TYPE),
    m_skipToEvent(0),
    m_useLArCaloHits(false),
    m_useLArMCParticles(true),
    m_pEventFileReader(nullptr)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

EventReadingAlgorithm::~EventReadingAlgorithm()
{
    delete m_pEventFileReader;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventReadingAlgorithm::Initialize()
{
    if (!m_geometryFileName.empty())
    {
        if (BINARY == m_geometryFileType)
        {
            BinaryFileReader fileReader(this->GetPandora(), m_geometryFileName);
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, fileReader.ReadGeometry());
        }
        else if (XML == m_geometryFileType)
        {
            XmlFileReader fileReader(this->GetPandora(), m_geometryFileName);
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, fileReader.ReadGeometry());
        }
        else
        {
            return STATUS_CODE_FAILURE;
        }
    }

    if (!m_eventFileName.empty())
    {
        if (BINARY == m_eventFileType)
        {
            m_pEventFileReader = new BinaryFileReader(this->GetPandora(), m_eventFileName);
        }
        else if (XML == m_eventFileType)
        {
            m_pEventFileReader = new XmlFileReader(this->GetPandora(), m_eventFileName);
        }
        else
        {
            return STATUS_CODE_FAILURE;
        }

        if (m_useLArCaloHits)
            m_pEventFileReader->SetFactory(new LArCaloHitFactory);

        if (m_useLArMCParticles)
            m_pEventFileReader->SetFactory(new LArMCParticleFactory);

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, m_pEventFileReader->GoToEvent(m_skipToEvent));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventReadingAlgorithm::Run()
{
    if ((nullptr != m_pEventFileReader) && !m_eventFileName.empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, m_pEventFileReader->ReadEvent());
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RepeatEventPreparation(*this));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

FileType EventReadingAlgorithm::GetFileType(const std::string &fileName) const
{
    std::string fileExtension(fileName.substr(fileName.find_last_of(".")));
    std::transform(fileExtension.begin(), fileExtension.end(), fileExtension.begin(), ::tolower);

    if (std::string(".xml") == fileExtension)
    {
        return XML;
    }
    else if (std::string(".pndr") == fileExtension)
    {
        return BINARY;
    }
    else
    {
        std::cout << "EventReadingAlgorithm: Unknown file type specified " << fileName << std::endl;
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventReadingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    ExternalEventReadingParameters *pExternalParameters(nullptr);

    if (this->ExternalParametersPresent())
    {
        this->RegisterParameterAccessAttempt();
        pExternalParameters = dynamic_cast<ExternalEventReadingParameters*>(this->GetExternalParameters());

        if (!pExternalParameters)
            return STATUS_CODE_FAILURE;
    }

    if (pExternalParameters && !pExternalParameters->m_geometryFileName.empty())
    {
        m_geometryFileName = pExternalParameters->m_geometryFileName;
    }
    else
    {
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "GeometryFileName", m_geometryFileName));
    }

    if (pExternalParameters && !pExternalParameters->m_eventFileName.empty())
    {
        m_eventFileName = pExternalParameters->m_eventFileName;
    }
    else
    {
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "EventFileName", m_eventFileName));
    }

    if (pExternalParameters && pExternalParameters->m_skipToEvent.IsInitialized())
    {
        m_skipToEvent = pExternalParameters->m_skipToEvent.Get();
    }
    else
    {
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SkipToEvent", m_skipToEvent));
    }

    if (!m_geometryFileName.empty())
        m_geometryFileType = this->GetFileType(m_geometryFileName);

    if (!m_eventFileName.empty())
        m_eventFileType = this->GetFileType(m_eventFileName);

    if (m_geometryFileName.empty() && m_eventFileName.empty())
    {
        std::cout << "EventReadingAlgorithm - nothing to do; neither geometry nor event file specified." << std::endl;
        return STATUS_CODE_NOT_INITIALIZED;
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseLArCaloHits", m_useLArCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseLArMCParticles", m_useLArMCParticles));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
