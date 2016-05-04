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

#include "larpandoracontent/LArContent/LArObjects/LArMCParticle.h"

#include "larpandoracontent/LArContent/LArPersistency/EventReadingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

EventReadingAlgorithm::EventReadingAlgorithm() :
    m_geometryFileType(UNKNOWN_FILE_TYPE),
    m_eventFileType(UNKNOWN_FILE_TYPE),
    m_shouldReadGeometry(false),
    m_shouldReadEvents(true),
    m_skipToEvent(0),
    m_pEventFileReader(NULL)
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
    if (m_shouldReadGeometry)
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

    if (m_shouldReadEvents)
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

        m_pEventFileReader->SetFactory(new LArMCParticleFactory);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, m_pEventFileReader->GoToEvent(m_skipToEvent));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventReadingAlgorithm::Run()
{
    if ((NULL != m_pEventFileReader) && m_shouldReadEvents)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, m_pEventFileReader->ReadEvent());
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RepeatEventPreparation(*this));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventReadingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShouldReadGeometry", m_shouldReadGeometry));

    if (m_shouldReadGeometry)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
            "GeometryFileName", m_geometryFileName));

        std::string fileExtension(m_geometryFileName.substr(m_geometryFileName.find_last_of(".")));
        std::transform(fileExtension.begin(), fileExtension.end(), fileExtension.begin(), ::tolower);

        if (std::string(".xml") == fileExtension)
        {
            m_geometryFileType = XML;
        }
        else if (std::string(".pndr") == fileExtension)
        {
            m_geometryFileType = BINARY;
        }
        else
        {
            std::cout << "EventReadingAlgorithm: Unknown geometry file type specified " << std::endl;
            return STATUS_CODE_INVALID_PARAMETER;
        }
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShouldReadEvents", m_shouldReadEvents));

    if (m_shouldReadEvents)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
            "EventFileName", m_eventFileName));

        std::string fileExtension(m_eventFileName.substr(m_eventFileName.find_last_of(".")));
        std::transform(fileExtension.begin(), fileExtension.end(), fileExtension.begin(), ::tolower);

        if (std::string(".xml") == fileExtension)
        {
            m_eventFileType = XML;
        }
        else if (std::string(".pndr") == fileExtension)
        {
            m_eventFileType = BINARY;
        }
        else
        {
            std::cout << "EventReadingAlgorithm: Unknown event file type specified " << std::endl;
            return STATUS_CODE_INVALID_PARAMETER;
        }
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SkipToEvent", m_skipToEvent));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
