/**
 *  @file   LArContent/include/LArPersistency/EventWritingAlgorithm.h
 * 
 *  @brief  Header file for the event writing algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_EVENT_WRITING_ALGORITHM_H
#define LAR_EVENT_WRITING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "Persistency/PandoraIO.h"

namespace pandora {class FileWriter;}

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_content
{

/**
 *  @brief  EventWritingAlgorithm class
 */
class EventWritingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

    /**
     *  @brief  Default constructor
     */
    EventWritingAlgorithm();

    /**
     *  @brief  Destructor
     */
    ~EventWritingAlgorithm();

private:
    pandora::StatusCode Initialize();
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::FileType       m_geometryFileType;             ///< The geometry file type
    pandora::FileType       m_eventFileType;                ///< The event file type

    bool                    m_shouldWriteGeometry;          ///< Whether to write geometry to a specified file
    std::string             m_geometryFileName;             ///< Name of the output geometry file

    bool                    m_shouldWriteEvents;            ///< Whether to write events to a specified file
    std::string             m_eventFileName;                ///< Name of the output event file

    bool                    m_shouldWriteMCRelationships;   ///< Whether to write mc relationship information to the events file
    bool                    m_shouldWriteTrackRelationships;///< Whether to write track relationship information to the events file

    bool                    m_shouldOverwriteEventFile;     ///< Whether to overwrite existing event file with specified name, or append
    bool                    m_shouldOverwriteGeometryFile;  ///< Whether to overwrite existing geometry file with specified name, or append

    pandora::FileWriter    *m_pEventFileWriter;             ///< Address of the event file writer
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *EventWritingAlgorithm::Factory::CreateAlgorithm() const
{
    return new EventWritingAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_EVENT_WRITING_ALGORITHM_H
