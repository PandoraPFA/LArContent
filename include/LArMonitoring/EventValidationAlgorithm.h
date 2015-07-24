/**
 *  @file   LArContent/include/LArMonitoring/EventValidationAlgorithm.h
 *
 *  @brief  Header file for the event validation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_EVENT_VALIDATION_ALGORITHM_H
#define LAR_EVENT_VALIDATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  EventValidationAlgorithm class
 */
class EventValidationAlgorithm: public pandora::Algorithm
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
    EventValidationAlgorithm();

    /**
     *  @brief  Destructor
     */
    ~EventValidationAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string     m_caloHitListName;          ///< Name of input calo hit list
    std::string     m_mcParticleListName;       ///< Name of input MC particle list
    std::string     m_pfoListName;              ///< Name of input Pfo list
    std::string     m_fileName;                 ///< Name of output file
    std::string     m_treeName;                 ///< Name of output tree

    bool            m_extractNeutrinoDaughters; ///< Whether to treat each neutrino pfo daughter as a standalone top-level pfo
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *EventValidationAlgorithm::Factory::CreateAlgorithm() const
{
    return new EventValidationAlgorithm();
}

} // namespace lar_content

#endif // LAR_EVENT_VALIDATION_ALGORITHM_H
