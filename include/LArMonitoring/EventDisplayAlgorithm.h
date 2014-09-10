/**
 *  @file   LArContent/include/LArMonitoring/EventDisplayAlgorithm.h
 *
 *  @brief  Header file for the event display algorithm
 *
 *  $Log: $
 */

#ifndef LAR_EVENT_DISPLAY_ALGORITHM_H
#define LAR_EVENT_DISPLAY_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#ifdef MONITORING
#include "PandoraMonitoringApi.h"
#endif

namespace lar_content
{

/**
 *  @brief  EventDisplayAlgorithm class
 */
class EventDisplayAlgorithm : public pandora::Algorithm
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

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

#ifdef MONITORING
    Color GetColor( unsigned int icolor );
#endif

    std::string        m_clusterListName;
    std::string        m_particleListName;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *EventDisplayAlgorithm::Factory::CreateAlgorithm() const
{
    return new EventDisplayAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_EVENT_DISPLAY_ALGORITHM_H
