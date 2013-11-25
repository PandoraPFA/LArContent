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

#include "PandoraMonitoringApi.h"

namespace lar
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

    std::string        m_seedClusterListName;
    std::string        m_nonSeedClusterListName;
    std::string        m_vertexName;
    std::string        m_particleListName;

};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *EventDisplayAlgorithm::Factory::CreateAlgorithm() const
{
    return new EventDisplayAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_EVENT_DISPLAY_ALGORITHM_H
