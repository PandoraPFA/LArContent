/**
 *  @file   LArContent/include/LArThreeDReco/LArEventBuilding/NeutrinoEventCreationAlgorithm.h
 * 
 *  @brief  Header file for the neutrino event creation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_NEUTRINO_EVENT_CREATION_ALGORITHM_H
#define LAR_NEUTRINO_EVENT_CREATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  NeutrinoEventCreationAlgorithm class
 */
class NeutrinoEventCreationAlgorithm : public pandora::Algorithm
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

    std::string        m_vertexListName;       ///< The name of the neutrino vertex list
    std::string        m_neutrinoPfoListName;  ///< The name of the neutrino pfo list
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *NeutrinoEventCreationAlgorithm::Factory::CreateAlgorithm() const
{
    return new NeutrinoEventCreationAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_NEUTRINO_EVENT_CREATION_ALGORITHM_H
