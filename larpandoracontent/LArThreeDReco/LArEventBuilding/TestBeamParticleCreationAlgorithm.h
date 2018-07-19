/**
 *  @file   larpandoracontent/LArThreeDReco/LArEventBuilding/TestBeamParticleCreationAlgorithm.h
 * 
 *  @brief  Header file for the test beam particle creation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_TEST_BEAM_PARTICLE_CLREATION_ALGORITHM_H
#define LAR_TEST_BEAM_PARTICLE_CLREATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  TestBeamParticleCreationAlgorithm class
 */
class TestBeamParticleCreationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Constructor
     */
    TestBeamParticleCreationAlgorithm();

private:
    pandora::StatusCode Run();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string    m_parentPfoListName;         ///< The parent pfo list name
    std::string    m_trackPfoListName;          ///< The track pfo list name
    std::string    m_showerPfoListName;         ///< The shower pfo list name
    std::string    m_parentVertexListName;      ///< Parent vertex list name
    std::string    m_daughterVertexListName;    ///< Daughter vertex list name
    bool           m_keepInteractionVertex;     ///< Retain the vertex for the test beam particle at the low z point
    bool           m_keepStartVertex;           ///< Retain the vertex for the test beam particle at the interaction point
};

} // namespace lar_content

#endif // #ifndef LAR_TEST_BEAM_PARTICLE_CLREATION_ALGORITHM_H
