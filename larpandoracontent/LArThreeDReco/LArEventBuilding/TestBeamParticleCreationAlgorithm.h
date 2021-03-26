/**
 *  @file   larpandoracontent/LArThreeDReco/LArEventBuilding/TestBeamParticleCreationAlgorithm.h
 *
 *  @brief  Header file for the test beam particle creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_TEST_BEAM_PARTICLE_CREATION_ALGORITHM_H
#define LAR_TEST_BEAM_PARTICLE_CREATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  TestBeamParticleCreationAlgorithm class
 */
class TestBeamParticleCreationAlgorithm : public pandora::Algorithm
{
private:
    pandora::StatusCode Run();

    /**
     *  @brief  Set up the test beam pfo
     *
     *  @param  pNuPfo the input neutrino-hypothesis pfo
     *  @param  pTestBeamPfo to receive the output test-beam-hypothesis pfo
     *  @param  testBeamStartVertex to receive the position of the test beam start-position vertex (position of hit at minimum z)
     *
     *  @return status code
     */
    pandora::StatusCode SetupTestBeamPfo(
        const pandora::Pfo *const pNuPfo, const pandora::Pfo *&pTestBeamPfo, pandora::CartesianVector &testBeamStartVertex) const;

    /**
     *  @brief  Set up the test beam vertex
     *
     *  @param  pNuPfo the input neutrino-hypothesis pfo
     *  @param  pTestBeamPfo the input test-beam-hypothesis pfo
     *  @param  testBeamStartVertex the input position of the test beam start-position vertex
     *
     *  @return status code
     */
    pandora::StatusCode SetupTestBeamVertex(
        const pandora::Pfo *const pNuPfo, const pandora::Pfo *const pTestBeamPfo, const pandora::CartesianVector &testBeamStartVertex) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_parentPfoListName; ///< The parent pfo list name
    std::string m_trackPfoListName;  ///< The track pfo list name
    std::string m_showerPfoListName; ///< The shower pfo list name

    std::string m_parentVertexListName;   ///< The parent vertex list name
    std::string m_daughterVertexListName; ///< The daughter vertex list name
};

} // namespace lar_content

#endif // #ifndef LAR_TEST_BEAM_PARTICLE_CREATION_ALGORITHM_H
