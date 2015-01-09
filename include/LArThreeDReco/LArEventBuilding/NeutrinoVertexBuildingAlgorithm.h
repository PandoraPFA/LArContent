/**
 *  @file   LArContent/include/LArThreeDReco/LArEventBuilding/NeutrinoVertexBuildingAlgorithm.h
 *
 *  @brief  Header file for the neutrino vertex building algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_NEUTRINO_VERTEX_BUILDING_ALGORITHM_H
#define LAR_NEUTRINO_VERTEX_BUILDING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArObjects/LArPointingCluster.h"

namespace lar_content
{

/**
 *  @brief  NeutrinoVertexBuildingAlgorithm class
 */
class NeutrinoVertexBuildingAlgorithm : public pandora::Algorithm
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
    NeutrinoVertexBuildingAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Get the vector of daughter pfos
     *
     *  @param  pfoList the input list of neutrino pfos
     *  @param  pfoVector to receive the vector of daughter pfos
     */
    void GetDaughterPfos(const pandora::PfoList *const pfoList, pandora::PfoVector &pfoVector) const;

    /**
     *  @brief  Build a map of 3D sliding fits from the input Pfos.
     *
     *  @param  pfoVector the input particle flow objects
     *  @param  pointingClusterMap the output map of 3D pointing clusters
     */
    void BuildPointingClusterMap(const pandora::PfoVector &pfoVector, LArPointingClusterMap &pointingClusterMap) const;

    /**
     *  @brief  Reconstruct the vertex and direction of daughter Pfos
     *
     *  @param  pointingClusterMap the input map of 3D pointing clusters
     *  @param  pfoVector the input list of Pfos
     */
    void BuildDaughterParticles(const LArPointingClusterMap &pointingClusterMap, const pandora::PfoVector &pfoVector) const;

    /**
     *  @brief  Reconstruct the vertex and direction of a track-like Pfos
     *
     *  @param  pointingClusterMap the input map of 3D pointing clusters
     *  @param  pPfo the input Pfo
     */
    void BuildDaughterTrack(const LArPointingClusterMap &pointingClusterMap, pandora::ParticleFlowObject *const pDaughterPfo) const;

    /**
     *  @brief  Reconstruct the vertex and direction of a shower-like Pfos
     *
     *  @param  pPfo the input Pfo
     */
    void BuildDaughterShower(pandora::ParticleFlowObject *const pDaughterPfo) const;

    /**
     *  @brief  Set the vertex and direction of the Pfos
     *
     *  @param  vtxPosition the input vertex position
     *  @param  vtxDirection the input vertex direction
     *  @param  pPfo the input Pfo
     */
    void SetParticleParameters(const pandora::CartesianVector &vtxPosition, const pandora::CartesianVector &vtxDirection,
        pandora::ParticleFlowObject *const pPfo) const;

    bool                    m_useParentShowerVertex;    ///< use the parent pfo for the shower vertices
    unsigned int            m_halfWindowLayers;         ///< number of layers to use for half-window of sliding fit
    std::string             m_neutrinoListName;         ///< The input list of pfo list names
    std::string             m_vertexListName;           ///< The name of the output cosmic-ray vertex list
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *NeutrinoVertexBuildingAlgorithm::Factory::CreateAlgorithm() const
{
    return new NeutrinoVertexBuildingAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_NEUTRINO_VERTEX_BUILDING_ALGORITHM_H
