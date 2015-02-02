/**
 *  @file   LArContent/include/LArThreeDReco/LArCosmicRay/CosmicRayVertexBuildingAlgorithm.h
 *
 *  @brief  Header file for the cosmic-ray vertex building algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_COSMIC_RAY_VERTEX_BUILDING_ALGORITHM_H
#define LAR_COSMIC_RAY_VERTEX_BUILDING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArObjects/LArPointingCluster.h"

namespace lar_content
{

/**
 *  @brief  CosmicRayVertexBuildingAlgorithm class
 */
class CosmicRayVertexBuildingAlgorithm : public pandora::Algorithm
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
    CosmicRayVertexBuildingAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Get the list of input pfos to this algorithm
     *
     *  @param  pfoList to receive the list of input pfos
     */
    void GetCosmicPfos(const pandora::PfoList *const pPfoList, pandora::PfoVector &pfoVector) const;

    /**
     *  @brief  Build a map of 3D sliding fits from the input Pfos.
     *
     *  @param  pfoList the input particle flow objects
     *  @param  pointingClusterMap the output map of 3D pointing clusters
     */
    void BuildPointingClusterMap(const pandora::PfoVector &pfoVector, LArPointingClusterMap &pointingClusterMap) const;

    /**
     *  @brief  Reconstruct the vertex and direction of a list of cosmic-ray Pfos
     *
     *  @param  pointingClusterMap the input map of 3D pointing clusters
     *  @param  pfoList the input list of Pfos
     */
    void BuildCosmicRayParticles(const LArPointingClusterMap &pointingClusterMap, const pandora::PfoVector &pfoVector) const;

    /**
     *  @brief  Reconstruct the vertex and direction of a parent cosmic-ray Pfo
     *
     *  @param  pointingClusterMap the input map of 3D pointing clusters
     *  @param  pPfo the input Pfo
     */
    void BuildCosmicRayParent(const LArPointingClusterMap &pointingClusterMap, pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Reconstruct the vertex and direction of a daughter cosmic-ray Pfo
     *
     *  @param  pPfo the daughter Pfo
     */
    void BuildCosmicRayDaughter(pandora::ParticleFlowObject *const pPfo) const;

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
    std::string             m_parentPfoListName;        ///< The name of the input pfo list
    std::string             m_vertexListName;           ///< The name of the output vertex list
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CosmicRayVertexBuildingAlgorithm::Factory::CreateAlgorithm() const
{
    return new CosmicRayVertexBuildingAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_COSMIC_RAY_VERTEX_BUILDING_ALGORITHM_H
