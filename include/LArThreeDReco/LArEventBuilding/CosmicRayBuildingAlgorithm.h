/**
 *  @file   LArContent/include/LArThreeDReco/LArEventBuilding/CosmicRayBuildingAlgorithm.h
 *
 *  @brief  Header file for the cosmic-ray building algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_COSMIC_RAY_BUILDING_ALGORITHM_H
#define LAR_COSMIC_RAY_BUILDING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArObjects/LArPointingCluster.h"

namespace lar_content
{

/**
 *  @brief  CosmicRayBuildingAlgorithm class
 */
class CosmicRayBuildingAlgorithm : public pandora::Algorithm
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
    CosmicRayBuildingAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Get the list of input pfos to this algorithm
     *
     *  @param  pfoList to receive the list of input pfos
     */
    void GetInputPfoList(pandora::PfoList &pfoList) const;

    /**
     *  @brief  Build a map of 3D sliding fits from the input Pfos.
     *
     *  @param  pfoList the input particle flow objects
     *  @param  pointingClusterMap the output map of 3D pointing clusters
     */
    void BuildPointingClusterMap(const pandora::PfoList &pfoList, LArPointingClusterMap &pointingClusterMap) const;

    /**
     *  @brief  Reconstruct the vertex and direction of a list of cosmic-ray Pfos
     *
     *  @param  pointingClusterMap the input map of 3D pointing clusters
     *  @param  pfoList the input list of Pfos
     */
    void BuildCosmicRayParticles(const LArPointingClusterMap &pointingClusterMap, const pandora::PfoList &pfoList) const;

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

    unsigned int            m_halfWindowLayers;    ///< number of layers to use for half-window of sliding fit
    pandora::StringVector   m_pfoListNames;        ///< The input list of pfo list names
    std::string             m_vertexListName;      ///< The name of the output cosmic-ray vertex list
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CosmicRayBuildingAlgorithm::Factory::CreateAlgorithm() const
{
    return new CosmicRayBuildingAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_COSMIC_RAY_BUILDING_ALGORITHM_H
