/**
 *  @file   larpandoracontent/LArThreeDReco/LArEventBuilding/NeutrinoPropertiesAlgorithm.h
 *
 *  @brief  Header file for the neutrino properties algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_NEUTRINO_PROPERTIES_ALGORITHM_H
#define LAR_NEUTRINO_PROPERTIES_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  NeutrinoPropertiesAlgorithm class
 */
class NeutrinoPropertiesAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    NeutrinoPropertiesAlgorithm();

private:
    pandora::StatusCode Run();

    /**
     *  @brief  identifying the primary daughter of a neutrino pfo and set the particle id accordingly
     *
     *  @param  pNeutrinoPfo address of the neutrino pfo
     */
    void SetNeutrinoId(const pandora::ParticleFlowObject *const pNeutrinoPfo) const;

    /**
     *  @brief  Get the number of two dimensional hits (TPC_VIEW_U, V or W) contained in clusters in a pfo and all its daughters
     *
     *  @param  pPfo address of the pfo
     *
     *  @return the number of two dimensional hits
     */
    unsigned int GetNTwoDHitsInPfoChain(const pandora::ParticleFlowObject *const pPfo) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_neutrinoPfoListName; ///< The name of the output neutrino pfo list

    bool m_includeIsolatedHits; ///< Whether to include isolated hits when counting 2d hits in pfo chain
};

} // namespace lar_content

#endif // #ifndef LAR_NEUTRINO_PROPERTIES_ALGORITHM_H
