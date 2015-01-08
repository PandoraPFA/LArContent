/**
 *  @file   LArContent/include/LArThreeDReco/LArEventBuilding/NeutrinoEventBuildingAlgorithm.h
 *
 *  @brief  Header file for the neutrino event building algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_NEUTRINO_EVENT_BUILDING_ALGORITHM_H
#define LAR_NEUTRINO_EVENT_BUILDING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  NeutrinoEventBuildingAlgorithm class
 */
class NeutrinoEventBuildingAlgorithm : public pandora::Algorithm
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

    /**
     *  @brief  Get the list of daughter pfos to be added to the new neutrino pfo
     *
     *  @param  pfoList to receive the list of daughter pfos
     */
    void GetNeutrinoPfo(pandora::ParticleFlowObject *&pNeutrinoPfo) const;

    /**
     *  @brief  Get the list of daughter pfos to be added to the new neutrino pfo
     *
     *  @param  pfoList to receive the list of daughter pfos
     */
    void GetDaughterPfoList(pandora::PfoList &pfoList) const;

    /**
     *  @brief  Add a provided list of daughters to the neutrino pfo
     *
     *  @param  pNeutrinoPfo address of the neutrino pfo
     *  @param  daughterPfoList the list of pfos to be added as daughters of the neutrino
     */
    void AddDaughters(pandora::ParticleFlowObject *const pNeutrinoPfo, const pandora::PfoList &daughterPfoList) const;

    /**
     *  @brief  identifying the primary daughter of a neutrino pfo and set the particle id accordingly
     *
     *  @param  pNeutrinoPfo address of the neutrino pfo
     */
    void SetNeutrinoId(pandora::ParticleFlowObject *const pNeutrinoPfo) const;

    /**
     *  @brief  Get the number of two dimensional hits (TPC_VIEW_U, V or W) contained in clusters in a pfo and all its daughters
     *
     *  @param  pPfo address of the pfo
     *
     *  @return the number of two dimensional hits
     */
    unsigned int GetNTwoDHitsInPfo(const pandora::ParticleFlowObject *const pPfo) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector   m_daughterPfoListNames;   ///< The list of daughter pfo list names
    std::string             m_neutrinoPfoListName;    ///< The name of the output neutrino pfo list
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *NeutrinoEventBuildingAlgorithm::Factory::CreateAlgorithm() const
{
    return new NeutrinoEventBuildingAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_NEUTRINO_EVENT_BUILDING_ALGORITHM_H
