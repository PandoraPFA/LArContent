/**
 *  @file   LArContent/include/LArThreeDReco/LArEventBuilding/NeutrinoBuildingAlgorithm.h
 * 
 *  @brief  Header file for the neutrino building algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_NEUTRINO_BUILDING_ALGORITHM_H
#define LAR_NEUTRINO_BUILDING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  NeutrinoBuildingAlgorithm class
 */
class NeutrinoBuildingAlgorithm : public pandora::Algorithm
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
    void GetDaughterPfoList(pandora::PfoList &pfoList) const;

    /**
     *  @brief  Create and save the (empty) neutrino pfo, with placeholder particle id
     * 
     *  @param  pNeutrinoPfo to receive the address of the neutrino pfo
     */
    void CreateNeutrinoPfo(pandora::ParticleFlowObject *&pNeutrinoPfo) const;

    /**
     *  @brief  Add a provided list of daughters to the neutrino pfo, identifying the primary daughter and setting particle id accordingly
     * 
     *  @param  pNeutrinoPfo address of the neutrino pfo
     *  @param  daughterPfoList the list of pfos to be added as daughters of the neutrino
     */
    void AddDaughtersAndSetId(pandora::ParticleFlowObject *const pNeutrinoPfo, const pandora::PfoList &daughterPfoList) const;

    /**
     *  @brief  Get the number of two dimensional hits (TPC_VIEW_U, V or W) contained in clusters in a pfo and all its daughters
     * 
     *  @param  pPfo address of the pfo
     * 
     *  @return the number of two dimensional hits
     */
    unsigned int GetNTwoDHitsInPfo(const pandora::ParticleFlowObject *const pPfo) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector   m_daughterPfoListNames; ///< The list of daughter pfo list names
    std::string             m_outputPfoListName;    ///< The name of the output neutrino pfo list
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *NeutrinoBuildingAlgorithm::Factory::CreateAlgorithm() const
{
    return new NeutrinoBuildingAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_NEUTRINO_BUILDING_ALGORITHM_H
