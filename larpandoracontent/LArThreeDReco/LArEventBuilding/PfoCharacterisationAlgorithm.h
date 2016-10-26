/**
 *  @file   larpandoracontent/LArThreeDReco/LArEventBuilding/PfoCharacterisationAlgorithm.h
 * 
 *  @brief  Header file for the pfo characterisation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_PFO_CHARACTERISATION_ALGORITHM_H
#define LAR_PFO_CHARACTERISATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  PfoCharacterisationAlgorithm class
 */
class PfoCharacterisationAlgorithm : public pandora::Algorithm
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
     *  @brief  Whether pfo is identified as a clear track
     *
     *  @param  pPfo address of the relevant pfo
     * 
     *  @return boolean
     */
    bool IsClearTrack(const pandora::ParticleFlowObject *const pPfo) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string             m_trackPfoListName;         ///< The track pfo list name
    std::string             m_showerPfoListName;        ///< The shower pfo list name
    pandora::StringVector   m_inputPfoListNames;        ///< The names of the input pfo lists
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *PfoCharacterisationAlgorithm::Factory::CreateAlgorithm() const
{
    return new PfoCharacterisationAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_PFO_CHARACTERISATION_ALGORITHM_H
