/**
 *  @file   larpandoracontent/LArCheating/CheatingPfoCharacterisationAlgorithm.h
 * 
 *  @brief  Header file for the cheating pfo characterisation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_CHEATING_PFO_CHARACTERISATION_ALGORITHM_H
#define LAR_CHEATING_PFO_CHARACTERISATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  CheatingPfoCharacterisationAlgorithm class
 */
class CheatingPfoCharacterisationAlgorithm : public pandora::Algorithm
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
    CheatingPfoCharacterisationAlgorithm();

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

    bool                    m_updateClusterIds;         ///< Whether to update daughter cluster particle id labels to match pfo id.
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CheatingPfoCharacterisationAlgorithm::Factory::CreateAlgorithm() const
{
    return new CheatingPfoCharacterisationAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_PFO_CHARACTERISATION_ALGORITHM_H
