/**
 *  @file   larpandoracontent/LArTrackShowerId/PfoCharacterisationAlgorithm.h
 * 
 *  @brief  Header file for the pfo characterisation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_PFO_CHARACTERISATION_ALGORITHM_H
#define LAR_PFO_CHARACTERISATION_ALGORITHM_H 1

#include "larpandoracontent/LArTrackShowerId/ShowerGrowingAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  PfoCharacterisationAlgorithm class
 */
class PfoCharacterisationAlgorithm : public ShowerGrowingAlgorithm
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
    PfoCharacterisationAlgorithm();

    /**
     *  @brief  Destructor
     */
    ~PfoCharacterisationAlgorithm();

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

    std::string             m_clusterListNameU;         ///< The u cluster list name
    std::string             m_clusterListNameV;         ///< The v cluster list name
    std::string             m_clusterListNameW;         ///< The w cluster list name

    bool                    m_updateClusterIds;         ///< Whether to update daughter cluster particle id labels to match pfo id
    bool                    m_postBranchAddition;       ///< Whether to use configuration for shower clusters post branch addition

    bool                    m_writeToTree;              ///< Whether to write monitoring details to tree
    std::string             m_treeName;                 ///< Name of output tree
    std::string             m_fileName;                 ///< Name of output file
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *PfoCharacterisationAlgorithm::Factory::CreateAlgorithm() const
{
    return new PfoCharacterisationAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_PFO_CHARACTERISATION_ALGORITHM_H
