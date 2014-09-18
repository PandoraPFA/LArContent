/**
 *  @file   LArContent/include/LArUtility/ListPreparationAlgorithm.h
 * 
 *  @brief  Header file for the list preparation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_LIST_PREPARATION_ALGORITHM_H
#define LAR_LIST_PREPARATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  ListPreparationAlgorithm class
 */
class ListPreparationAlgorithm : public pandora::Algorithm
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
    ListPreparationAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief Build separate CaloHitLists for each view
     */
    void ProcessCaloHits();

    /**
     *  @brief Clean up the input CaloHitList
     *
     *  @param inputList the input CaloHitList
     *  @param outputList the output CaloHitList
     */
    void GetFilteredCaloHitList(const pandora::CaloHitList &inputList, pandora::CaloHitList &outputList);

    /**
     *  @brief Build separate MCParticleLists for each view
     */
    void ProcessMCParticles();

    float           m_mipEquivalentCut;                 ///< Minimum mip equivalent energy for calo hit
    bool            m_onlyAvailableCaloHits;            ///< Whether to only include available calo hits
    std::string     m_inputCaloHitListName;             ///< The input calo hit list name
    std::string     m_outputCaloHitListNameU;           ///< The output calo hit list name for TPC_VIEW_U hits
    std::string     m_outputCaloHitListNameV;           ///< The output calo hit list name for TPC_VIEW_V hits
    std::string     m_outputCaloHitListNameW;           ///< The output calo hit list name for TPC_VIEW_W hits
    std::string     m_filteredCaloHitListName;          ///< The output calo hit list name for all U, V and W hits
    std::string     m_currentCaloHitListReplacement;    ///< The name of the calo hit list to replace the current list (optional)

    bool            m_mcNeutrinoSelection;              ///< Whether to only select MC particles associated with a neutrino
    std::string     m_inputMCParticleListName;          ///< The input MC particle list name
    std::string     m_outputMCParticleListNameU;        ///< The output MC particle list name for MC_VIEW_U particles
    std::string     m_outputMCParticleListNameV;        ///< The output MC particle list name for MC_VIEW_V particles
    std::string     m_outputMCParticleListNameW;        ///< The output MC particle list name for MC_VIEW_W particles
    std::string     m_outputMCParticleListName3D;       ///< The output MC particle list name for 3D particles
    std::string     m_currentMCParticleListReplacement; ///< The name of the MC particle list to replace the current list (optional)
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ListPreparationAlgorithm::Factory::CreateAlgorithm() const
{
    return new ListPreparationAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_LIST_PREPARATION_ALGORITHM_H
