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

namespace lar
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

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief Build separate MCParticleLists for each view
     */
    pandora::StatusCode ProcessMCParticles();

    /**
     *  @brief Build separate CaloHitLists for each view
     */
    pandora::StatusCode ProcessCaloHits();


    float           m_mipEquivalentCut;                 ///< Minimum mip equivalent energy for calo hit

    std::string     m_outputCaloHitListNameU;           ///< The output calo hit list name for VIEW_U hits
    std::string     m_outputCaloHitListNameV;           ///< The output calo hit list name for VIEW_V hits
    std::string     m_outputCaloHitListNameW;           ///< The output calo hit list name for VIEW_W hits
    std::string     m_filteredCaloHitListName;          ///< The output calo hit list name for all U, V and W hits
    std::string     m_currentCaloHitListReplacement;    ///< The name of the calo hit list to replace the current list (optional)

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

} // namespace lar

#endif // #ifndef LAR_LIST_PREPARATION_ALGORITHM_H
