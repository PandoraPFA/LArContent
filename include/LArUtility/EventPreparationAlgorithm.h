/**
 *  @file   LArContent/include/LArUtility/EventPreparationAlgorithm.h
 * 
 *  @brief  Header file for the event preparation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_EVENT_PREPARATION_ALGORITHM_H
#define LAR_EVENT_PREPARATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  EventPreparationAlgorithm class
 */
class EventPreparationAlgorithm : public pandora::Algorithm
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

    /**
     *  @brief  Obtain a final version of a specified hit list; if there are multiple hits in the same location, take only most energetic
     * 
     *  @param  originalList the original hit list
     *  @param  finalList to receive the final hit list
     */
    pandora::StatusCode FinalizeHitList(const pandora::CaloHitList &originalList, pandora::CaloHitList &finalList) const;

    float           m_mipEquivalentCut;                 ///< Minimum mip equivalent energy for calo hit
    float           m_minCaloHitSeparationSquared;      ///< Square of minimum calo hit separation

    std::string     m_outputCaloHitListNameU;           ///< The output calo hit list name for VIEW_U hits
    std::string     m_outputCaloHitListNameV;           ///< The output calo hit list name for VIEW_V hits
    std::string     m_outputCaloHitListNameW;           ///< The output calo hit list name for VIEW_W hits
    std::string     m_filteredCaloHitListName;          ///< The output calo hit list name for all U, V and W hits
    std::string     m_currentCaloHitListReplacement;    ///< The name of the calo hit list to replace the current list (optional)

    std::string     m_outputMCParticleListNameU;        ///< The output MC particle list name for MC_VIEW_U hits
    std::string     m_outputMCParticleListNameV;        ///< The output MC particle list name for MC_VIEW_V hits
    std::string     m_outputMCParticleListNameW;        ///< The output MC particle list name for MC_VIEW_W hits
    std::string     m_outputMCParticleListName3D;       ///< The output MC particle list name for MC_VIEW_U hits
    std::string     m_currentMCParticleListReplacement; ///< The name of the MC particle list to replace the current list (optional)
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *EventPreparationAlgorithm::Factory::CreateAlgorithm() const
{
    return new EventPreparationAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_EVENT_PREPARATION_ALGORITHM_H
