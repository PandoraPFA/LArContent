/**
 *  @file   larpandoracontent/LArHelpers/LArMonitoringHelper.h
 *
 *  @brief  Header file for the lar monitoring helper helper class.
 *
 *  $Log: $
 */
#ifndef LAR_MONITORING_HELPER_H
#define LAR_MONITORING_HELPER_H 1

#include "Pandora/PandoraInternal.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content
{

/**
 *  @brief  LArMonitoringHelper class
 */
class LArMonitoringHelper
{
public:
    /**
     *  @brief  Count the number of calo hits, in a provided list, of a specified type
     *
     *  @param  hitType the hit type
     *  @param  caloHitList the calo hit list
     *
     *  @return the number of calo hits of the specified type
     */
    static unsigned int CountHitsByType(const pandora::HitType hitType, const pandora::CaloHitList &caloHitList);

    /**
     *  @brief  Order input MCParticles by their number of hits.
     *
     *  @param  selectedMCParticleToGoodHitsMaps the input vector of mappings from selected reconstructable MCParticles to their good hits
     *  @param  orderedMCParticleVector the output vector of ordered MCParticles
     */
    static void GetOrderedMCParticleVector(const LArMCParticleHelper::MCContributionMapVector &selectedMCParticleToGoodHitsMaps, pandora::MCParticleVector &orderedMCParticleVector);

    /**
     *  @brief  Order input Pfos by their number of hits.
     *
     *  @param  pfoToReconstructable2DHitsMap the input vector of mappings from Pfos to their reconstructable hits
     *  @param  orderedPfoVector the output vector of ordered Pfos
     */
    static void GetOrderedPfoVector(const LArMCParticleHelper::PfoContributionMap &pfoToReconstructable2DHitsMap, pandora::PfoVector &orderedPfoVector);

    /**
     *  @brief  Print details of selected MCParticles to the terminal in a table.
     *
     *  @param  selectedMCParticleToGoodHitsMap the input mapping from selected reconstructable MCParticles to their good hits
     *  @param  orderedMCParticleVector the input vector of ordered MCParticles
     */
    static void PrintMCParticleTable(const LArMCParticleHelper::MCContributionMap &selectedMCParticleToGoodHitsMaps, const pandora::MCParticleVector &orderedMCParticleVector);

    /**
     *  @brief  Print details of input Pfos to the terminal in a table.
     *
     *  @param  pfoToReconstructable2DHitsMap the input vector of mappings from Pfos to their reconstructable hits
     *  @param  orderedPfoVector the input vector of ordered Pfos
     */
    static void PrintPfoTable(const LArMCParticleHelper::PfoContributionMap &pfoToReconstructable2DHitsMap, const pandora::PfoVector &orderedPfoVector);

    /**
     *  @brief  Print the shared good hits between all Pfos and MCParticles
     *
     *  @param  orderedPfoVector the input vector of ordered Pfos
     *  @param  orderedMCParticleVector the input vector of ordered MCParticles
     *  @param  mcParticleToPfoHitSharingMap the output mapping from selected reconstructable MCParticles to Pfos and the number hits shared
     *  @param  nMatches the maximum number of Pfo matches to show
     */
    static void PrintMatchingTable(const pandora::PfoVector &orderedPfoVector, const pandora::MCParticleVector &orderedMCParticleVector, const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcParticleToPfoHitSharingMap, const unsigned int nMatches);
};

} // namespace lar_content

#endif // #ifndef LAR_MONITORING_HELPER_H
